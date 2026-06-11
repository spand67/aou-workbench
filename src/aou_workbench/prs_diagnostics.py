"""Diagnostics for existing PRS and clinical+PRS outputs."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import fisher_exact

from .clinical_model import _average_precision, _roc_auc
from .clinical_prs_model import clinical_prs_model_predictions_path
from .cohort_summary import clinical_model_input_path
from .config import ProjectConfig
from .io_utils import read_table, slugify, write_dataframe, write_json, write_text
from .microarray_plink_prs import microarray_prs_scores_path
from .paths import ProjectPaths, join_path
from .statistics import run_binary_logistic_regression


def prs_diagnostics_dir(paths: ProjectPaths, label: str) -> str:
    return join_path(paths.run_root, "clinical", "prs_diagnostics", slugify(label))


def prs_diagnostics_overall_metrics_path(paths: ProjectPaths, label: str) -> str:
    return join_path(prs_diagnostics_dir(paths, label), "overall_metrics.tsv")


def prs_diagnostics_bootstrap_ci_path(paths: ProjectPaths, label: str) -> str:
    return join_path(prs_diagnostics_dir(paths, label), "bootstrap_ci.tsv")


def prs_diagnostics_deciles_path(paths: ProjectPaths, label: str) -> str:
    return join_path(prs_diagnostics_dir(paths, label), "prs_deciles.tsv")


def prs_diagnostics_ancestry_path(paths: ProjectPaths, label: str) -> str:
    return join_path(prs_diagnostics_dir(paths, label), "ancestry_stratified_metrics.tsv")


def prs_diagnostics_definite_path(paths: ProjectPaths, label: str) -> str:
    return join_path(prs_diagnostics_dir(paths, label), "definite_case_sensitivity.tsv")


def prs_diagnostics_cofactor_path(paths: ProjectPaths, label: str) -> str:
    return join_path(prs_diagnostics_dir(paths, label), "cofactor_stratified_metrics.tsv")


def prs_diagnostics_calibration_path(paths: ProjectPaths, label: str) -> str:
    return join_path(prs_diagnostics_dir(paths, label), "calibration.tsv")


def prs_diagnostics_qc_path(paths: ProjectPaths, label: str) -> str:
    return join_path(prs_diagnostics_dir(paths, label), "diagnostics_qc.json")


def prs_diagnostics_report_path(paths: ProjectPaths, label: str) -> str:
    return join_path(prs_diagnostics_dir(paths, label), "report.md")


def _sigmoid(values: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(values, -40, 40)))


def _logit(values: pd.Series | np.ndarray) -> np.ndarray:
    clipped = np.clip(np.asarray(values, dtype=float), 1e-6, 1 - 1e-6)
    return np.log(clipped / (1 - clipped))


def _metric_values(y: np.ndarray, score: np.ndarray, *, probability: bool) -> dict[str, float]:
    row = {
        "roc_auc": _roc_auc(y, score),
        "average_precision": _average_precision(y, score),
        "brier_score": float("nan"),
        "log_loss": float("nan"),
    }
    if probability:
        clipped = np.clip(score, 1e-12, 1 - 1e-12)
        row["brier_score"] = float(np.mean((score - y) ** 2)) if len(y) else float("nan")
        row["log_loss"] = float(-np.mean(y * np.log(clipped) + (1 - y) * np.log(1 - clipped))) if len(y) else float("nan")
    return row


def _or_per_sd(frame: pd.DataFrame, score_column: str, outcome_column: str) -> dict[str, float]:
    work = frame[[outcome_column, score_column]].copy()
    work["_row_id"] = np.arange(len(work)).astype(str)
    work = work.dropna(subset=[outcome_column, score_column]).copy()
    score = pd.to_numeric(work[score_column], errors="coerce")
    sd = float(score.std(ddof=0))
    if not math.isfinite(sd) or sd == 0:
        return {"beta": float("nan"), "se": float("nan"), "odds_ratio": float("nan"), "regression_p": float("nan"), "n_samples": int(len(work))}
    standardized = (score - float(score.mean())) / sd
    exposure = pd.Series(standardized.to_numpy(), index=work["_row_id"].astype(str))
    sample_df = work.rename(columns={"_row_id": "person_id"})[["person_id", outcome_column]]
    return run_binary_logistic_regression(sample_df, exposure, outcome_column=outcome_column, covariates=())


def _score_metric_row(frame: pd.DataFrame, *, score_name: str, score_column: str, outcome_column: str, probability: bool) -> dict[str, object]:
    used = frame.dropna(subset=[score_column, outcome_column]).copy()
    y = pd.to_numeric(used[outcome_column], errors="coerce").astype(int).to_numpy()
    score = pd.to_numeric(used[score_column], errors="coerce").astype(float).to_numpy()
    regression = _or_per_sd(used, score_column, outcome_column)
    return {
        "score_name": score_name,
        "n": int(len(used)),
        "cases": int(y.sum()),
        "controls": int(len(y) - y.sum()),
        "prevalence": float(y.mean()) if len(y) else float("nan"),
        **_metric_values(y, score, probability=probability),
        "or_per_score_sd": regression["odds_ratio"],
        "beta_per_score_sd": regression["beta"],
        "se_per_score_sd": regression["se"],
        "p_value": regression["regression_p"],
    }


def _fit_recalibration(train: pd.DataFrame, *, probability_column: str, outcome_column: str) -> dict[str, float]:
    used = train.dropna(subset=[probability_column, outcome_column]).copy()
    y = pd.to_numeric(used[outcome_column], errors="coerce").astype(float).to_numpy()
    x_raw = _logit(pd.to_numeric(used[probability_column], errors="coerce"))
    x = np.column_stack([np.ones(len(used), dtype=float), x_raw])
    if len(used) <= 2 or len(np.unique(y)) < 2:
        return {"intercept": 0.0, "slope": 1.0, "n": int(len(used))}

    def objective(beta: np.ndarray) -> tuple[float, np.ndarray]:
        eta = x @ beta
        loss = float(np.sum(np.logaddexp(0, eta) - y * eta))
        p = _sigmoid(eta)
        return loss, x.T @ (p - y)

    opt = minimize(
        lambda beta: objective(beta)[0],
        np.array([0.0, 1.0], dtype=float),
        jac=lambda beta: objective(beta)[1],
        method="BFGS",
        options={"maxiter": 1000},
    )
    if not opt.success:
        return {"intercept": 0.0, "slope": 1.0, "n": int(len(used))}
    return {"intercept": float(opt.x[0]), "slope": float(opt.x[1]), "n": int(len(used))}


def _apply_recalibration(frame: pd.DataFrame, *, probability_column: str, intercept: float, slope: float) -> pd.Series:
    return pd.Series(_sigmoid(intercept + slope * _logit(pd.to_numeric(frame[probability_column], errors="coerce"))), index=frame.index)


def _load_analysis_frame(
    config: ProjectConfig,
    paths: ProjectPaths,
    *,
    gwas_label: str,
    prs_label: str,
    clinical_prs_label: str | None,
    threshold_label: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    model_input = read_table(clinical_model_input_path(paths))
    outcome = config.analysis.matched_outcome_column
    model_input = model_input.copy()
    model_input["person_id"] = model_input["person_id"].astype(str)
    model_input["analysis_case"] = pd.to_numeric(model_input[outcome], errors="coerce")

    prs_path = Path(microarray_prs_scores_path(paths, gwas_label, prs_label))
    if not prs_path.exists():
        raise RuntimeError(f"PRS scores are missing: {prs_path}")
    prs = read_table(str(prs_path))
    prs = prs[prs["threshold_label"].astype(str) == threshold_label].copy()
    prs["person_id"] = prs["person_id"].astype(str)
    prs["prs_score"] = pd.to_numeric(prs["prs_score"], errors="coerce")
    merged = model_input.merge(prs[["person_id", "threshold_label", "prs_score"]], on="person_id", how="left")

    qc: dict[str, object] = {
        "prs_scores_path": str(prs_path),
        "clinical_model_input_rows": int(model_input.shape[0]),
        "prs_threshold_rows": int(prs.shape[0]),
        "merged_rows": int(merged.shape[0]),
        "merged_prs_nonmissing_rows": int(merged["prs_score"].notna().sum()),
    }

    if clinical_prs_label:
        clinical_prs_path = Path(clinical_prs_model_predictions_path(paths, clinical_prs_label))
        qc["clinical_prs_predictions_path"] = str(clinical_prs_path)
        if clinical_prs_path.exists():
            pred = read_table(str(clinical_prs_path))
            pred = pred[pred["threshold_label"].astype(str) == threshold_label].copy() if "threshold_label" in pred.columns else pred.copy()
            pred["person_id"] = pred["person_id"].astype(str)
            pred["clinical_prs_probability"] = pd.to_numeric(pred["predicted_probability"], errors="coerce")
            merged = merged.merge(
                pred[["person_id", "clinical_prs_probability"]],
                on="person_id",
                how="left",
            )
            qc["clinical_prs_prediction_rows"] = int(pred.shape[0])
        else:
            qc["clinical_prs_prediction_rows"] = 0
    return merged, qc


def _overall_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    test = frame[(frame["analysis_split"].astype(str) == "test") & frame["analysis_case"].isin([0, 1])].copy()
    rows = [_score_metric_row(test, score_name="prs_only", score_column="prs_score", outcome_column="analysis_case", probability=False)]
    if "clinical_prs_probability" in test.columns and test["clinical_prs_probability"].notna().any():
        rows.append(
            _score_metric_row(
                test,
                score_name="clinical_prs_raw",
                score_column="clinical_prs_probability",
                outcome_column="analysis_case",
                probability=True,
            )
        )
        train = frame[(frame["analysis_split"].astype(str) == "train") & frame["analysis_case"].isin([0, 1])].copy()
        recal = _fit_recalibration(train, probability_column="clinical_prs_probability", outcome_column="analysis_case")
        test["clinical_prs_recalibrated_probability"] = _apply_recalibration(
            test,
            probability_column="clinical_prs_probability",
            intercept=recal["intercept"],
            slope=recal["slope"],
        )
        recal_row = _score_metric_row(
            test,
            score_name="clinical_prs_recalibrated",
            score_column="clinical_prs_recalibrated_probability",
            outcome_column="analysis_case",
            probability=True,
        )
        recal_row["recalibration_intercept"] = recal["intercept"]
        recal_row["recalibration_slope"] = recal["slope"]
        recal_row["recalibration_train_n"] = recal["n"]
        rows.append(recal_row)
    return pd.DataFrame(rows)


def _decile_table(frame: pd.DataFrame) -> pd.DataFrame:
    test = frame[(frame["analysis_split"].astype(str) == "test") & frame["analysis_case"].isin([0, 1])].dropna(subset=["prs_score"]).copy()
    if test.empty:
        return pd.DataFrame()
    ranks = test["prs_score"].rank(method="first")
    test["prs_decile"] = pd.qcut(ranks, q=min(10, len(test)), labels=False, duplicates="drop") + 1
    rows = []
    ref = test[test["prs_decile"] == int(test["prs_decile"].min())]
    ref_cases = int(ref["analysis_case"].eq(1).sum())
    ref_controls = int(ref["analysis_case"].eq(0).sum())
    for decile, group in test.groupby("prs_decile", sort=True):
        cases = int(group["analysis_case"].eq(1).sum())
        controls = int(group["analysis_case"].eq(0).sum())
        a, b, c, d = cases + 0.5, controls + 0.5, ref_cases + 0.5, ref_controls + 0.5
        odds_ratio = (a * d) / (b * c)
        se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        low = math.exp(math.log(odds_ratio) - 1.96 * se)
        high = math.exp(math.log(odds_ratio) + 1.96 * se)
        fisher_p = fisher_exact([[cases, controls], [ref_cases, ref_controls]])[1] if cases + controls and ref_cases + ref_controls else float("nan")
        rows.append(
            {
                "prs_decile": int(decile),
                "n": int(len(group)),
                "cases": cases,
                "controls": controls,
                "prevalence": cases / len(group) if len(group) else float("nan"),
                "prs_min": float(group["prs_score"].min()),
                "prs_median": float(group["prs_score"].median()),
                "prs_max": float(group["prs_score"].max()),
                "odds_ratio_vs_lowest_decile": float(odds_ratio),
                "or_ci_low": float(low),
                "or_ci_high": float(high),
                "fisher_p_vs_lowest_decile": float(fisher_p),
            }
        )
    return pd.DataFrame(rows)


def _stratified_metrics(frame: pd.DataFrame, *, group_column: str, score_columns: dict[str, tuple[str, bool]]) -> pd.DataFrame:
    test = frame[(frame["analysis_split"].astype(str) == "test") & frame["analysis_case"].isin([0, 1])].copy()
    if group_column not in test.columns:
        return pd.DataFrame()
    rows = []
    test[group_column] = test[group_column].astype("string").fillna("missing").astype(str)
    for group_value, group in test.groupby(group_column, sort=True):
        for score_name, (score_column, probability) in score_columns.items():
            if score_column not in group.columns or group[score_column].notna().sum() == 0:
                continue
            metric = _score_metric_row(group, score_name=score_name, score_column=score_column, outcome_column="analysis_case", probability=probability)
            rows.append({"stratum_variable": group_column, "stratum": group_value, **metric})
    return pd.DataFrame(rows)


def _definite_sensitivity(frame: pd.DataFrame, score_columns: dict[str, tuple[str, bool]]) -> pd.DataFrame:
    if "definite_rhabdo_case" not in frame.columns or "case_tier" not in frame.columns:
        return pd.DataFrame()
    test = frame[frame["analysis_split"].astype(str) == "test"].copy()
    tier = test["case_tier"].astype("string").str.lower()
    test = test[(tier.eq("control")) | (pd.to_numeric(test["definite_rhabdo_case"], errors="coerce").fillna(0).eq(1))].copy()
    test["definite_case"] = pd.to_numeric(test["definite_rhabdo_case"], errors="coerce").fillna(0).astype(int)
    rows = []
    for score_name, (score_column, probability) in score_columns.items():
        if score_column not in test.columns or test[score_column].notna().sum() == 0:
            continue
        metric = _score_metric_row(test, score_name=score_name, score_column=score_column, outcome_column="definite_case", probability=probability)
        rows.append(metric)
    return pd.DataFrame(rows)


def _cofactor_strata(frame: pd.DataFrame, score_columns: dict[str, tuple[str, bool]]) -> pd.DataFrame:
    candidates = [
        "preindex_sepsis",
        "periindex_sepsis",
        "postindex_sepsis",
        "remote_preindex_sepsis",
        "near_preindex_sepsis",
        "preindex_renal_injury",
        "periindex_renal_injury",
        "postindex_renal_injury",
        "remote_preindex_renal_injury",
        "near_preindex_renal_injury",
    ]
    rows = []
    for column in candidates:
        if column not in frame.columns:
            continue
        work = frame.copy()
        work[column] = (pd.to_numeric(work[column], errors="coerce").fillna(0) >= 1).astype(int).astype(str)
        metrics = _stratified_metrics(work, group_column=column, score_columns=score_columns)
        if not metrics.empty:
            rows.append(metrics)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _bootstrap_ci(
    frame: pd.DataFrame,
    *,
    score_columns: dict[str, tuple[str, bool]],
    iterations: int,
    seed: int,
) -> pd.DataFrame:
    test = frame[(frame["analysis_split"].astype(str) == "test") & frame["analysis_case"].isin([0, 1])].copy()
    if test.empty or iterations <= 0:
        return pd.DataFrame()
    rng = np.random.default_rng(seed)
    rows = []
    for score_name, (score_column, probability) in score_columns.items():
        if score_column not in test.columns:
            continue
        used = test.dropna(subset=[score_column, "analysis_case"]).reset_index(drop=True)
        if used.empty:
            continue
        values: dict[str, list[float]] = {"roc_auc": [], "average_precision": [], "or_per_score_sd": []}
        for _ in range(iterations):
            sample = used.iloc[rng.integers(0, len(used), len(used))].copy()
            if sample["analysis_case"].nunique() < 2:
                continue
            y = sample["analysis_case"].astype(int).to_numpy()
            score = pd.to_numeric(sample[score_column], errors="coerce").astype(float).to_numpy()
            values["roc_auc"].append(_roc_auc(y, score))
            values["average_precision"].append(_average_precision(y, score))
            values["or_per_score_sd"].append(_or_per_sd(sample, score_column, "analysis_case")["odds_ratio"])
        point = _score_metric_row(used, score_name=score_name, score_column=score_column, outcome_column="analysis_case", probability=probability)
        for metric, metric_values in values.items():
            clean = np.asarray([value for value in metric_values if pd.notna(value) and np.isfinite(value)], dtype=float)
            rows.append(
                {
                    "score_name": score_name,
                    "metric": metric,
                    "point_estimate": point[metric] if metric in point else point["or_per_score_sd"],
                    "ci_low": float(np.quantile(clean, 0.025)) if len(clean) else float("nan"),
                    "ci_high": float(np.quantile(clean, 0.975)) if len(clean) else float("nan"),
                    "bootstrap_iterations_requested": int(iterations),
                    "bootstrap_iterations_used": int(len(clean)),
                }
            )
    return pd.DataFrame(rows)


def _calibration_table(frame: pd.DataFrame) -> pd.DataFrame:
    test = frame[(frame["analysis_split"].astype(str) == "test") & frame["analysis_case"].isin([0, 1])].copy()
    if "clinical_prs_probability" not in test.columns or test["clinical_prs_probability"].notna().sum() == 0:
        return pd.DataFrame()
    train = frame[(frame["analysis_split"].astype(str) == "train") & frame["analysis_case"].isin([0, 1])].copy()
    recal = _fit_recalibration(train, probability_column="clinical_prs_probability", outcome_column="analysis_case")
    test["clinical_prs_recalibrated_probability"] = _apply_recalibration(
        test,
        probability_column="clinical_prs_probability",
        intercept=recal["intercept"],
        slope=recal["slope"],
    )
    rows = []
    for model, column in (
        ("clinical_prs_raw", "clinical_prs_probability"),
        ("clinical_prs_recalibrated", "clinical_prs_recalibrated_probability"),
    ):
        used = test.dropna(subset=[column]).copy()
        ranks = used[column].rank(method="first")
        used["calibration_bin"] = pd.qcut(ranks, q=min(10, len(used)), labels=False, duplicates="drop") + 1
        for bin_id, group in used.groupby("calibration_bin", sort=True):
            rows.append(
                {
                    "model": model,
                    "bin": int(bin_id),
                    "n": int(len(group)),
                    "mean_predicted_probability": float(group[column].mean()),
                    "observed_case_fraction": float(group["analysis_case"].mean()),
                    "recalibration_intercept": recal["intercept"] if model == "clinical_prs_recalibrated" else float("nan"),
                    "recalibration_slope": recal["slope"] if model == "clinical_prs_recalibrated" else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def _write_report(
    *,
    overall: pd.DataFrame,
    deciles: pd.DataFrame,
    ancestry: pd.DataFrame,
    definite: pd.DataFrame,
    cofactor: pd.DataFrame,
    bootstrap: pd.DataFrame,
    calibration: pd.DataFrame,
    path: str,
) -> None:
    def md(frame: pd.DataFrame, limit: int | None = None) -> str:
        if frame.empty:
            return "_No rows available._"
        view = frame.head(limit) if limit else frame
        try:
            return view.to_markdown(index=False)
        except ImportError:
            return view.to_string(index=False)

    body = [
        "# PRS Diagnostics",
        "",
        "Diagnostics generated from existing clinical, PRS, and clinical+PRS outputs. No GWAS or PRS scoring is rerun.",
        "",
        "## Overall Test Metrics",
        "",
        md(overall),
        "",
        "## Bootstrap Confidence Intervals",
        "",
        md(bootstrap),
        "",
        "## PRS Deciles",
        "",
        md(deciles),
        "",
        "## Ancestry-Stratified Metrics",
        "",
        md(ancestry),
        "",
        "## Definite-Case Sensitivity",
        "",
        md(definite),
        "",
        "## Sepsis And Renal-Injury Strata",
        "",
        md(cofactor, limit=80),
        "",
        "## Calibration",
        "",
        md(calibration),
        "",
    ]
    write_text("\n".join(body), path)


def run_prs_diagnostics(
    config: ProjectConfig,
    paths: ProjectPaths,
    *,
    gwas_label: str,
    prs_label: str,
    clinical_prs_label: str | None = None,
    threshold_label: str = "p0_01",
    label: str | None = None,
    bootstrap_iterations: int = 200,
    seed: int = 20260611,
) -> dict[str, pd.DataFrame | str]:
    output_label = label or f"prs_diagnostics_{prs_label}"
    frame, qc = _load_analysis_frame(
        config,
        paths,
        gwas_label=gwas_label,
        prs_label=prs_label,
        clinical_prs_label=clinical_prs_label,
        threshold_label=threshold_label,
    )
    score_columns: dict[str, tuple[str, bool]] = {"prs_only": ("prs_score", False)}
    if "clinical_prs_probability" in frame.columns and frame["clinical_prs_probability"].notna().any():
        train = frame[(frame["analysis_split"].astype(str) == "train") & frame["analysis_case"].isin([0, 1])].copy()
        recal = _fit_recalibration(train, probability_column="clinical_prs_probability", outcome_column="analysis_case")
        frame["clinical_prs_recalibrated_probability"] = _apply_recalibration(
            frame,
            probability_column="clinical_prs_probability",
            intercept=recal["intercept"],
            slope=recal["slope"],
        )
        score_columns["clinical_prs_raw"] = ("clinical_prs_probability", True)
        score_columns["clinical_prs_recalibrated"] = ("clinical_prs_recalibrated_probability", True)
        qc["clinical_prs_recalibration"] = recal

    overall = _overall_metrics(frame)
    deciles = _decile_table(frame)
    ancestry = _stratified_metrics(frame, group_column="model_ancestry_pred", score_columns=score_columns)
    definite = _definite_sensitivity(frame, score_columns)
    cofactor = _cofactor_strata(frame, score_columns)
    bootstrap = _bootstrap_ci(frame, score_columns=score_columns, iterations=bootstrap_iterations, seed=seed)
    calibration = _calibration_table(frame)

    write_dataframe(overall, prs_diagnostics_overall_metrics_path(paths, output_label))
    write_dataframe(deciles, prs_diagnostics_deciles_path(paths, output_label))
    write_dataframe(ancestry, prs_diagnostics_ancestry_path(paths, output_label))
    write_dataframe(definite, prs_diagnostics_definite_path(paths, output_label))
    write_dataframe(cofactor, prs_diagnostics_cofactor_path(paths, output_label))
    write_dataframe(bootstrap, prs_diagnostics_bootstrap_ci_path(paths, output_label))
    write_dataframe(calibration, prs_diagnostics_calibration_path(paths, output_label))
    write_json(
        {
            **qc,
            "gwas_label": gwas_label,
            "prs_label": prs_label,
            "clinical_prs_label": clinical_prs_label,
            "threshold_label": threshold_label,
            "bootstrap_iterations": int(bootstrap_iterations),
            "seed": int(seed),
            "score_columns": list(score_columns),
        },
        prs_diagnostics_qc_path(paths, output_label),
    )
    _write_report(
        overall=overall,
        deciles=deciles,
        ancestry=ancestry,
        definite=definite,
        cofactor=cofactor,
        bootstrap=bootstrap,
        calibration=calibration,
        path=prs_diagnostics_report_path(paths, output_label),
    )
    return {
        "overall": overall,
        "deciles": deciles,
        "ancestry": ancestry,
        "definite": definite,
        "cofactor": cofactor,
        "bootstrap": bootstrap,
        "calibration": calibration,
        "overall_path": prs_diagnostics_overall_metrics_path(paths, output_label),
        "deciles_path": prs_diagnostics_deciles_path(paths, output_label),
        "ancestry_path": prs_diagnostics_ancestry_path(paths, output_label),
        "definite_path": prs_diagnostics_definite_path(paths, output_label),
        "cofactor_path": prs_diagnostics_cofactor_path(paths, output_label),
        "bootstrap_path": prs_diagnostics_bootstrap_ci_path(paths, output_label),
        "calibration_path": prs_diagnostics_calibration_path(paths, output_label),
        "qc_path": prs_diagnostics_qc_path(paths, output_label),
        "report_path": prs_diagnostics_report_path(paths, output_label),
    }


__all__ = [
    "prs_diagnostics_ancestry_path",
    "prs_diagnostics_bootstrap_ci_path",
    "prs_diagnostics_calibration_path",
    "prs_diagnostics_cofactor_path",
    "prs_diagnostics_deciles_path",
    "prs_diagnostics_definite_path",
    "prs_diagnostics_overall_metrics_path",
    "prs_diagnostics_qc_path",
    "prs_diagnostics_report_path",
    "run_prs_diagnostics",
]
