"""Pragmatic clinical plus PRS model for held-out evaluation."""

from __future__ import annotations

from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from .clinical_model import (
    ClinicalModel,
    FeatureSpec,
    _calibration_table,
    _choose_threshold,
    _coefficient_table,
    _curve_points,
    _feature_matrix,
    _fit_feature_specs,
    _fit_logistic,
    _metrics,
    _predict,
    _write_calibration_svg,
    _write_line_svg,
)
from .cohort_summary import clinical_model_input_path
from .config import ProjectConfig
from .io_utils import ensure_parent_dir, read_table, slugify, write_dataframe, write_json, write_text
from .microarray_plink_gwas import (
    _expand_local_prefix,
    _plink_common_options,
    _plink_files_exist,
    _run_command,
    read_plink_fam,
)
from .microarray_plink_prs import (
    _score_column,
    microarray_prs_scores_path,
    microarray_prs_weights_path,
)
from .paths import ProjectPaths, join_path
from .stage4_hail_gwas import _pilot_case_control_definition_mask


def clinical_prs_model_default_label(prs_label: str) -> str:
    return f"clinical_prs_{prs_label}"


def clinical_prs_model_dir(paths: ProjectPaths, label: str) -> str:
    return join_path(paths.run_root, "clinical", "clinical_prs_model", slugify(label))


def clinical_prs_model_metrics_path(paths: ProjectPaths, label: str) -> str:
    return join_path(clinical_prs_model_dir(paths, label), "metrics.tsv")


def clinical_prs_model_coefficients_path(paths: ProjectPaths, label: str) -> str:
    return join_path(clinical_prs_model_dir(paths, label), "coefficients.tsv")


def clinical_prs_model_predictions_path(paths: ProjectPaths, label: str) -> str:
    return join_path(clinical_prs_model_dir(paths, label), "predictions.tsv")


def clinical_prs_model_calibration_path(paths: ProjectPaths, label: str) -> str:
    return join_path(clinical_prs_model_dir(paths, label), "calibration.tsv")


def clinical_prs_model_report_path(paths: ProjectPaths, label: str) -> str:
    return join_path(clinical_prs_model_dir(paths, label), "report.md")


def clinical_prs_model_roc_svg_path(paths: ProjectPaths, label: str) -> str:
    return join_path(clinical_prs_model_dir(paths, label), "roc.svg")


def clinical_prs_model_pr_svg_path(paths: ProjectPaths, label: str) -> str:
    return join_path(clinical_prs_model_dir(paths, label), "precision_recall.svg")


def clinical_prs_model_calibration_svg_path(paths: ProjectPaths, label: str) -> str:
    return join_path(clinical_prs_model_dir(paths, label), "calibration.svg")


def clinical_prs_model_train_keep_path(paths: ProjectPaths, label: str) -> str:
    return join_path(clinical_prs_model_dir(paths, label), "train_keep.tsv")


def clinical_prs_model_train_weights_path(paths: ProjectPaths, label: str) -> str:
    return join_path(clinical_prs_model_dir(paths, label), "train_prs_weights.tsv")


def clinical_prs_model_train_extract_path(paths: ProjectPaths, label: str) -> str:
    return join_path(clinical_prs_model_dir(paths, label), "train_prs_variant_ids.txt")


def _train_score_prefix(paths: ProjectPaths, label: str) -> str:
    return join_path(clinical_prs_model_dir(paths, label), "train_prs")


def _train_score_file(paths: ProjectPaths, label: str) -> Path | None:
    prefix = Path(_train_score_prefix(paths, label))
    matches = sorted(prefix.parent.glob(f"{prefix.name}.sscore"))
    return matches[0] if matches else None


def _write_space_delimited(frame: pd.DataFrame, path: str, *, header: bool) -> None:
    ensure_parent_dir(path)
    frame.to_csv(path, sep=" ", index=False, header=header)


def _write_train_keep(
    config: ProjectConfig,
    model_input: pd.DataFrame,
    fam_df: pd.DataFrame,
    paths: ProjectPaths,
    *,
    label: str,
    eligibility_flag: str,
) -> tuple[pd.DataFrame, dict[str, int]]:
    outcome = config.analysis.matched_outcome_column
    required = {"person_id", outcome, "analysis_split", eligibility_flag}
    missing = sorted(required.difference(model_input.columns))
    if missing:
        raise RuntimeError("Clinical model input is missing clinical+PRS columns: " + ", ".join(missing))

    sample = model_input.copy()
    sample["person_id"] = sample["person_id"].astype(str)
    counts = {"clinical_model_input_participants": int(sample["person_id"].nunique())}
    sample = sample[sample["analysis_split"].astype(str) == "train"].copy()
    counts["after_train_split_participants"] = int(sample["person_id"].nunique())
    sample = sample[pd.to_numeric(sample[eligibility_flag], errors="coerce").fillna(0) == 1].copy()
    counts["after_eligibility_participants"] = int(sample["person_id"].nunique())
    sample = sample[_pilot_case_control_definition_mask(sample, outcome)].copy()
    sample[outcome] = pd.to_numeric(sample[outcome], errors="coerce").astype(int)
    counts["after_case_control_definition_participants"] = int(sample["person_id"].nunique())
    counts["after_case_control_definition_cases"] = int(sample[outcome].eq(1).sum())
    counts["after_case_control_definition_controls"] = int(sample[outcome].eq(0).sum())
    fam_ids = set(fam_df["IID"].astype(str))
    sample = sample[sample["person_id"].isin(fam_ids)].copy()
    sample["FID"] = "0"
    sample["IID"] = sample["person_id"].astype(str)
    counts["after_microarray_fam_overlap_participants"] = int(sample["person_id"].nunique())
    counts["after_microarray_fam_overlap_cases"] = int(sample[outcome].eq(1).sum())
    counts["after_microarray_fam_overlap_controls"] = int(sample[outcome].eq(0).sum())
    _write_space_delimited(sample[["FID", "IID"]], clinical_prs_model_train_keep_path(paths, label), header=False)
    return sample, counts


def _load_threshold_test_prs(paths: ProjectPaths, gwas_label: str, prs_label: str, threshold_label: str) -> pd.DataFrame:
    path = Path(microarray_prs_scores_path(paths, gwas_label, prs_label))
    if not path.exists():
        raise RuntimeError(f"PRS scores are missing: {path}")
    scores = read_table(str(path))
    scores = scores[scores["threshold_label"].astype(str) == threshold_label].copy()
    if scores.empty:
        raise RuntimeError(f"No PRS scores found for threshold label `{threshold_label}` in {path}.")
    scores["person_id"] = scores["person_id"].astype(str)
    scores["prs_score"] = pd.to_numeric(scores["prs_score"], errors="coerce")
    return scores[["person_id", "threshold_label", "prs_score"]].dropna(subset=["prs_score"])


def _load_and_write_threshold_weights(
    paths: ProjectPaths,
    *,
    gwas_label: str,
    prs_label: str,
    label: str,
    p_threshold: float,
) -> pd.DataFrame:
    path = Path(microarray_prs_weights_path(paths, gwas_label, prs_label))
    if not path.exists():
        raise RuntimeError(f"PRS weights are missing: {path}")
    weights = read_table(str(path))
    required = {"ID", "A1", "BETA", "P"}
    missing = sorted(required.difference(weights.columns))
    if missing:
        raise RuntimeError("PRS weights are missing columns: " + ", ".join(missing))
    weights = weights.copy()
    weights["P"] = pd.to_numeric(weights["P"], errors="coerce")
    weights["BETA"] = pd.to_numeric(weights["BETA"], errors="coerce")
    weights = weights.dropna(subset=["ID", "A1", "BETA", "P"]).copy()
    weights = weights[weights["P"] <= float(p_threshold)].copy()
    if weights.empty:
        raise RuntimeError(f"No PRS weights remain at p <= {p_threshold}.")
    write_dataframe(weights[["ID", "A1", "BETA", "P"]], clinical_prs_model_train_weights_path(paths, label))
    ensure_parent_dir(clinical_prs_model_train_extract_path(paths, label))
    weights["ID"].astype(str).to_csv(clinical_prs_model_train_extract_path(paths, label), index=False, header=False)
    return weights.reset_index(drop=True)


def _build_train_score_command(
    *,
    plink2_bin: str,
    plink_prefix: str,
    paths: ProjectPaths,
    label: str,
    threads: int | None,
    memory_mb: int | None,
) -> list[str]:
    return [
        plink2_bin,
        "--bfile",
        plink_prefix,
        "--keep",
        clinical_prs_model_train_keep_path(paths, label),
        "--extract",
        clinical_prs_model_train_extract_path(paths, label),
        "--score",
        clinical_prs_model_train_weights_path(paths, label),
        "1",
        "2",
        "3",
        "header-read",
        "cols=+scoresums",
        "list-variants",
        "--out",
        _train_score_prefix(paths, label),
        *_plink_common_options(threads, memory_mb),
    ]


def _parse_train_scores(paths: ProjectPaths, label: str, threshold_label: str) -> pd.DataFrame:
    path = _train_score_file(paths, label)
    if path is None:
        return pd.DataFrame(columns=["person_id", "threshold_label", "prs_score"])
    frame = pd.read_csv(path, sep=r"\s+", dtype=str)
    if frame.empty:
        return pd.DataFrame(columns=["person_id", "threshold_label", "prs_score"])
    iid_column = "IID" if "IID" in frame.columns else "#IID" if "#IID" in frame.columns else None
    if iid_column is None:
        raise RuntimeError(f"Could not find IID column in {path}.")
    score_column = _score_column(frame)
    return pd.DataFrame(
        {
            "person_id": frame[iid_column].astype(str),
            "threshold_label": threshold_label,
            "prs_score": pd.to_numeric(frame[score_column], errors="coerce"),
        }
    ).dropna(subset=["prs_score"])


def _prs_feature_spec(train: pd.DataFrame) -> FeatureSpec:
    values = pd.to_numeric(train["prs_score"], errors="coerce")
    median = float(values.median())
    filled = values.fillna(median)
    sd = float(filled.std(ddof=0))
    if not np.isfinite(sd) or sd == 0:
        raise RuntimeError("Training PRS scores have zero or undefined variance.")
    return FeatureSpec(
        name="prs_score_per_sd",
        source_column="prs_score",
        kind="continuous",
        median=median,
        mean=float(filled.mean()),
        sd=sd,
    )


def _fit_clinical_prs_model(
    train: pd.DataFrame,
    *,
    outcome_column: str,
    l2_penalty: float,
) -> tuple[ClinicalModel, np.ndarray]:
    specs = (*_fit_feature_specs(train), _prs_feature_spec(train))
    x = _feature_matrix(train, specs)
    y = pd.to_numeric(train[outcome_column], errors="coerce").astype(int).to_numpy()
    beta = _fit_logistic(x, y, l2_penalty=l2_penalty)
    train_prob = 1.0 / (1.0 + np.exp(-np.clip(x @ beta, -40, 40)))
    threshold = _choose_threshold(y, train_prob)
    return ClinicalModel(beta=beta, features=specs, threshold=threshold, l2_penalty=l2_penalty), train_prob


def _write_report(
    *,
    metrics: pd.DataFrame,
    coefficients: pd.DataFrame,
    gwas_label: str,
    prs_label: str,
    threshold_label: str,
    p_threshold: float,
    path: str,
) -> None:
    try:
        metrics_md = metrics.to_markdown(index=False)
        coefficient_md = coefficients.to_markdown(index=False)
    except ImportError:
        metrics_md = metrics.to_string(index=False)
        coefficient_md = coefficients.to_string(index=False)
    body = [
        "# Clinical Plus PRS Model",
        "",
        f"- GWAS label: `{gwas_label}`",
        f"- PRS label: `{prs_label}`",
        f"- PRS threshold: `{threshold_label}` (`p <= {p_threshold:g}`)",
        "- Model: L2-regularized logistic regression trained on the train split and evaluated once on the held-out test split.",
        "- Predictors: age, sex category, ancestry category, and standardized PRS.",
        "- Note: this is the pragmatic option 1 model. The training PRS is derived from a GWAS that included the training participants, so the PRS coefficient and train metrics can be optimistic. The held-out test metrics are the main readout.",
        "",
        "## Metrics",
        "",
        metrics_md,
        "",
        "## Coefficients",
        "",
        coefficient_md,
        "",
    ]
    write_text("\n".join(body), path)


def run_clinical_prs_model(
    config: ProjectConfig,
    paths: ProjectPaths,
    *,
    gwas_label: str,
    prs_label: str,
    plink_prefix: str | None = None,
    plink2_bin: str = "plink2",
    threshold_label: str = "p0_01",
    p_threshold: float = 0.01,
    label: str | None = None,
    eligibility_flag: str = "primary_model_eligible",
    l2_penalty: float = 1.0,
    threads: int | None = None,
    memory_mb: int | None = None,
) -> dict[str, pd.DataFrame | str]:
    output_label = label or clinical_prs_model_default_label(prs_label)
    local_prefix = _expand_local_prefix(plink_prefix or str(Path.home() / "plink_microarray" / "arrays"))
    if not _plink_files_exist(local_prefix):
        raise RuntimeError(f"Missing local PLINK BED/BIM/FAM at prefix {local_prefix}.")
    if shutil.which(plink2_bin) is None and not Path(plink2_bin).exists():
        raise RuntimeError(f"Could not find PLINK2 binary `{plink2_bin}` on PATH.")

    model_input = read_table(clinical_model_input_path(paths))
    outcome = config.analysis.matched_outcome_column
    fam_df = read_plink_fam(local_prefix)
    train_sample, train_counts = _write_train_keep(
        config,
        model_input,
        fam_df,
        paths,
        label=output_label,
        eligibility_flag=eligibility_flag,
    )
    weights = _load_and_write_threshold_weights(
        paths,
        gwas_label=gwas_label,
        prs_label=prs_label,
        label=output_label,
        p_threshold=p_threshold,
    )
    if _train_score_file(paths, output_label):
        print(f"Using existing train PRS score output for {output_label}.", flush=True)
    elif not train_sample.empty:
        _run_command(
            _build_train_score_command(
                plink2_bin=plink2_bin,
                plink_prefix=local_prefix,
                paths=paths,
                label=output_label,
                threads=threads,
                memory_mb=memory_mb,
            )
        )
    train_prs = _parse_train_scores(paths, output_label, threshold_label)
    test_prs = _load_threshold_test_prs(paths, gwas_label, prs_label, threshold_label)

    analysis = model_input.copy()
    analysis["person_id"] = analysis["person_id"].astype(str)
    analysis[outcome] = pd.to_numeric(analysis[outcome], errors="coerce")
    analysis = analysis[
        (pd.to_numeric(analysis[eligibility_flag], errors="coerce").fillna(0) == 1)
        & analysis[outcome].isin([0, 1])
        & analysis["analysis_split"].isin(["train", "test"])
    ].copy()
    prs = pd.concat([train_prs, test_prs], ignore_index=True)
    analysis = analysis.merge(prs[["person_id", "threshold_label", "prs_score"]], on="person_id", how="inner")
    analysis["prs_score"] = pd.to_numeric(analysis["prs_score"], errors="coerce")
    analysis = analysis.dropna(subset=["prs_score"]).copy()
    train = analysis[analysis["analysis_split"] == "train"].copy()
    test = analysis[analysis["analysis_split"] == "test"].copy()
    if train.empty or test.empty:
        raise RuntimeError("Clinical+PRS model requires non-empty train and test rows with PRS scores.")
    if train[outcome].nunique() < 2 or test[outcome].nunique() < 2:
        raise RuntimeError("Clinical+PRS model requires both cases and controls in train and test rows.")

    model, train_prob = _fit_clinical_prs_model(train, outcome_column=outcome, l2_penalty=l2_penalty)
    test_prob = _predict(model, test)
    y_train = train[outcome].astype(int).to_numpy()
    y_test = test[outcome].astype(int).to_numpy()
    metrics = pd.DataFrame(
        [
            _metrics(y_train, train_prob, threshold=model.threshold, label="train"),
            _metrics(y_test, test_prob, threshold=model.threshold, label="test"),
        ]
    )
    metrics.insert(1, "threshold_label", threshold_label)
    coefficients = _coefficient_table(model)
    calibration = _calibration_table(y_test, test_prob)
    predictions = pd.concat(
        [
            train.assign(predicted_probability=train_prob, predicted_label=(train_prob >= model.threshold).astype(int)),
            test.assign(predicted_probability=test_prob, predicted_label=(test_prob >= model.threshold).astype(int)),
        ],
        ignore_index=True,
    )
    prediction_columns = [
        column
        for column in (
            "person_id",
            outcome,
            "analysis_split",
            "cv_fold",
            "match_group_id",
            eligibility_flag,
            "threshold_label",
            "prs_score",
            "predicted_probability",
            "predicted_label",
        )
        if column in predictions.columns
    ]
    predictions = predictions[prediction_columns].copy()

    write_dataframe(metrics, clinical_prs_model_metrics_path(paths, output_label))
    write_dataframe(coefficients, clinical_prs_model_coefficients_path(paths, output_label))
    write_dataframe(predictions, clinical_prs_model_predictions_path(paths, output_label))
    write_dataframe(calibration, clinical_prs_model_calibration_path(paths, output_label))
    write_json(
        {
            "gwas_label": gwas_label,
            "prs_label": prs_label,
            "threshold_label": threshold_label,
            "p_threshold": float(p_threshold),
            "eligibility_flag": eligibility_flag,
            "l2_penalty": float(l2_penalty),
            "train_sample_counts": train_counts,
            "n_weights": int(weights.shape[0]),
            "threshold": model.threshold,
            "features": [spec.__dict__ for spec in model.features],
            "method_note": "Pragmatic option 1: train PRS uses train-derived GWAS weights; held-out test metrics are the primary readout.",
        },
        join_path(clinical_prs_model_dir(paths, output_label), "model_spec.json"),
    )
    _write_line_svg(_curve_points(y_test, test_prob, "roc"), clinical_prs_model_roc_svg_path(paths, output_label), title="Clinical + PRS ROC", x_label="False positive rate", y_label="True positive rate", diagonal=True)
    _write_line_svg(_curve_points(y_test, test_prob, "pr"), clinical_prs_model_pr_svg_path(paths, output_label), title="Clinical + PRS Precision-Recall", x_label="Recall", y_label="Precision")
    _write_calibration_svg(calibration, clinical_prs_model_calibration_svg_path(paths, output_label))
    _write_report(
        metrics=metrics,
        coefficients=coefficients,
        gwas_label=gwas_label,
        prs_label=prs_label,
        threshold_label=threshold_label,
        p_threshold=p_threshold,
        path=clinical_prs_model_report_path(paths, output_label),
    )
    return {
        "metrics": metrics,
        "coefficients": coefficients,
        "predictions": predictions,
        "calibration": calibration,
        "metrics_path": clinical_prs_model_metrics_path(paths, output_label),
        "coefficients_path": clinical_prs_model_coefficients_path(paths, output_label),
        "predictions_path": clinical_prs_model_predictions_path(paths, output_label),
        "calibration_path": clinical_prs_model_calibration_path(paths, output_label),
        "report_path": clinical_prs_model_report_path(paths, output_label),
        "roc_svg_path": clinical_prs_model_roc_svg_path(paths, output_label),
        "pr_svg_path": clinical_prs_model_pr_svg_path(paths, output_label),
        "calibration_svg_path": clinical_prs_model_calibration_svg_path(paths, output_label),
    }


__all__ = [
    "clinical_prs_model_calibration_path",
    "clinical_prs_model_calibration_svg_path",
    "clinical_prs_model_coefficients_path",
    "clinical_prs_model_default_label",
    "clinical_prs_model_metrics_path",
    "clinical_prs_model_predictions_path",
    "clinical_prs_model_pr_svg_path",
    "clinical_prs_model_report_path",
    "clinical_prs_model_roc_svg_path",
    "run_clinical_prs_model",
]
