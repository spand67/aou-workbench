"""Held-out clinical and PRS model comparison utilities."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from .clinical_model import (
    _average_precision,
    _roc_auc,
    clinical_model_predictions_path,
    run_clinical_model,
)
from .config import ProjectConfig
from .io_utils import read_table, slugify, write_dataframe, write_text
from .microarray_plink_prs import microarray_prs_scores_path
from .paths import ProjectPaths, join_path
from .statistics import run_binary_logistic_regression


def model_comparison_dir(paths: ProjectPaths, label: str) -> str:
    return join_path(paths.run_root, "clinical", "model_comparison", slugify(label))


def model_comparison_metrics_path(paths: ProjectPaths, label: str) -> str:
    return join_path(model_comparison_dir(paths, label), "metrics.tsv")


def model_comparison_predictions_path(paths: ProjectPaths, label: str) -> str:
    return join_path(model_comparison_dir(paths, label), "predictions.tsv")


def model_comparison_report_path(paths: ProjectPaths, label: str) -> str:
    return join_path(model_comparison_dir(paths, label), "report.md")


def _default_label(gwas_label: str, prs_label: str) -> str:
    return f"{gwas_label}_{prs_label}_heldout"


def _metric_row(
    *,
    model: str,
    threshold_label: str,
    y: np.ndarray,
    score: np.ndarray,
    score_is_probability: bool,
) -> dict[str, object]:
    row: dict[str, object] = {
        "model": model,
        "threshold_label": threshold_label,
        "n_participants": int(len(y)),
        "n_cases": int(y.sum()),
        "n_controls": int(len(y) - y.sum()),
        "roc_auc": _roc_auc(y, score),
        "average_precision": _average_precision(y, score),
        "brier_score": float("nan"),
        "log_loss": float("nan"),
        "or_per_score_sd": float("nan"),
        "beta_per_score_sd": float("nan"),
        "se_per_score_sd": float("nan"),
        "p_value": float("nan"),
    }
    if score_is_probability:
        clipped = np.clip(score, 1e-12, 1 - 1e-12)
        row["brier_score"] = float(np.mean((score - y) ** 2)) if len(y) else float("nan")
        row["log_loss"] = float(-np.mean(y * np.log(clipped) + (1 - y) * np.log(1 - clipped))) if len(y) else float("nan")
    return row


def _prs_regression_row(frame: pd.DataFrame, threshold_label: str) -> dict[str, object]:
    raw = pd.to_numeric(frame["prs_score"], errors="coerce")
    sd = float(raw.std(ddof=0))
    if math.isfinite(sd) and sd > 0:
        standardized = (raw - float(raw.mean())) / sd
    else:
        standardized = raw * np.nan
    exposure = pd.Series(standardized.to_numpy(), index=frame["person_id"].astype(str))
    regression = run_binary_logistic_regression(
        frame[["person_id", "analysis_case"]],
        exposure,
        outcome_column="analysis_case",
        covariates=(),
    )
    y = frame["analysis_case"].astype(int).to_numpy()
    row = _metric_row(
        model="prs_only",
        threshold_label=threshold_label,
        y=y,
        score=raw.astype(float).to_numpy(),
        score_is_probability=False,
    )
    row.update(
        {
            "or_per_score_sd": regression["odds_ratio"],
            "beta_per_score_sd": regression["beta"],
            "se_per_score_sd": regression["se"],
            "p_value": regression["regression_p"],
        }
    )
    return row


def _write_report(
    *,
    metrics: pd.DataFrame,
    gwas_label: str,
    prs_label: str,
    path: str,
) -> None:
    try:
        metrics_md = metrics.to_markdown(index=False)
    except ImportError:
        metrics_md = metrics.to_string(index=False)
    body = [
        "# Held-Out Clinical And PRS Model Comparison",
        "",
        f"- GWAS label: `{gwas_label}`",
        f"- PRS label: `{prs_label}`",
        "- Evaluation set: held-out test participants with both clinical predictions and PRS scores.",
        "- Clinical-only predictions are trained on the training split by `run-clinical-model`.",
        "- PRS-only scores use train-split GWAS weights and held-out test scoring from `run-microarray-plink-prs`.",
        "- No clinical+PRS combiner is trained here; that requires cross-fitted or validation PRS scores to avoid test-set tuning.",
        "",
        "## Metrics",
        "",
        metrics_md,
        "",
    ]
    write_text("\n".join(body), path)


def run_heldout_model_comparison(
    config: ProjectConfig,
    paths: ProjectPaths,
    *,
    gwas_label: str,
    prs_label: str,
    label: str | None = None,
    rerun_clinical_model: bool = False,
    eligibility_flag: str = "primary_model_eligible",
    l2_penalty: float = 1.0,
) -> dict[str, pd.DataFrame | str]:
    output_label = label or _default_label(gwas_label, prs_label)
    clinical_path = Path(clinical_model_predictions_path(paths))
    if rerun_clinical_model or not clinical_path.exists():
        run_clinical_model(config, paths, eligibility_flag=eligibility_flag, l2_penalty=l2_penalty)
    if not clinical_path.exists():
        raise RuntimeError("Clinical predictions are missing. Run `aou-workbench run-clinical-model` first.")

    prs_path = Path(microarray_prs_scores_path(paths, gwas_label, prs_label))
    if not prs_path.exists():
        raise RuntimeError(
            f"PRS scores are missing for GWAS label `{gwas_label}` and PRS label `{prs_label}`. "
            "Run `aou-workbench run-microarray-plink-prs` first."
        )

    clinical = read_table(str(clinical_path))
    prs = read_table(str(prs_path))
    if clinical.empty or prs.empty:
        raise RuntimeError("Clinical predictions and PRS scores must both be non-empty.")

    outcome = config.analysis.matched_outcome_column
    clinical = clinical[clinical["analysis_split"].astype(str) == "test"].copy()
    clinical["person_id"] = clinical["person_id"].astype(str)
    clinical["analysis_case"] = pd.to_numeric(clinical[outcome], errors="coerce").astype("Int64")
    clinical["clinical_probability"] = pd.to_numeric(clinical["predicted_probability"], errors="coerce")

    prs = prs.copy()
    prs["person_id"] = prs["person_id"].astype(str)
    prs["analysis_case"] = pd.to_numeric(prs["analysis_case"], errors="coerce").astype("Int64")
    prs["prs_score"] = pd.to_numeric(prs["prs_score"], errors="coerce")

    merged = prs.merge(
        clinical[["person_id", "analysis_case", "clinical_probability"]],
        on=["person_id", "analysis_case"],
        how="inner",
    )
    merged = merged.dropna(subset=["analysis_case", "clinical_probability", "prs_score"]).copy()
    if merged.empty:
        raise RuntimeError("No held-out participants have both clinical predictions and PRS scores.")
    merged["analysis_case"] = merged["analysis_case"].astype(int)

    rows: list[dict[str, object]] = []
    for threshold_label, group in merged.groupby("threshold_label", sort=False):
        y = group["analysis_case"].astype(int).to_numpy()
        rows.append(
            _metric_row(
                model="clinical_only",
                threshold_label=str(threshold_label),
                y=y,
                score=group["clinical_probability"].astype(float).to_numpy(),
                score_is_probability=True,
            )
        )
        rows.append(_prs_regression_row(group, str(threshold_label)))
    metrics = pd.DataFrame(rows)

    prediction_columns = [
        "person_id",
        "analysis_case",
        "threshold_label",
        "clinical_probability",
        "prs_score",
    ]
    predictions = merged[prediction_columns].sort_values(["threshold_label", "person_id"]).reset_index(drop=True)
    write_dataframe(metrics, model_comparison_metrics_path(paths, output_label))
    write_dataframe(predictions, model_comparison_predictions_path(paths, output_label))
    _write_report(
        metrics=metrics,
        gwas_label=gwas_label,
        prs_label=prs_label,
        path=model_comparison_report_path(paths, output_label),
    )
    return {
        "metrics": metrics,
        "predictions": predictions,
        "metrics_path": model_comparison_metrics_path(paths, output_label),
        "predictions_path": model_comparison_predictions_path(paths, output_label),
        "report_path": model_comparison_report_path(paths, output_label),
    }


__all__ = [
    "model_comparison_metrics_path",
    "model_comparison_predictions_path",
    "model_comparison_report_path",
    "run_heldout_model_comparison",
]
