"""Clinical-only prediction model for the matched rhabdomyolysis cohort."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from .cohort_summary import clinical_model_input_path
from .config import ProjectConfig
from .io_utils import read_table, write_dataframe, write_json, write_text
from .paths import ProjectPaths, join_path


_CONTINUOUS_FEATURES: tuple[tuple[str, str, Callable[[pd.Series], pd.Series]], ...] = (
    ("age_at_index", "age_at_index_per_sd", lambda values: values),
    ("observation_days", "log_observation_days_per_sd", np.log1p),
    ("omop_condition_record_dates", "log_condition_record_dates_per_sd", np.log1p),
)
_BINARY_FEATURES: tuple[str, ...] = ()
_CATEGORICAL_FEATURES = {
    "model_sex_category": ("female", "male", "other_or_unknown", "missing"),
    "model_ancestry_pred": ("eur", "afr", "amr", "eas", "sas", "mid", "other", "missing"),
}


@dataclass(frozen=True)
class FeatureSpec:
    name: str
    source_column: str
    kind: str
    reference: str = ""
    median: float = 0.0
    mean: float = 0.0
    sd: float = 1.0
    transform: str = "identity"


@dataclass(frozen=True)
class ClinicalModel:
    beta: np.ndarray
    features: tuple[FeatureSpec, ...]
    threshold: float
    l2_penalty: float


def clinical_model_dir(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "clinical", "model")


def clinical_model_metrics_path(paths: ProjectPaths) -> str:
    return join_path(clinical_model_dir(paths), "metrics.tsv")


def clinical_model_coefficients_path(paths: ProjectPaths) -> str:
    return join_path(clinical_model_dir(paths), "coefficients.tsv")


def clinical_model_predictions_path(paths: ProjectPaths) -> str:
    return join_path(clinical_model_dir(paths), "predictions.tsv")


def clinical_model_cv_metrics_path(paths: ProjectPaths) -> str:
    return join_path(clinical_model_dir(paths), "cv_metrics.tsv")


def clinical_model_calibration_path(paths: ProjectPaths) -> str:
    return join_path(clinical_model_dir(paths), "calibration.tsv")


def clinical_model_report_path(paths: ProjectPaths) -> str:
    return join_path(clinical_model_dir(paths), "report.md")


def clinical_model_roc_svg_path(paths: ProjectPaths) -> str:
    return join_path(clinical_model_dir(paths), "roc.svg")


def clinical_model_pr_svg_path(paths: ProjectPaths) -> str:
    return join_path(clinical_model_dir(paths), "precision_recall.svg")


def clinical_model_calibration_svg_path(paths: ProjectPaths) -> str:
    return join_path(clinical_model_dir(paths), "calibration.svg")


def _sigmoid(values: np.ndarray) -> np.ndarray:
    clipped = np.clip(values, -40, 40)
    return 1.0 / (1.0 + np.exp(-clipped))


def _stable_category(series: pd.Series) -> pd.Series:
    values = series.astype("string").str.lower().str.strip()
    missing = values.isna() | values.isin({"", "nan", "none", "null", "na", "n/a"})
    return values.where(~missing, "missing").astype(str)


def _numeric_values(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def _fit_feature_specs(train_df: pd.DataFrame) -> tuple[FeatureSpec, ...]:
    specs: list[FeatureSpec] = []
    for source_column, feature_name, transform in _CONTINUOUS_FEATURES:
        if source_column not in train_df.columns:
            continue
        raw = _numeric_values(train_df, source_column)
        transformed = transform(raw.clip(lower=0)) if transform is np.log1p else transform(raw)
        if transformed.notna().sum() == 0:
            continue
        median = float(transformed.median())
        filled = transformed.fillna(median)
        sd = float(filled.std(ddof=0))
        if not np.isfinite(sd) or sd == 0:
            continue
        specs.append(
            FeatureSpec(
                name=feature_name,
                source_column=source_column,
                kind="continuous",
                median=median,
                mean=float(filled.mean()),
                sd=sd,
                transform="log1p" if transform is np.log1p else "identity",
            )
        )
    for source_column in _BINARY_FEATURES:
        if source_column not in train_df.columns:
            continue
        values = (_numeric_values(train_df, source_column).fillna(0) >= 1).astype(float)
        if values.nunique(dropna=False) <= 1:
            continue
        specs.append(FeatureSpec(name=source_column, source_column=source_column, kind="binary"))
    for source_column, preferred in _CATEGORICAL_FEATURES.items():
        if source_column not in train_df.columns:
            continue
        values = _stable_category(train_df[source_column])
        observed = [category for category in preferred if category in set(values)]
        observed.extend(sorted(set(values).difference(observed)))
        if len(observed) <= 1:
            continue
        reference = observed[0]
        for category in observed[1:]:
            specs.append(
                FeatureSpec(
                    name=f"{source_column}={category}",
                    source_column=source_column,
                    kind="categorical",
                    reference=reference,
                )
            )
    return tuple(specs)


def _feature_matrix(frame: pd.DataFrame, specs: tuple[FeatureSpec, ...]) -> np.ndarray:
    columns = [np.ones(len(frame), dtype=float)]
    for spec in specs:
        if spec.kind == "continuous":
            values = _numeric_values(frame, spec.source_column)
            if spec.transform == "log1p":
                values = np.log1p(values.clip(lower=0))
            values = values.fillna(spec.median)
            columns.append(((values - spec.mean) / spec.sd).astype(float).to_numpy())
        elif spec.kind == "binary":
            columns.append((_numeric_values(frame, spec.source_column).fillna(0) >= 1).astype(float).to_numpy())
        elif spec.kind == "categorical":
            category = spec.name.split("=", 1)[1]
            values = _stable_category(frame[spec.source_column]) if spec.source_column in frame.columns else pd.Series("missing", index=frame.index)
            columns.append((values == category).astype(float).to_numpy())
        else:  # pragma: no cover - defensive branch
            raise ValueError(f"Unknown feature kind: {spec.kind}")
    return np.column_stack(columns)


def _fit_logistic(x: np.ndarray, y: np.ndarray, *, l2_penalty: float) -> np.ndarray:
    def objective(beta: np.ndarray) -> tuple[float, np.ndarray]:
        eta = x @ beta
        loss = float(np.sum(np.logaddexp(0, eta) - y * eta))
        if l2_penalty:
            loss += 0.5 * l2_penalty * float(np.sum(beta[1:] ** 2))
        p = _sigmoid(eta)
        grad = x.T @ (p - y)
        if l2_penalty:
            grad[1:] += l2_penalty * beta[1:]
        return loss, grad

    opt = minimize(
        lambda beta: objective(beta)[0],
        np.zeros(x.shape[1], dtype=float),
        jac=lambda beta: objective(beta)[1],
        method="BFGS",
        options={"maxiter": 1000},
    )
    if not opt.success:
        opt = minimize(
            lambda beta: objective(beta)[0],
            np.zeros(x.shape[1], dtype=float),
            jac=lambda beta: objective(beta)[1],
            method="L-BFGS-B",
            options={"maxiter": 1000},
        )
    if not opt.success:
        raise RuntimeError(f"Clinical logistic model failed to converge: {opt.message}")
    return np.asarray(opt.x, dtype=float)


def _roc_auc(y_true: np.ndarray, score: np.ndarray) -> float:
    positives = y_true == 1
    n_pos = int(positives.sum())
    n_neg = int((~positives).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = pd.Series(score).rank(method="average").to_numpy()
    rank_sum = float(ranks[positives].sum())
    return (rank_sum - (n_pos * (n_pos + 1) / 2.0)) / (n_pos * n_neg)


def _average_precision(y_true: np.ndarray, score: np.ndarray) -> float:
    order = np.argsort(-score)
    y_sorted = y_true[order]
    positives = int(y_sorted.sum())
    if positives == 0:
        return float("nan")
    precision = np.cumsum(y_sorted) / (np.arange(len(y_sorted)) + 1)
    return float((precision * y_sorted).sum() / positives)


def _roc_curve(y_true: np.ndarray, score: np.ndarray) -> pd.DataFrame:
    thresholds = np.r_[np.inf, np.sort(np.unique(score))[::-1], -np.inf]
    rows = []
    positives = max(int((y_true == 1).sum()), 1)
    negatives = max(int((y_true == 0).sum()), 1)
    for threshold in thresholds:
        pred = score >= threshold
        tp = int(((pred == 1) & (y_true == 1)).sum())
        fp = int(((pred == 1) & (y_true == 0)).sum())
        rows.append({"threshold": threshold, "sensitivity": tp / positives, "fpr": fp / negatives})
    return pd.DataFrame(rows)


def _choose_threshold(y_true: np.ndarray, score: np.ndarray) -> float:
    curve = _roc_curve(y_true, score)
    finite = curve[np.isfinite(curve["threshold"])].copy()
    if finite.empty:
        return 0.5
    finite["youden"] = finite["sensitivity"] - finite["fpr"]
    return float(finite.sort_values(["youden", "threshold"], ascending=[False, False]).iloc[0]["threshold"])


def _metrics(y_true: np.ndarray, score: np.ndarray, *, threshold: float, label: str) -> dict[str, object]:
    pred = score >= threshold
    tp = int(((pred == 1) & (y_true == 1)).sum())
    tn = int(((pred == 0) & (y_true == 0)).sum())
    fp = int(((pred == 1) & (y_true == 0)).sum())
    fn = int(((pred == 0) & (y_true == 1)).sum())
    eps = 1e-12
    clipped = np.clip(score, eps, 1 - eps)
    return {
        "evaluation_set": label,
        "n": int(len(y_true)),
        "cases": int(y_true.sum()),
        "controls": int(len(y_true) - y_true.sum()),
        "prevalence": float(y_true.mean()) if len(y_true) else float("nan"),
        "roc_auc": _roc_auc(y_true, score),
        "average_precision": _average_precision(y_true, score),
        "brier_score": float(np.mean((score - y_true) ** 2)) if len(y_true) else float("nan"),
        "log_loss": float(-np.mean(y_true * np.log(clipped) + (1 - y_true) * np.log(1 - clipped))) if len(y_true) else float("nan"),
        "threshold": threshold,
        "sensitivity": tp / (tp + fn) if (tp + fn) else float("nan"),
        "specificity": tn / (tn + fp) if (tn + fp) else float("nan"),
        "ppv": tp / (tp + fp) if (tp + fp) else float("nan"),
        "npv": tn / (tn + fn) if (tn + fn) else float("nan"),
        "accuracy": (tp + tn) / len(y_true) if len(y_true) else float("nan"),
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
    }


def _calibration_table(y_true: np.ndarray, score: np.ndarray, *, bins: int = 10) -> pd.DataFrame:
    if len(y_true) == 0:
        return pd.DataFrame(columns=["bin", "n", "mean_predicted_probability", "observed_case_fraction"])
    order = np.argsort(score)
    groups = np.array_split(order, min(bins, len(order)))
    rows = []
    for index, group in enumerate(groups, start=1):
        rows.append(
            {
                "bin": index,
                "n": int(len(group)),
                "mean_predicted_probability": float(score[group].mean()),
                "observed_case_fraction": float(y_true[group].mean()),
            }
        )
    return pd.DataFrame(rows)


def _curve_points(y_true: np.ndarray, score: np.ndarray, curve: str) -> pd.DataFrame:
    order = np.argsort(-score)
    y_sorted = y_true[order]
    tp = np.cumsum(y_sorted)
    fp = np.cumsum(1 - y_sorted)
    positives = max(int(y_true.sum()), 1)
    negatives = max(int(len(y_true) - y_true.sum()), 1)
    if curve == "roc":
        return pd.DataFrame(
            {
                "x": np.r_[0.0, fp / negatives],
                "y": np.r_[0.0, tp / positives],
            }
        )
    precision = tp / np.maximum(tp + fp, 1)
    recall = tp / positives
    return pd.DataFrame({"x": np.r_[0.0, recall], "y": np.r_[1.0, precision]})


def _write_line_svg(df: pd.DataFrame, path: str, *, title: str, x_label: str, y_label: str, diagonal: bool = False) -> None:
    width, height = 560, 460
    left, right, top, bottom = 70, 30, 45, 65
    inner_w = width - left - right
    inner_h = height - top - bottom
    if df.empty:
        points = ""
    else:
        x = pd.to_numeric(df["x"], errors="coerce").fillna(0).clip(0, 1)
        y = pd.to_numeric(df["y"], errors="coerce").fillna(0).clip(0, 1)
        coords = [
            f"{left + float(xv) * inner_w:.2f},{top + inner_h - float(yv) * inner_h:.2f}"
            for xv, yv in zip(x, y, strict=False)
        ]
        points = " ".join(coords)
    diagonal_svg = (
        f'<line x1="{left}" y1="{top + inner_h}" x2="{left + inner_w}" y2="{top}" stroke="#9ca3af" stroke-dasharray="5 5" />'
        if diagonal
        else ""
    )
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
        '<rect width="100%" height="100%" fill="#fffdf7" />'
        f'<text x="{left}" y="26" font-size="22" font-family="Menlo, monospace">{title}</text>'
        f'<line x1="{left}" y1="{top + inner_h}" x2="{left + inner_w}" y2="{top + inner_h}" stroke="#1f2937" />'
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + inner_h}" stroke="#1f2937" />'
        f"{diagonal_svg}"
        f'<polyline points="{points}" fill="none" stroke="#0f766e" stroke-width="3" />'
        f'<text x="{width / 2:.1f}" y="{height - 20}" text-anchor="middle" font-size="13" font-family="Menlo, monospace">{x_label}</text>'
        f'<text x="18" y="{height / 2:.1f}" font-size="13" transform="rotate(-90 18 {height / 2:.1f})" font-family="Menlo, monospace">{y_label}</text>'
        "</svg>"
    )
    write_text(svg, path)


def _write_calibration_svg(calibration: pd.DataFrame, path: str) -> None:
    df = calibration.rename(columns={"mean_predicted_probability": "x", "observed_case_fraction": "y"})
    _write_line_svg(df, path, title="Clinical Model Calibration", x_label="Mean predicted probability", y_label="Observed case fraction", diagonal=True)


def _fit_model(train_df: pd.DataFrame, *, outcome_column: str, l2_penalty: float) -> tuple[ClinicalModel, np.ndarray]:
    specs = _fit_feature_specs(train_df)
    if not specs:
        raise ValueError("No usable clinical model features were found.")
    x = _feature_matrix(train_df, specs)
    y = pd.to_numeric(train_df[outcome_column], errors="coerce").astype(int).to_numpy()
    beta = _fit_logistic(x, y, l2_penalty=l2_penalty)
    train_prob = _sigmoid(x @ beta)
    threshold = _choose_threshold(y, train_prob)
    return ClinicalModel(beta=beta, features=specs, threshold=threshold, l2_penalty=l2_penalty), train_prob


def _predict(model: ClinicalModel, frame: pd.DataFrame) -> np.ndarray:
    return _sigmoid(_feature_matrix(frame, model.features) @ model.beta)


def _coefficient_table(model: ClinicalModel) -> pd.DataFrame:
    rows = [
        {
            "feature": "intercept",
            "source_column": "",
            "kind": "intercept",
            "reference": "",
            "beta": float(model.beta[0]),
            "odds_ratio": float(math.exp(model.beta[0])) if np.isfinite(model.beta[0]) else np.nan,
        }
    ]
    for beta, spec in zip(model.beta[1:], model.features, strict=False):
        rows.append(
            {
                "feature": spec.name,
                "source_column": spec.source_column,
                "kind": spec.kind,
                "reference": spec.reference,
                "beta": float(beta),
                "odds_ratio": float(math.exp(beta)) if np.isfinite(beta) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _cv_metrics(train_df: pd.DataFrame, *, outcome_column: str, l2_penalty: float) -> pd.DataFrame:
    if "cv_fold" not in train_df.columns:
        return pd.DataFrame()
    rows = []
    for fold in sorted(pd.to_numeric(train_df["cv_fold"], errors="coerce").dropna().astype(int).unique().tolist()):
        fit_df = train_df[pd.to_numeric(train_df["cv_fold"], errors="coerce") != fold].copy()
        valid_df = train_df[pd.to_numeric(train_df["cv_fold"], errors="coerce") == fold].copy()
        if fit_df.empty or valid_df.empty or fit_df[outcome_column].nunique() < 2 or valid_df[outcome_column].nunique() < 2:
            continue
        model, _ = _fit_model(fit_df, outcome_column=outcome_column, l2_penalty=l2_penalty)
        valid_prob = _predict(model, valid_df)
        y_valid = pd.to_numeric(valid_df[outcome_column], errors="coerce").astype(int).to_numpy()
        metric = _metrics(y_valid, valid_prob, threshold=model.threshold, label=f"cv_fold_{fold}")
        metric["fold"] = fold
        rows.append(metric)
    return pd.DataFrame(rows)


def _write_report(
    *,
    metrics: pd.DataFrame,
    cv_metrics: pd.DataFrame,
    coefficients: pd.DataFrame,
    eligibility_flag: str,
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
        "# Clinical-Only Prediction Model",
        "",
        f"- Eligibility flag: `{eligibility_flag}`",
        "- Model: L2-regularized logistic regression",
        "- Primary predictors: age, observation depth, condition-record depth, sex category, ancestry category, pre-index sepsis, pre-index renal injury",
        "- Peri-index and post-index variables are not used as predictors.",
        "",
        "## Metrics",
        "",
        md(metrics),
        "",
        "## Cross-Validation Metrics",
        "",
        md(cv_metrics),
        "",
        "## Coefficients",
        "",
        md(coefficients, limit=30),
        "",
    ]
    write_text("\n".join(body), path)


def run_clinical_model(
    config: ProjectConfig,
    paths: ProjectPaths,
    *,
    eligibility_flag: str = "primary_model_eligible",
    l2_penalty: float = 1.0,
) -> dict[str, pd.DataFrame | str]:
    model_input = read_table(clinical_model_input_path(paths))
    if eligibility_flag not in model_input.columns:
        raise ValueError(f"Clinical model input is missing eligibility flag: {eligibility_flag}")
    outcome_column = config.analysis.matched_outcome_column
    analysis = model_input[pd.to_numeric(model_input[eligibility_flag], errors="coerce").fillna(0) == 1].copy()
    analysis[outcome_column] = pd.to_numeric(analysis[outcome_column], errors="coerce")
    analysis = analysis[analysis[outcome_column].isin([0, 1])].copy()
    train = analysis[analysis["analysis_split"] == "train"].copy()
    test = analysis[analysis["analysis_split"] == "test"].copy()
    if train.empty or test.empty:
        raise ValueError("Clinical model requires non-empty train and test splits.")
    if train[outcome_column].nunique() < 2 or test[outcome_column].nunique() < 2:
        raise ValueError("Clinical model requires both cases and controls in train and test splits.")

    model, train_prob = _fit_model(train, outcome_column=outcome_column, l2_penalty=l2_penalty)
    test_prob = _predict(model, test)
    y_train = train[outcome_column].astype(int).to_numpy()
    y_test = test[outcome_column].astype(int).to_numpy()
    metrics = pd.DataFrame(
        [
            _metrics(y_train, train_prob, threshold=model.threshold, label="train"),
            _metrics(y_test, test_prob, threshold=model.threshold, label="test"),
        ]
    )
    cv_metrics = _cv_metrics(train, outcome_column=outcome_column, l2_penalty=l2_penalty)
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
            outcome_column,
            "analysis_split",
            "cv_fold",
            "match_group_id",
            eligibility_flag,
            "predicted_probability",
            "predicted_label",
        )
        if column in predictions.columns
    ]
    predictions = predictions[prediction_columns].copy()

    write_dataframe(metrics, clinical_model_metrics_path(paths))
    write_dataframe(cv_metrics, clinical_model_cv_metrics_path(paths))
    write_dataframe(coefficients, clinical_model_coefficients_path(paths))
    write_dataframe(predictions, clinical_model_predictions_path(paths))
    write_dataframe(calibration, clinical_model_calibration_path(paths))
    write_json(
        {
            "eligibility_flag": eligibility_flag,
            "l2_penalty": l2_penalty,
            "threshold": model.threshold,
            "features": [spec.__dict__ for spec in model.features],
        },
        join_path(clinical_model_dir(paths), "model_spec.json"),
    )
    _write_line_svg(_curve_points(y_test, test_prob, "roc"), clinical_model_roc_svg_path(paths), title="Clinical Model ROC", x_label="False positive rate", y_label="True positive rate", diagonal=True)
    _write_line_svg(_curve_points(y_test, test_prob, "pr"), clinical_model_pr_svg_path(paths), title="Clinical Model Precision-Recall", x_label="Recall", y_label="Precision")
    _write_calibration_svg(calibration, clinical_model_calibration_svg_path(paths))
    _write_report(
        metrics=metrics,
        cv_metrics=cv_metrics,
        coefficients=coefficients,
        eligibility_flag=eligibility_flag,
        path=clinical_model_report_path(paths),
    )
    return {
        "metrics": metrics,
        "cv_metrics": cv_metrics,
        "coefficients": coefficients,
        "predictions": predictions,
        "calibration": calibration,
        "metrics_path": clinical_model_metrics_path(paths),
        "cv_metrics_path": clinical_model_cv_metrics_path(paths),
        "coefficients_path": clinical_model_coefficients_path(paths),
        "predictions_path": clinical_model_predictions_path(paths),
        "calibration_path": clinical_model_calibration_path(paths),
        "report_path": clinical_model_report_path(paths),
        "roc_svg_path": clinical_model_roc_svg_path(paths),
        "pr_svg_path": clinical_model_pr_svg_path(paths),
        "calibration_svg_path": clinical_model_calibration_svg_path(paths),
    }


__all__ = [
    "clinical_model_calibration_path",
    "clinical_model_calibration_svg_path",
    "clinical_model_coefficients_path",
    "clinical_model_cv_metrics_path",
    "clinical_model_metrics_path",
    "clinical_model_predictions_path",
    "clinical_model_pr_svg_path",
    "clinical_model_report_path",
    "clinical_model_roc_svg_path",
    "run_clinical_model",
]
