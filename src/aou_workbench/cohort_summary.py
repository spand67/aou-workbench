"""Clinical and demographic summary tables for unmatched and matched analyses."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, ttest_ind

from .config import ProjectConfig
from .io_utils import write_dataframe, write_text
from .paths import ProjectPaths, join_path
from .sample_restriction import restrict_frame_for_gwas


def cohort_summary_table_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "clinical_demographic_summary.tsv")


def cohort_summary_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "clinical_demographic_summary.md")


_UNMATCHED_CASES = "unmatched_rhabdo_wgs"
_UNMATCHED_CONTROLS = "unmatched_non_rhabdo_wgs"
_MATCHED_CASES = "matched_rhabdo_wgs"
_MATCHED_CONTROLS = "matched_controls_wgs"


def _summary_columns() -> list[str]:
    return [
        "section",
        "variable",
        _UNMATCHED_CASES,
        _UNMATCHED_CONTROLS,
        "unmatched_smd",
        "unmatched_p_value",
        _MATCHED_CASES,
        _MATCHED_CONTROLS,
        "matched_smd",
        "matched_p_value",
    ]


def _empty_value() -> str:
    return ""


def _format_count(frame: pd.DataFrame) -> str:
    return str(int(frame["person_id"].astype(str).nunique()))


def _numeric_series(frame: pd.DataFrame, column: str) -> pd.Series:
    return pd.to_numeric(frame.get(column), errors="coerce").dropna()


def _binary_series(frame: pd.DataFrame, column: str) -> pd.Series:
    return (_numeric_series(frame, column) >= 1).astype(float)


def _categorical_series(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame.get(column, pd.Series(dtype=object)).dropna().astype(str)


def _format_mean_sd(frame: pd.DataFrame, column: str) -> str:
    values = _numeric_series(frame, column)
    if values.empty:
        return _empty_value()
    sd = float(values.std(ddof=1)) if len(values) > 1 else 0.0
    return f"{values.mean():.1f} ({sd:.1f})"


def _format_binary(frame: pd.DataFrame, column: str) -> str:
    values = _binary_series(frame, column)
    if values.empty:
        return _empty_value()
    count = int(values.sum())
    return f"{count} ({(100.0 * count / len(values)):.1f}%)"


def _format_category(frame: pd.DataFrame, column: str, value: str) -> str:
    series = _categorical_series(frame, column)
    if series.empty:
        return _empty_value()
    count = int((series == value).sum())
    return f"{count} ({(100.0 * count / len(series)):.1f}%)"


def _format_p_value(value: float) -> str:
    if pd.isna(value):
        return _empty_value()
    if value < 1e-4:
        return f"{value:.2e}"
    return f"{value:.4f}"


def _format_smd(value: float) -> str:
    if pd.isna(value):
        return _empty_value()
    return f"{value:.3f}"


def _continuous_smd(left: pd.Series, right: pd.Series) -> float:
    if left.empty or right.empty:
        return np.nan
    left_sd = float(left.std(ddof=1)) if len(left) > 1 else 0.0
    right_sd = float(right.std(ddof=1)) if len(right) > 1 else 0.0
    pooled = np.sqrt((left_sd**2 + right_sd**2) / 2.0)
    if pooled == 0:
        return 0.0
    return float((left.mean() - right.mean()) / pooled)


def _binary_smd(left: pd.Series, right: pd.Series) -> float:
    if left.empty or right.empty:
        return np.nan
    left_p = float(left.mean())
    right_p = float(right.mean())
    pooled = np.sqrt((left_p * (1.0 - left_p) + right_p * (1.0 - right_p)) / 2.0)
    if pooled == 0:
        return 0.0
    return float((left_p - right_p) / pooled)


def _categorical_smd(left: pd.Series, right: pd.Series, categories: list[str]) -> float:
    if left.empty or right.empty or not categories:
        return np.nan
    smds = []
    for category in categories:
        smds.append(
            abs(
                _binary_smd(
                    (left == category).astype(float),
                    (right == category).astype(float),
                )
            )
        )
    return float(max(value for value in smds if pd.notna(value))) if smds else np.nan


def _continuous_p_value(left: pd.Series, right: pd.Series) -> float:
    if len(left) < 2 or len(right) < 2:
        return np.nan
    return float(ttest_ind(left, right, equal_var=False, nan_policy="omit").pvalue)


def _binary_p_value(left: pd.Series, right: pd.Series) -> float:
    if left.empty or right.empty:
        return np.nan
    left_count = int(left.sum())
    right_count = int(right.sum())
    table = np.array(
        [
            [left_count, int(len(left) - left_count)],
            [right_count, int(len(right) - right_count)],
        ]
    )
    if (table < 5).any():
        return float(fisher_exact(table)[1])
    return float(chi2_contingency(table, correction=False)[1])


def _categorical_p_value(left: pd.Series, right: pd.Series, categories: list[str]) -> float:
    if left.empty or right.empty or not categories:
        return np.nan
    table = np.array(
        [
            [int((left == category).sum()) for category in categories],
            [int((right == category).sum()) for category in categories],
        ]
    )
    if table.sum() == 0:
        return np.nan
    return float(chi2_contingency(table, correction=False)[1])


def _render_markdown(summary_df: pd.DataFrame, *, has_max_unrelated: bool) -> str:
    if summary_df.empty:
        return "# Clinical and Demographic Summary\n\n_No rows available._\n"
    try:
        table = summary_df.to_markdown(index=False)
    except ImportError:
        table = summary_df.to_string(index=False)
    notes = [
        "# Clinical and Demographic Summary",
        "",
        "Columns represent the two WGS-restricted analyses:",
        "- Unmatched OMOP rhabdo vs non-rhabdo",
        "- Matched OMOP rhabdo vs matched controls",
    ]
    if has_max_unrelated:
        notes.append("- Max-unrelated restriction was applied in addition to WGS presence")
    notes.extend(
        [
            "",
            "Balance diagnostics emphasize standardized mean difference (SMD); values with |SMD| > 0.1 merit attention.",
            "",
            table,
            "",
        ]
    )
    return "\n".join(notes)


def summarize_clinical_demographics(
    config: ProjectConfig,
    built_df: pd.DataFrame,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
) -> pd.DataFrame:
    built_wgs = restrict_frame_for_gwas(config, built_df, require_wgs=True)
    matched_wgs = restrict_frame_for_gwas(config, matched_df, require_wgs=True)

    groups = {
        _UNMATCHED_CASES: built_wgs[built_wgs["rhabdo_case"].fillna(0).astype(int) == 1].copy(),
        _UNMATCHED_CONTROLS: built_wgs[built_wgs["rhabdo_case"].fillna(0).astype(int) == 0].copy(),
        _MATCHED_CASES: matched_wgs[matched_wgs[config.analysis.matched_outcome_column].fillna(0).astype(int) == 1].copy(),
        _MATCHED_CONTROLS: matched_wgs[matched_wgs[config.analysis.matched_outcome_column].fillna(0).astype(int) == 0].copy(),
    }

    rows: list[dict[str, str]] = []

    def add_row(
        section: str,
        variable: str,
        formatter: Callable[[pd.DataFrame], str],
        unmatched_metrics: tuple[float, float] | None = None,
        matched_metrics: tuple[float, float] | None = None,
    ) -> None:
        row = {
            "section": section,
            "variable": variable,
            _UNMATCHED_CASES: formatter(groups[_UNMATCHED_CASES]),
            _UNMATCHED_CONTROLS: formatter(groups[_UNMATCHED_CONTROLS]),
            "unmatched_smd": _format_smd(unmatched_metrics[0]) if unmatched_metrics else _empty_value(),
            "unmatched_p_value": _format_p_value(unmatched_metrics[1]) if unmatched_metrics else _empty_value(),
            _MATCHED_CASES: formatter(groups[_MATCHED_CASES]),
            _MATCHED_CONTROLS: formatter(groups[_MATCHED_CONTROLS]),
            "matched_smd": _format_smd(matched_metrics[0]) if matched_metrics else _empty_value(),
            "matched_p_value": _format_p_value(matched_metrics[1]) if matched_metrics else _empty_value(),
        }
        rows.append(row)

    add_row("Sample size", "N", _format_count)

    age_column = "age_at_index" if "age_at_index" in built_wgs.columns else "age_raw"
    unmatched_age_left = _numeric_series(groups[_UNMATCHED_CASES], age_column)
    unmatched_age_right = _numeric_series(groups[_UNMATCHED_CONTROLS], age_column)
    matched_age_left = _numeric_series(groups[_MATCHED_CASES], age_column)
    matched_age_right = _numeric_series(groups[_MATCHED_CONTROLS], age_column)
    add_row(
        "Demographics",
        "Age at index, mean (SD)",
        lambda frame: _format_mean_sd(frame, age_column),
        unmatched_metrics=(
            _continuous_smd(unmatched_age_left, unmatched_age_right),
            _continuous_p_value(unmatched_age_left, unmatched_age_right),
        ),
        matched_metrics=(
            _continuous_smd(matched_age_left, matched_age_right),
            _continuous_p_value(matched_age_left, matched_age_right),
        ),
    )

    unmatched_female_left = _binary_series(groups[_UNMATCHED_CASES], "is_female")
    unmatched_female_right = _binary_series(groups[_UNMATCHED_CONTROLS], "is_female")
    matched_female_left = _binary_series(groups[_MATCHED_CASES], "is_female")
    matched_female_right = _binary_series(groups[_MATCHED_CONTROLS], "is_female")
    add_row(
        "Demographics",
        "Female sex, n (%)",
        lambda frame: _format_binary(frame, "is_female"),
        unmatched_metrics=(
            _binary_smd(unmatched_female_left, unmatched_female_right),
            _binary_p_value(unmatched_female_left, unmatched_female_right),
        ),
        matched_metrics=(
            _binary_smd(matched_female_left, matched_female_right),
            _binary_p_value(matched_female_left, matched_female_right),
        ),
    )

    ancestry_values = sorted(
        {
            value
            for value in pd.concat(
                [
                    groups[_UNMATCHED_CASES].get("ancestry_pred", pd.Series(dtype=object)).astype(str),
                    groups[_UNMATCHED_CONTROLS].get("ancestry_pred", pd.Series(dtype=object)).astype(str),
                    groups[_MATCHED_CASES].get("ancestry_pred", pd.Series(dtype=object)).astype(str),
                    groups[_MATCHED_CONTROLS].get("ancestry_pred", pd.Series(dtype=object)).astype(str),
                ],
                ignore_index=True,
            ).dropna()
            if str(value).strip()
        }
    )
    unmatched_ancestry_left = _categorical_series(groups[_UNMATCHED_CASES], "ancestry_pred")
    unmatched_ancestry_right = _categorical_series(groups[_UNMATCHED_CONTROLS], "ancestry_pred")
    matched_ancestry_left = _categorical_series(groups[_MATCHED_CASES], "ancestry_pred")
    matched_ancestry_right = _categorical_series(groups[_MATCHED_CONTROLS], "ancestry_pred")
    add_row(
        "Demographics",
        "Ancestry distribution, overall",
        lambda _frame: _empty_value(),
        unmatched_metrics=(
            _categorical_smd(unmatched_ancestry_left, unmatched_ancestry_right, ancestry_values),
            _categorical_p_value(unmatched_ancestry_left, unmatched_ancestry_right, ancestry_values),
        ),
        matched_metrics=(
            _categorical_smd(matched_ancestry_left, matched_ancestry_right, ancestry_values),
            _categorical_p_value(matched_ancestry_left, matched_ancestry_right, ancestry_values),
        ),
    )
    for ancestry in ancestry_values:
        add_row(
            "Demographics",
            f"Ancestry: {ancestry}, n (%)",
            lambda frame, ancestry=ancestry: _format_category(frame, "ancestry_pred", ancestry),
            unmatched_metrics=(
                _binary_smd(
                    (unmatched_ancestry_left == ancestry).astype(float),
                    (unmatched_ancestry_right == ancestry).astype(float),
                ),
                np.nan,
            ),
            matched_metrics=(
                _binary_smd(
                    (matched_ancestry_left == ancestry).astype(float),
                    (matched_ancestry_right == ancestry).astype(float),
                ),
                np.nan,
            ),
        )

    clinical_columns = [
        column
        for column in config.phenotype.clinical_cofactor_columns
        if column in built_wgs.columns or column in matched_wgs.columns
    ]
    for column in clinical_columns:
        unmatched_binary_left = _binary_series(groups[_UNMATCHED_CASES], column)
        unmatched_binary_right = _binary_series(groups[_UNMATCHED_CONTROLS], column)
        matched_binary_left = _binary_series(groups[_MATCHED_CASES], column)
        matched_binary_right = _binary_series(groups[_MATCHED_CONTROLS], column)
        add_row(
            "Clinical",
            column.replace("_", " ").title() + ", n (%)",
            lambda frame, column=column: _format_binary(frame, column),
            unmatched_metrics=(
                _binary_smd(unmatched_binary_left, unmatched_binary_right),
                _binary_p_value(unmatched_binary_left, unmatched_binary_right),
            ),
            matched_metrics=(
                _binary_smd(matched_binary_left, matched_binary_right),
                _binary_p_value(matched_binary_left, matched_binary_right),
            ),
        )

    summary = pd.DataFrame(rows, columns=_summary_columns())
    table_path = cohort_summary_table_path(paths)
    report_path = cohort_summary_report_path(paths)
    write_dataframe(summary, table_path)
    write_text(
        _render_markdown(summary, has_max_unrelated=bool(config.workbench.max_unrelated_path)),
        report_path,
    )
    return summary


__all__ = [
    "cohort_summary_report_path",
    "cohort_summary_table_path",
    "summarize_clinical_demographics",
]
