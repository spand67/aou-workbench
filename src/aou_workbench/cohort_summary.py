"""Clinical and demographic summary tables for unmatched and matched analyses."""

from __future__ import annotations

import hashlib
import html
from collections.abc import Callable

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, t, ttest_ind

from .config import ProjectConfig
from .cohort import load_clinical_cofactor_events
from .io_utils import write_dataframe, write_text
from .paths import ProjectPaths, join_path
from .sample_restriction import restrict_frame_for_gwas


def cohort_summary_table_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "clinical_demographic_summary.tsv")


def cohort_summary_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "clinical_demographic_summary.md")


def consort_counts_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "consort_counts.tsv")


def consort_counts_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "consort_counts.md")


def matched_table1_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "table1_clinical_matched.tsv")


def matched_table1_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "table1_clinical_matched.md")


def critical_illness_summary_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "critical_illness_summary.tsv")


def critical_illness_summary_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "critical_illness_summary.md")


def case_cofactor_prior_timing_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "case_cofactor_prior_timing.tsv")


def case_cofactor_prior_timing_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "case_cofactor_prior_timing.md")


def case_cofactor_prior_timing_histogram_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "case_cofactor_prior_timing_histogram.svg")


def missingness_summary_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "missingness_summary.tsv")


def missingness_summary_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "missingness_summary.md")


def model_split_summary_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "model_split_summary.tsv")


def model_split_summary_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "model_split_summary.md")


def model_eligibility_summary_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "model_eligibility_summary.tsv")


def model_eligibility_summary_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "model_eligibility_summary.md")


def split_table1_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "table1_clinical_by_split.tsv")


def split_table1_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "table1_clinical_by_split.md")


def clinical_model_input_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "clinical_model_input.tsv")


def clinical_characterization_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "clinical_characterization_report.md")


_UNMATCHED_CASES = "unmatched_rhabdo_wgs"
_UNMATCHED_CONTROLS = "unmatched_non_rhabdo_wgs"
_MATCHED_CASES = "matched_rhabdo_wgs"
_MATCHED_CONTROLS = "matched_controls_wgs"
_TRAIN_FRACTION = 0.8
_CV_FOLDS = 5
_PRIOR_TIMING_BINS = [
    ("same_day", 0, 0),
    ("1_to_7_days_before", 1, 7),
    ("8_to_30_days_before", 8, 30),
    ("31_to_90_days_before", 31, 90),
    ("91_to_365_days_before", 91, 365),
    ("more_than_365_days_before", 366, None),
]
_MISSING_LABELS = {
    "",
    "nan",
    "none",
    "null",
    "na",
    "n/a",
    "missing",
    "no matching concept",
    "unknown",
    "skip",
    "pmi: skip",
    "prefer not to answer",
    "pmi: prefer not to answer",
}


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


def _normalize_category_values(series: pd.Series) -> pd.Series:
    normalized = series.astype(str).str.lower().str.strip()
    normalized = normalized.replace({value: "missing" for value in _MISSING_LABELS})
    normalized = normalized.replace(
        {
            "f": "female",
            "2": "female",
            "m": "male",
            "1": "male",
            "other": "other_or_unknown",
            "other/unknown": "other_or_unknown",
        }
    )
    return normalized


def _sex_category_series(frame: pd.DataFrame) -> pd.Series:
    if "sex_category" in frame.columns:
        values = _normalize_category_values(frame["sex_category"])
    elif "gender_concept_name" in frame.columns:
        values = _normalize_category_values(frame["gender_concept_name"])
    elif "is_female" in frame.columns:
        numeric = pd.to_numeric(frame["is_female"], errors="coerce")
        values = pd.Series("missing", index=frame.index, dtype=object)
        values.loc[numeric == 1] = "female"
        values.loc[numeric == 0] = "male"
    else:
        values = pd.Series("missing", index=frame.index, dtype=object)
    values = values.where(values.isin({"female", "male", "missing"}), "other_or_unknown")
    return values


def _has_sex_information(frame: pd.DataFrame) -> bool:
    return any(column in frame.columns for column in ("sex_category", "gender_concept_name", "is_female"))


def _ordered_categories(series_list: list[pd.Series], preferred: list[str]) -> list[str]:
    observed = {
        value
        for series in series_list
        for value in series.dropna().astype(str).tolist()
        if value.strip()
    }
    return [value for value in preferred if value in observed] + sorted(observed.difference(preferred))


def _ensure_categories(categories: list[str], required: list[str]) -> list[str]:
    return categories + [category for category in required if category not in categories]


def _format_category_from_series(series: pd.Series, value: str) -> str:
    if series.empty:
        return _empty_value()
    count = int((series == value).sum())
    return f"{count} ({(100.0 * count / len(series)):.1f}%)"


def _display_category(value: str) -> str:
    return value.replace("_or_", "/").replace("_", " ")


def _category_series_with_missing(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series("missing", index=frame.index, dtype=object)
    values = frame[column].astype("string").str.strip()
    missing = values.isna() | values.str.lower().isin(_MISSING_LABELS)
    return values.where(~missing, "missing").astype(str)


def _missing_mask(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(True, index=frame.index)
    values = frame[column]
    mask = values.isna()
    if pd.api.types.is_object_dtype(values) or pd.api.types.is_string_dtype(values):
        text = values.astype("string").str.strip().str.lower()
        mask = mask | text.isin(_MISSING_LABELS)
    return mask


def _stable_hash_score(value: object) -> float:
    digest = hashlib.sha256(str(value).encode("utf-8")).hexdigest()
    return int(digest[:16], 16) / float(16**16 - 1)


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


def _format_effect(value: float) -> str:
    if pd.isna(value):
        return _empty_value()
    if abs(value) >= 100 or (abs(value) < 0.01 and value != 0):
        return f"{value:.2e}"
    return f"{value:.3f}"


def _format_ci(low: float, high: float) -> str:
    if pd.isna(low) or pd.isna(high):
        return _empty_value()
    return f"{_format_effect(low)}, {_format_effect(high)}"


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


def _mean_difference_ci(left: pd.Series, right: pd.Series) -> tuple[float, float, float]:
    if len(left) < 2 or len(right) < 2:
        return np.nan, np.nan, np.nan
    diff = float(left.mean() - right.mean())
    left_var = float(left.var(ddof=1))
    right_var = float(right.var(ddof=1))
    se_sq = left_var / len(left) + right_var / len(right)
    if se_sq <= 0:
        return diff, diff, diff
    se = float(np.sqrt(se_sq))
    numerator = se_sq**2
    denominator = ((left_var / len(left)) ** 2 / (len(left) - 1)) + ((right_var / len(right)) ** 2 / (len(right) - 1))
    df = numerator / denominator if denominator > 0 else min(len(left), len(right)) - 1
    critical = float(t.ppf(0.975, df)) if df > 0 else np.nan
    if pd.isna(critical):
        return diff, np.nan, np.nan
    return diff, diff - critical * se, diff + critical * se


def _binary_p_value(left: pd.Series, right: pd.Series) -> float:
    if left.empty or right.empty:
        return np.nan
    left_count = int(left.sum())
    right_count = int(right.sum())
    if left_count == 0 and right_count == 0:
        return np.nan
    table = np.array(
        [
            [left_count, int(len(left) - left_count)],
            [right_count, int(len(right) - right_count)],
        ]
    )
    if (table < 5).any():
        return float(fisher_exact(table)[1])
    return float(chi2_contingency(table, correction=False)[1])


def _odds_ratio_ci(left: pd.Series, right: pd.Series) -> tuple[float, float, float]:
    if left.empty or right.empty:
        return np.nan, np.nan, np.nan
    a = float(left.sum())
    b = float(len(left) - left.sum())
    c = float(right.sum())
    d = float(len(right) - right.sum())
    if a == 0 and c == 0:
        return np.nan, np.nan, np.nan
    if min(a, b, c, d) == 0:
        a += 0.5
        b += 0.5
        c += 0.5
        d += 0.5
    odds_ratio = (a * d) / (b * c) if b * c else np.nan
    if pd.isna(odds_ratio) or odds_ratio <= 0:
        return odds_ratio, np.nan, np.nan
    se = np.sqrt((1 / a) + (1 / b) + (1 / c) + (1 / d))
    low = float(np.exp(np.log(odds_ratio) - 1.96 * se))
    high = float(np.exp(np.log(odds_ratio) + 1.96 * se))
    return float(odds_ratio), low, high


def _categorical_p_value(left: pd.Series, right: pd.Series, categories: list[str]) -> float:
    if left.empty or right.empty or not categories:
        return np.nan
    table = np.array(
        [
            [int((left == category).sum()) for category in categories],
            [int((right == category).sum()) for category in categories],
        ]
    )
    table = table[:, table.sum(axis=0) > 0]
    if table.sum() == 0 or table.shape[1] < 2 or (table.sum(axis=1) == 0).any():
        return np.nan
    try:
        return float(chi2_contingency(table, correction=False)[1])
    except ValueError:
        return np.nan


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


def _markdown_table(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False)
    except ImportError:
        return df.to_string(index=False)


def _write_markdown_table(df: pd.DataFrame, path: str, title: str) -> None:
    write_text(f"# {title}\n\n{_markdown_table(df)}\n", path)


def build_consort_counts(built_df: pd.DataFrame, matched_df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        [
            ("Built cohort rows", len(built_df)),
            (">=2 OMOP condition dates", int(built_df.get("eligible_ehr_denominator", pd.Series(dtype=int)).sum())),
            ("Broad rhabdo cases", int(built_df.get("broad_rhabdo_case", pd.Series(dtype=int)).sum())),
            ("Definite rhabdo cases", int(built_df.get("definite_rhabdo_case", pd.Series(dtype=int)).sum())),
            ("CK >5000 without rhabdo, excluded", int((built_df.get("case_tier", pd.Series(dtype=str)) == "indeterminate_ck_only").sum())),
            ("Eligible controls", int(built_df.get("eligible_control", pd.Series(dtype=int)).sum())),
            ("Matched analytic rows", len(matched_df)),
            ("Matched cases", int((matched_df.get("analysis_case", pd.Series(dtype=int)) == 1).sum())),
            ("Matched controls", int((matched_df.get("analysis_case", pd.Series(dtype=int)) == 0).sum())),
        ],
        columns=["step", "n"],
    )


def _matched_groups(config: ProjectConfig, matched_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    outcome = config.analysis.matched_outcome_column
    cases = matched_df[matched_df[outcome].fillna(0).astype(int) == 1].copy()
    controls = matched_df[matched_df[outcome].fillna(0).astype(int) == 0].copy()
    return cases, controls


def build_matched_table1(config: ProjectConfig, matched_df: pd.DataFrame) -> pd.DataFrame:
    cases, controls = _matched_groups(config, matched_df)
    rows: list[dict[str, object]] = [
        {
            "section": "Sample size",
            "variable": "N",
            "matched_cases": _format_count(cases),
            "matched_controls": _format_count(controls),
            "effect_measure": "",
            "effect_size": "",
            "ci_95": "",
            "smd": "",
            "p_value": "",
        }
    ]

    def add_continuous(section: str, label: str, column: str) -> None:
        if column not in matched_df.columns:
            return
        left = _numeric_series(cases, column)
        right = _numeric_series(controls, column)
        effect, low, high = _mean_difference_ci(left, right)
        rows.append(
            {
                "section": section,
                "variable": label,
                "matched_cases": _format_mean_sd(cases, column),
                "matched_controls": _format_mean_sd(controls, column),
                "effect_measure": "mean_difference",
                "effect_size": _format_effect(effect),
                "ci_95": _format_ci(low, high),
                "smd": _format_smd(_continuous_smd(left, right)),
                "p_value": _format_p_value(_continuous_p_value(left, right)),
            }
        )

    def add_binary(section: str, label: str, column: str) -> None:
        if column not in matched_df.columns:
            return
        left = _binary_series(cases, column)
        right = _binary_series(controls, column)
        add_binary_values(
            section,
            label,
            left,
            right,
            _format_binary(cases, column),
            _format_binary(controls, column),
        )

    def add_binary_values(
        section: str,
        label: str,
        left: pd.Series,
        right: pd.Series,
        case_text: str,
        control_text: str,
    ) -> None:
        effect, low, high = _odds_ratio_ci(left, right)
        rows.append(
            {
                "section": section,
                "variable": label,
                "matched_cases": case_text,
                "matched_controls": control_text,
                "effect_measure": "odds_ratio",
                "effect_size": _format_effect(effect),
                "ci_95": _format_ci(low, high),
                "smd": _format_smd(_binary_smd(left, right)),
                "p_value": _format_p_value(_binary_p_value(left, right)),
            }
        )

    add_continuous("Demographics", "Age at index, mean (SD)", "age_at_index")
    if _has_sex_information(matched_df):
        sex_cases = _sex_category_series(cases)
        sex_controls = _sex_category_series(controls)
        sex_required = ["female", "male", "other_or_unknown", "missing"]
        sex_categories = _ensure_categories(_ordered_categories([sex_cases, sex_controls], sex_required), sex_required)
        rows.append(
            {
                "section": "Demographics",
                "variable": "Sex category distribution, overall",
                "matched_cases": "",
                "matched_controls": "",
                "effect_measure": "max_category_smd",
                "effect_size": _format_smd(_categorical_smd(sex_cases, sex_controls, sex_categories)),
                "ci_95": "",
                "smd": _format_smd(_categorical_smd(sex_cases, sex_controls, sex_categories)),
                "p_value": _format_p_value(_categorical_p_value(sex_cases, sex_controls, sex_categories)),
            }
        )
        for category in sex_categories:
            left = (sex_cases == category).astype(float)
            right = (sex_controls == category).astype(float)
            add_binary_values(
                "Demographics",
                f"Sex: {_display_category(category)}, n (%)",
                left,
                right,
                _format_category_from_series(sex_cases, category),
                _format_category_from_series(sex_controls, category),
            )

    ancestry = "ancestry_pred"
    if ancestry in matched_df.columns:
        left = _category_series_with_missing(cases, ancestry)
        right = _category_series_with_missing(controls, ancestry)
        categories = _ordered_categories(
            [left, right],
            ["AFR", "AMR", "EAS", "EUR", "MID", "SAS", "other", "missing"],
        )
        categories = _ensure_categories(categories, ["missing"])
        rows.append(
            {
                "section": "Demographics",
                "variable": "Ancestry distribution, overall",
                "matched_cases": "",
                "matched_controls": "",
                "effect_measure": "max_category_smd",
                "effect_size": _format_smd(_categorical_smd(left, right, categories)),
                "ci_95": "",
                "smd": _format_smd(_categorical_smd(left, right, categories)),
                "p_value": _format_p_value(_categorical_p_value(left, right, categories)),
            }
        )
        for category in categories:
            left_binary = (left == category).astype(float)
            right_binary = (right == category).astype(float)
            add_binary_values(
                "Demographics",
                f"Ancestry: {category}, n (%)",
                left_binary,
                right_binary,
                _format_category_from_series(left, category),
                _format_category_from_series(right, category),
            )

    clinical_columns = [
        "crush_injury",
        "preindex_crush_injury",
        "periindex_crush_injury",
        "sepsis",
        "remote_preindex_sepsis",
        "near_preindex_sepsis",
        "preindex_sepsis",
        "periindex_sepsis",
        "postindex_sepsis",
        "renal_injury",
        "remote_preindex_renal_injury",
        "near_preindex_renal_injury",
        "preindex_renal_injury",
        "periindex_renal_injury",
        "postindex_renal_injury",
    ]
    for column in clinical_columns:
        add_binary("Clinical", column.replace("_", " ").title() + ", n (%)", column)

    return pd.DataFrame(rows)


def build_critical_illness_summary(config: ProjectConfig, built_df: pd.DataFrame, matched_df: pd.DataFrame) -> pd.DataFrame:
    cases, controls = _matched_groups(config, matched_df)
    broad_flag = (
        built_df["broad_rhabdo_case"].fillna(0).astype(int)
        if "broad_rhabdo_case" in built_df.columns
        else pd.Series(0, index=built_df.index)
    )
    groups = {
        "built_broad_cases": built_df[broad_flag == 1].copy(),
        "matched_cases": cases,
        "matched_controls": controls,
    }
    variables = [
        "definite_rhabdo_case",
        "sepsis",
        "remote_preindex_sepsis",
        "near_preindex_sepsis",
        "preindex_sepsis",
        "periindex_sepsis",
        "postindex_sepsis",
        "renal_injury",
        "remote_preindex_renal_injury",
        "near_preindex_renal_injury",
        "preindex_renal_injury",
        "periindex_renal_injury",
        "postindex_renal_injury",
    ]
    rows: list[dict[str, object]] = []
    for group_name, frame in groups.items():
        n = int(frame["person_id"].astype(str).nunique()) if "person_id" in frame.columns else len(frame)
        rows.append({"group": group_name, "variable": "N", "n": n, "pct": ""})
        for variable in variables:
            if variable not in frame.columns:
                continue
            values = _binary_series(frame, variable)
            count = int(values.sum())
            rows.append(
                {
                    "group": group_name,
                    "variable": variable,
                    "n": count,
                    "pct": f"{(100.0 * count / len(values)):.1f}" if len(values) else "",
                }
            )
        for label, columns in {
            "remote_preindex_sepsis_or_renal_injury": ("remote_preindex_sepsis", "remote_preindex_renal_injury"),
            "near_preindex_sepsis_or_renal_injury": ("near_preindex_sepsis", "near_preindex_renal_injury"),
            "preindex_sepsis_or_renal_injury": ("preindex_sepsis", "preindex_renal_injury"),
            "periindex_sepsis_or_renal_injury": ("periindex_sepsis", "periindex_renal_injury"),
        }.items():
            if all(column in frame.columns for column in columns):
                values = ((_binary_series(frame, columns[0]) >= 1) | (_binary_series(frame, columns[1]) >= 1)).astype(float)
                count = int(values.sum())
                rows.append(
                    {
                        "group": group_name,
                        "variable": label,
                        "n": count,
                        "pct": f"{(100.0 * count / len(values)):.1f}" if len(values) else "",
                    }
                )
    return pd.DataFrame(rows)


def _quantile_text(values: pd.Series, q: float) -> str:
    if values.empty:
        return ""
    return f"{float(values.quantile(q)):.1f}"


def _svg_text(value: object) -> str:
    return html.escape(str(value), quote=True)


def _prior_timing_placeholder_svg(message: str) -> str:
    return (
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="240">'
        '<rect width="100%" height="100%" fill="#fffdf7" />'
        '<text x="40" y="80" font-size="26" font-family="Menlo, monospace">'
        "Case Cofactor Prior Timing"
        "</text>"
        f'<text x="40" y="130" font-size="16" font-family="Menlo, monospace">{_svg_text(message)}</text>'
        "</svg>"
    )


def _prior_window_label(window: str) -> str:
    return {
        "same_day": "0",
        "1_to_7_days_before": "1-7",
        "8_to_30_days_before": "8-30",
        "31_to_90_days_before": "31-90",
        "91_to_365_days_before": "91-365",
        "more_than_365_days_before": ">365",
    }.get(window, window)


def _prior_group_label(group: str) -> str:
    return {
        "built_broad_cases": "Built broad cases",
        "matched_cases": "Matched cases",
    }.get(group, group.replace("_", " ").title())


def _prior_cofactor_label(cofactor: str) -> str:
    return {"sepsis": "Sepsis", "renal_injury": "Renal injury"}.get(cofactor, cofactor.replace("_", " ").title())


def _case_index_frame(config: ProjectConfig, frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty or "person_id" not in frame.columns or "index_date" not in frame.columns:
        return pd.DataFrame(columns=["person_id", "index_date"])
    outcome = config.analysis.matched_outcome_column
    if outcome in frame.columns:
        case_mask = pd.to_numeric(frame[outcome], errors="coerce").fillna(0) == 1
    elif "broad_rhabdo_case" in frame.columns:
        case_mask = pd.to_numeric(frame["broad_rhabdo_case"], errors="coerce").fillna(0) == 1
    else:
        case_mask = pd.Series(False, index=frame.index)
    cases = frame[case_mask][["person_id", "index_date"]].copy()
    cases["person_id"] = cases["person_id"].astype(str)
    cases["index_date"] = pd.to_datetime(cases["index_date"], errors="coerce")
    return cases.dropna(subset=["index_date"]).drop_duplicates("person_id")


def _empty_prior_timing_row(
    group: str,
    cofactor: str,
    window: str,
    n_cases: int,
    n_cases_with_prior: int,
    n_cases_in_window: int,
) -> dict[str, object]:
    pct_all = f"{(100.0 * n_cases_in_window / n_cases):.1f}" if n_cases else ""
    pct_prior = f"{(100.0 * n_cases_in_window / n_cases_with_prior):.1f}" if n_cases_with_prior else ""
    return {
        "group": group,
        "cofactor": cofactor,
        "timing_basis": "nearest_prior_event_per_case",
        "window": window,
        "n_cases": n_cases,
        "n_cases_with_prior_event": n_cases_with_prior,
        "n_cases_in_window": n_cases_in_window,
        "pct_of_all_cases": pct_all,
        "pct_of_cases_with_prior_event": pct_prior,
        "median_days_before_index": "",
        "q1_days_before_index": "",
        "q3_days_before_index": "",
    }


def build_case_cofactor_prior_timing_summary(
    config: ProjectConfig,
    built_df: pd.DataFrame,
    matched_df: pd.DataFrame,
    events: pd.DataFrame | None = None,
) -> pd.DataFrame:
    cofactor_names = {rule.name for rule in config.phenotype.clinical_cofactors}
    targets = [name for name in ("sepsis", "renal_injury") if name in cofactor_names]
    if not targets:
        return pd.DataFrame()
    events = load_clinical_cofactor_events(config) if events is None else events.copy()
    groups = {
        "built_broad_cases": _case_index_frame(config, built_df),
        "matched_cases": _case_index_frame(config, matched_df),
    }
    if events.empty:
        rows = []
        for group, cases in groups.items():
            n_cases = int(cases["person_id"].nunique())
            for cofactor in targets:
                rows.append(_empty_prior_timing_row(group, cofactor, "any_prior", n_cases, 0, 0))
                for window, _, _ in _PRIOR_TIMING_BINS:
                    rows.append(_empty_prior_timing_row(group, cofactor, window, n_cases, 0, 0))
        return pd.DataFrame(rows)
    events = events.copy()
    events["person_id"] = events["person_id"].astype(str)
    events["cofactor"] = events["cofactor"].astype(str)
    events["condition_date"] = pd.to_datetime(events["condition_date"], errors="coerce")
    events = events[events["cofactor"].isin(targets)].dropna(subset=["condition_date"]).copy()
    rows: list[dict[str, object]] = []
    for group_name, cases in groups.items():
        n_cases = int(cases["person_id"].nunique()) if "person_id" in cases.columns else 0
        for cofactor in targets:
            subset = events[events["cofactor"] == cofactor].merge(cases, on="person_id", how="inner")
            if subset.empty:
                rows.append(_empty_prior_timing_row(group_name, cofactor, "any_prior", n_cases, 0, 0))
                for window, _, _ in _PRIOR_TIMING_BINS:
                    rows.append(_empty_prior_timing_row(group_name, cofactor, window, n_cases, 0, 0))
                continue
            subset["days_before_index"] = (subset["index_date"] - subset["condition_date"]).dt.days
            prior = subset[subset["days_before_index"] >= 0].copy()
            nearest = (
                prior.sort_values(["person_id", "days_before_index"])
                .drop_duplicates("person_id", keep="first")
                .copy()
            )
            days = pd.to_numeric(nearest.get("days_before_index"), errors="coerce").dropna()
            n_prior = int(nearest["person_id"].nunique()) if not nearest.empty else 0
            any_row = _empty_prior_timing_row(group_name, cofactor, "any_prior", n_cases, n_prior, n_prior)
            any_row["median_days_before_index"] = _quantile_text(days, 0.5)
            any_row["q1_days_before_index"] = _quantile_text(days, 0.25)
            any_row["q3_days_before_index"] = _quantile_text(days, 0.75)
            rows.append(any_row)
            for window, low, high in _PRIOR_TIMING_BINS:
                mask = days >= low
                if high is not None:
                    mask &= days <= high
                count = int(mask.sum())
                rows.append(_empty_prior_timing_row(group_name, cofactor, window, n_cases, n_prior, count))
    return pd.DataFrame(rows)


def write_case_cofactor_prior_timing_histogram(summary: pd.DataFrame, path: str) -> None:
    windows = [window for window, _, _ in _PRIOR_TIMING_BINS]
    if summary.empty:
        write_text(_prior_timing_placeholder_svg("No timing rows available."), path)
        return
    plot = summary[summary["window"].isin(windows)].copy()
    if plot.empty:
        write_text(_prior_timing_placeholder_svg("No binned prior timing rows available."), path)
        return
    plot["n_cases_in_window"] = pd.to_numeric(plot["n_cases_in_window"], errors="coerce").fillna(0)
    plot["pct_of_all_cases"] = pd.to_numeric(plot["pct_of_all_cases"], errors="coerce")
    groups = [group for group in ("matched_cases", "built_broad_cases") if group in set(plot["group"])]
    groups.extend(sorted(set(plot["group"]).difference(groups)))
    cofactors = [cofactor for cofactor in ("sepsis", "renal_injury") if cofactor in set(plot["cofactor"])]
    cofactors.extend(sorted(set(plot["cofactor"]).difference(cofactors)))
    if not groups or not cofactors:
        write_text(_prior_timing_placeholder_svg("No sepsis or renal-injury timing rows available."), path)
        return

    width = 1120
    panel_height = 310
    top_margin = 82
    bottom_margin = 66
    left = 86
    right = 36
    inner_width = width - left - right
    height = top_margin + panel_height * len(groups) + bottom_margin
    axis_color = "#1f2937"
    colors = {"sepsis": "#0f766e", "renal_injury": "#b45309"}
    fallback_colors = ["#2563eb", "#7c3aed", "#be123c", "#4b5563"]

    legend_parts = []
    legend_x = left
    for index, cofactor in enumerate(cofactors):
        color = colors.get(cofactor, fallback_colors[index % len(fallback_colors)])
        x = legend_x + index * 170
        legend_parts.append(
            f'<rect x="{x}" y="48" width="18" height="18" fill="{color}" />'
            f'<text x="{x + 26}" y="62" font-size="14" font-family="Menlo, monospace">{_svg_text(_prior_cofactor_label(cofactor))}</text>'
        )

    panels: list[str] = []
    for group_index, group in enumerate(groups):
        panel_top = top_margin + group_index * panel_height
        axis_y = panel_top + panel_height - 58
        plot_top = panel_top + 34
        plot_height = axis_y - plot_top
        group_df = plot[plot["group"] == group].copy()
        max_y = max(float(group_df["n_cases_in_window"].max()), 1.0)
        max_y *= 1.08
        bin_width = inner_width / len(windows)
        bar_width = min(30.0, (bin_width * 0.72) / max(len(cofactors), 1))
        panel_parts = [
            f'<text x="{left}" y="{panel_top + 18}" font-size="18" font-weight="700" font-family="Menlo, monospace">{_svg_text(_prior_group_label(group))}</text>',
            f'<line x1="{left}" y1="{axis_y}" x2="{width - right}" y2="{axis_y}" stroke="{axis_color}" />',
            f'<line x1="{left}" y1="{plot_top}" x2="{left}" y2="{axis_y}" stroke="{axis_color}" />',
            f'<text x="{left - 10}" y="{plot_top + 4}" text-anchor="end" font-size="11" font-family="Menlo, monospace">{int(round(max_y))}</text>',
            f'<text x="{left - 10}" y="{axis_y + 4}" text-anchor="end" font-size="11" font-family="Menlo, monospace">0</text>',
        ]
        for window_index, window in enumerate(windows):
            center = left + window_index * bin_width + bin_width / 2
            panel_parts.append(
                f'<text x="{center:.1f}" y="{axis_y + 22}" text-anchor="middle" font-size="12" font-family="Menlo, monospace">{_svg_text(_prior_window_label(window))}</text>'
            )
            for cofactor_index, cofactor in enumerate(cofactors):
                row = group_df[(group_df["window"] == window) & (group_df["cofactor"] == cofactor)]
                count = float(row["n_cases_in_window"].iloc[0]) if not row.empty else 0.0
                pct = row["pct_of_all_cases"].iloc[0] if not row.empty else np.nan
                bar_height = 0.0 if max_y <= 0 else (count / max_y) * plot_height
                x = center - (bar_width * len(cofactors)) / 2 + cofactor_index * bar_width
                y = axis_y - bar_height
                color = colors.get(cofactor, fallback_colors[cofactor_index % len(fallback_colors)])
                tooltip = (
                    f"{_prior_group_label(group)}; {_prior_cofactor_label(cofactor)}; "
                    f"{_prior_window_label(window)} days before index: {int(count)} cases"
                )
                if pd.notna(pct):
                    tooltip += f" ({float(pct):.1f}% of all cases)"
                panel_parts.append(
                    f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width - 2:.1f}" height="{bar_height:.1f}" fill="{color}" opacity="0.88">'
                    f"<title>{_svg_text(tooltip)}</title>"
                    "</rect>"
                )
                if count > 0:
                    label_y = max(plot_top + 12, y - 5)
                    panel_parts.append(
                        f'<text x="{x + (bar_width - 2) / 2:.1f}" y="{label_y:.1f}" text-anchor="middle" font-size="10" font-family="Menlo, monospace">{int(count)}</text>'
                    )
        panel_parts.append(
            f'<text x="{left + inner_width / 2:.1f}" y="{axis_y + 46}" text-anchor="middle" font-size="12" font-family="Menlo, monospace">Days before first rhabdomyolysis condition record</text>'
        )
        panel_parts.append(
            f'<text x="24" y="{plot_top + plot_height / 2:.1f}" font-size="12" transform="rotate(-90 24 {plot_top + plot_height / 2:.1f})" font-family="Menlo, monospace">Cases</text>'
        )
        panels.append("".join(panel_parts))

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
        '<rect width="100%" height="100%" fill="#fffdf7" />'
        f'<text x="{left}" y="30" font-size="24" font-family="Menlo, monospace">Prior sepsis/renal injury timing before rhabdo index</text>'
        f'{"".join(legend_parts)}'
        f'{"".join(panels)}'
        "</svg>"
    )
    write_text(svg, path)


def add_model_splits(
    matched_df: pd.DataFrame,
    *,
    train_fraction: float = _TRAIN_FRACTION,
    cv_folds: int = _CV_FOLDS,
) -> pd.DataFrame:
    output = matched_df.copy()
    if output.empty:
        output["analysis_split"] = pd.Series(dtype=object)
        output["cv_fold"] = pd.Series(dtype=object)
        return output
    group_column = "match_group_id" if "match_group_id" in output.columns else "person_id"
    output["_split_group"] = output[group_column].fillna(output["person_id"]).astype(str)
    groups = (
        pd.DataFrame({"_split_group": output["_split_group"].drop_duplicates().tolist()})
        .assign(_split_score=lambda frame: frame["_split_group"].map(_stable_hash_score))
        .sort_values(["_split_score", "_split_group"])
        .reset_index(drop=True)
    )
    if len(groups) <= 1:
        groups["analysis_split"] = "train"
    else:
        test_fraction = max(0.0, min(1.0, 1.0 - train_fraction))
        test_count = max(1, int(round(len(groups) * test_fraction)))
        test_count = min(test_count, len(groups) - 1)
        groups["analysis_split"] = "train"
        groups.loc[: test_count - 1, "analysis_split"] = "test"
    train_order = groups[groups["analysis_split"] == "train"].reset_index(drop=True)
    if not train_order.empty and cv_folds > 0:
        fold_count = min(cv_folds, len(train_order))
        train_order["cv_fold"] = (train_order.index % fold_count) + 1
        groups = groups.merge(train_order[["_split_group", "cv_fold"]], on="_split_group", how="left")
    else:
        groups["cv_fold"] = pd.NA
    output = output.merge(groups[["_split_group", "analysis_split", "cv_fold"]], on="_split_group", how="left")
    output = output.drop(columns=["_split_group"])
    output["cv_fold"] = output["cv_fold"].astype("Int64")
    return output


def build_model_split_summary(config: ProjectConfig, split_df: pd.DataFrame) -> pd.DataFrame:
    outcome = config.analysis.matched_outcome_column
    rows: list[dict[str, object]] = []

    def add_group(label: str, frame: pd.DataFrame) -> None:
        if frame.empty:
            rows.append(
                {
                    "group": label,
                    "n_rows": 0,
                    "n_match_groups": 0,
                    "n_cases": 0,
                    "n_controls": 0,
                    "controls_per_case": "",
                }
            )
            return
        cases = int((pd.to_numeric(frame.get(outcome), errors="coerce").fillna(0) == 1).sum())
        controls = int((pd.to_numeric(frame.get(outcome), errors="coerce").fillna(0) == 0).sum())
        group_column = "match_group_id" if "match_group_id" in frame.columns else "person_id"
        rows.append(
            {
                "group": label,
                "n_rows": int(len(frame)),
                "n_match_groups": int(frame[group_column].astype(str).nunique()),
                "n_cases": cases,
                "n_controls": controls,
                "controls_per_case": f"{controls / cases:.2f}" if cases else "",
            }
        )

    add_group("all_matched", split_df)
    if "analysis_split" in split_df.columns:
        for split in ("train", "test"):
            add_group(split, split_df[split_df["analysis_split"] == split])
    if "cv_fold" in split_df.columns:
        train = split_df[split_df.get("analysis_split", pd.Series(index=split_df.index, dtype=object)) == "train"]
        for fold in sorted(train["cv_fold"].dropna().astype(int).unique().tolist()):
            add_group(f"train_cv_fold_{fold}", train[train["cv_fold"].astype("Int64") == fold])
    return pd.DataFrame(rows)


def _flag_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(0, index=frame.index, dtype=int)
    return (pd.to_numeric(frame[column], errors="coerce").fillna(0) >= 1).astype(int)


def _any_flag(frame: pd.DataFrame, columns: tuple[str, ...]) -> pd.Series:
    if not columns:
        return pd.Series(0, index=frame.index, dtype=int)
    values = pd.DataFrame({column: _flag_series(frame, column) for column in columns}, index=frame.index)
    return (values.max(axis=1) >= 1).astype(int)


def _matched_case_group_flag(
    config: ProjectConfig,
    frame: pd.DataFrame,
    columns: tuple[str, ...],
) -> pd.Series:
    if frame.empty:
        return pd.Series(dtype=int)
    group_column = "match_group_id" if "match_group_id" in frame.columns else "person_id"
    outcome = config.analysis.matched_outcome_column
    outcome_values = pd.to_numeric(frame.get(outcome), errors="coerce").fillna(0)
    case_rows = frame[outcome_values == 1].copy()
    if case_rows.empty:
        return pd.Series(0, index=frame.index, dtype=int)
    case_rows["_group_flag"] = _any_flag(case_rows, columns)
    flags_by_group = case_rows.groupby(case_rows[group_column].astype(str))["_group_flag"].max()
    return frame[group_column].astype(str).map(flags_by_group).fillna(0).astype(int)


def add_model_eligibility_flags(config: ProjectConfig, split_df: pd.DataFrame) -> pd.DataFrame:
    output = split_df.copy()
    crush_columns = ("preindex_crush_injury", "periindex_crush_injury")
    output["row_preperi_crush_injury"] = _any_flag(output, crush_columns)
    output["matched_case_preperi_crush_injury"] = _matched_case_group_flag(config, output, crush_columns)
    output["row_preindex_sepsis"] = _flag_series(output, "preindex_sepsis")
    output["row_periindex_sepsis"] = _flag_series(output, "periindex_sepsis")
    output["matched_case_preindex_sepsis"] = _matched_case_group_flag(config, output, ("preindex_sepsis",))
    output["matched_case_periindex_sepsis"] = _matched_case_group_flag(config, output, ("periindex_sepsis",))

    required_numeric = ["age_at_index", "observation_days", "omop_condition_record_dates"]
    numeric_complete = pd.Series(True, index=output.index)
    for column in required_numeric:
        if column in output.columns:
            numeric_complete &= pd.to_numeric(output[column], errors="coerce").notna()
        else:
            numeric_complete &= False
    output["clinical_model_complete_core"] = numeric_complete.astype(int)

    nontraumatic = (
        (output["row_preperi_crush_injury"] == 0)
        & (output["matched_case_preperi_crush_injury"] == 0)
        & (output["clinical_model_complete_core"] == 1)
    )
    output["nontraumatic_model_eligible"] = nontraumatic.astype(int)
    output["primary_model_eligible"] = output["nontraumatic_model_eligible"]
    output["sensitivity_no_periindex_sepsis_eligible"] = (
        nontraumatic
        & (output["row_periindex_sepsis"] == 0)
        & (output["matched_case_periindex_sepsis"] == 0)
    ).astype(int)
    output["sensitivity_include_trauma_eligible"] = output["clinical_model_complete_core"]
    return output


def build_model_eligibility_summary(config: ProjectConfig, model_df: pd.DataFrame) -> pd.DataFrame:
    outcome = config.analysis.matched_outcome_column
    group_column = "match_group_id" if "match_group_id" in model_df.columns else "person_id"
    sets = [
        ("primary_model_eligible", "primary_nontraumatic"),
        ("sensitivity_no_periindex_sepsis_eligible", "sensitivity_no_periindex_sepsis"),
        ("sensitivity_include_trauma_eligible", "sensitivity_include_trauma"),
    ]
    splits: list[tuple[str, pd.DataFrame]] = [("all", model_df)]
    if "analysis_split" in model_df.columns:
        splits.extend(
            [
                ("train", model_df[model_df["analysis_split"] == "train"]),
                ("test", model_df[model_df["analysis_split"] == "test"]),
            ]
        )
    rows: list[dict[str, object]] = []
    for split_label, split in splits:
        split_outcome = pd.to_numeric(split.get(outcome), errors="coerce").fillna(0)
        split_cases = int((split_outcome == 1).sum())
        split_controls = int((split_outcome == 0).sum())
        for flag, label in sets:
            if flag not in split.columns:
                continue
            eligible = split[pd.to_numeric(split[flag], errors="coerce").fillna(0) == 1].copy()
            eligible_outcome = pd.to_numeric(eligible.get(outcome), errors="coerce").fillna(0)
            cases = int((eligible_outcome == 1).sum())
            controls = int((eligible_outcome == 0).sum())
            rows.append(
                {
                    "eligibility_set": label,
                    "split": split_label,
                    "n_rows": int(len(eligible)),
                    "n_match_groups": int(eligible[group_column].astype(str).nunique()) if not eligible.empty else 0,
                    "n_cases": cases,
                    "n_controls": controls,
                    "controls_per_case": f"{controls / cases:.2f}" if cases else "",
                    "excluded_rows": int(len(split) - len(eligible)),
                    "excluded_cases": int(split_cases - cases),
                    "excluded_controls": int(split_controls - controls),
                }
            )
    return pd.DataFrame(rows)


def _missingness_variables(frame: pd.DataFrame) -> list[str]:
    preferred = [
        "person_id",
        "analysis_case",
        "analysis_split",
        "cv_fold",
        "primary_model_eligible",
        "nontraumatic_model_eligible",
        "sensitivity_no_periindex_sepsis_eligible",
        "sensitivity_include_trauma_eligible",
        "clinical_model_complete_core",
        "eligible_control",
        "eligible_ehr_denominator",
        "matched_case_preperi_crush_injury",
        "row_preperi_crush_injury",
        "matched_case_preindex_sepsis",
        "matched_case_periindex_sepsis",
        "row_preindex_sepsis",
        "row_periindex_sepsis",
        "match_group_id",
        "index_date",
        "age_at_index",
        "observation_days",
        "omop_condition_record_dates",
        "gender_concept_name",
        "sex_category",
        "is_female",
        "ancestry_pred",
        *[f"pc{i}" for i in range(1, 11)],
        "preindex_crush_injury",
        "remote_preindex_sepsis",
        "near_preindex_sepsis",
        "preindex_sepsis",
        "remote_preindex_renal_injury",
        "near_preindex_renal_injury",
        "preindex_renal_injury",
        "periindex_crush_injury",
        "periindex_sepsis",
        "periindex_renal_injury",
        "postindex_sepsis",
        "postindex_renal_injury",
    ]
    return [column for column in preferred if column in frame.columns]


def build_missingness_summary(config: ProjectConfig, split_df: pd.DataFrame) -> pd.DataFrame:
    outcome = config.analysis.matched_outcome_column
    groups: list[tuple[str, pd.DataFrame]] = [("all_matched", split_df)]
    if outcome in split_df.columns:
        outcome_values = pd.to_numeric(split_df[outcome], errors="coerce")
        groups.extend(
            [
                ("matched_cases", split_df[outcome_values == 1]),
                ("matched_controls", split_df[outcome_values == 0]),
            ]
        )
    if "analysis_split" in split_df.columns:
        groups.extend(
            [
                ("train", split_df[split_df["analysis_split"] == "train"]),
                ("test", split_df[split_df["analysis_split"] == "test"]),
            ]
        )
    rows: list[dict[str, object]] = []
    for variable in _missingness_variables(split_df):
        for group_name, frame in groups:
            n = int(len(frame))
            missing = int(_missing_mask(frame, variable).sum()) if n else 0
            rows.append(
                {
                    "variable": variable,
                    "group": group_name,
                    "n": n,
                    "observed_n": n - missing,
                    "missing_n": missing,
                    "missing_pct": f"{(100.0 * missing / n):.1f}" if n else "",
                }
            )
    return pd.DataFrame(rows)


def build_split_table1(config: ProjectConfig, split_df: pd.DataFrame) -> pd.DataFrame:
    if "analysis_split" not in split_df.columns:
        return pd.DataFrame()
    tables = []
    for split in ("train", "test"):
        subset = split_df[split_df["analysis_split"] == split].copy()
        if subset.empty:
            continue
        table = build_matched_table1(config, subset)
        table.insert(0, "analysis_split", split)
        tables.append(table)
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


def build_clinical_model_input(config: ProjectConfig, split_df: pd.DataFrame) -> pd.DataFrame:
    output = split_df.copy()
    if "primary_model_eligible" not in output.columns:
        output = add_model_eligibility_flags(config, output)
    outcome = pd.to_numeric(
        output.get(config.analysis.matched_outcome_column, pd.Series(np.nan, index=output.index)),
        errors="coerce",
    )
    tier = output.get("case_tier", pd.Series("", index=output.index)).astype("string").str.lower()
    if "eligible_control" not in output.columns:
        output["eligible_control"] = (outcome.eq(0) & tier.eq("control")).astype(int)
    if "eligible_ehr_denominator" not in output.columns:
        condition_dates = pd.to_numeric(output.get("omop_condition_record_dates", 0), errors="coerce").fillna(0)
        output["eligible_ehr_denominator"] = condition_dates.ge(2).astype(int)
    if "broad_rhabdo_case" not in output.columns:
        output["broad_rhabdo_case"] = (outcome.eq(1) & tier.isin({"broad", "definite"})).astype(int)
    if "definite_rhabdo_case" not in output.columns:
        output["definite_rhabdo_case"] = (outcome.eq(1) & tier.eq("definite")).astype(int)
    if "high_ck_without_rhabdo" not in output.columns:
        output["high_ck_without_rhabdo"] = tier.eq("indeterminate_ck_only").astype(int)
    if _has_sex_information(output):
        output["model_sex_category"] = _sex_category_series(output)
    if "ancestry_pred" in output.columns:
        output["model_ancestry_pred"] = _category_series_with_missing(output, "ancestry_pred")
    preferred = [
        "person_id",
        "analysis_case",
        "analysis_split",
        "cv_fold",
        "primary_model_eligible",
        "nontraumatic_model_eligible",
        "sensitivity_no_periindex_sepsis_eligible",
        "sensitivity_include_trauma_eligible",
        "clinical_model_complete_core",
        "eligible_control",
        "eligible_ehr_denominator",
        "matched_case_preperi_crush_injury",
        "row_preperi_crush_injury",
        "matched_case_preindex_sepsis",
        "matched_case_periindex_sepsis",
        "row_preindex_sepsis",
        "row_periindex_sepsis",
        "match_group_id",
        "match_role",
        "matched_case_person_id",
        "case_tier",
        "broad_rhabdo_case",
        "definite_rhabdo_case",
        "high_ck_without_rhabdo",
        "index_date",
        "age_at_index",
        "observation_days",
        "omop_condition_record_dates",
        "model_sex_category",
        "sex_category",
        "is_female",
        "model_ancestry_pred",
        "ancestry_pred",
        *[f"pc{i}" for i in range(1, 11)],
        "preindex_crush_injury",
        "remote_preindex_sepsis",
        "near_preindex_sepsis",
        "preindex_sepsis",
        "remote_preindex_renal_injury",
        "near_preindex_renal_injury",
        "preindex_renal_injury",
        "periindex_crush_injury",
        "periindex_sepsis",
        "periindex_renal_injury",
        "postindex_sepsis",
        "postindex_renal_injury",
    ]
    columns = [column for column in preferred if column in output.columns]
    return output[columns].copy()


def characterize_case_control_cohort(
    config: ProjectConfig,
    built_df: pd.DataFrame,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
) -> dict[str, pd.DataFrame]:
    split_df = add_model_eligibility_flags(config, add_model_splits(matched_df))
    consort = build_consort_counts(built_df, matched_df)
    table1 = build_matched_table1(config, split_df)
    split_table1 = build_split_table1(config, split_df)
    critical = build_critical_illness_summary(config, built_df, split_df)
    prior_timing = build_case_cofactor_prior_timing_summary(config, built_df, split_df)
    split_summary = build_model_split_summary(config, split_df)
    eligibility = build_model_eligibility_summary(config, split_df)
    missingness = build_missingness_summary(config, split_df)
    model_input = build_clinical_model_input(config, split_df)

    write_dataframe(consort, consort_counts_path(paths))
    write_dataframe(table1, matched_table1_path(paths))
    write_dataframe(split_table1, split_table1_path(paths))
    write_dataframe(critical, critical_illness_summary_path(paths))
    write_dataframe(prior_timing, case_cofactor_prior_timing_path(paths))
    write_case_cofactor_prior_timing_histogram(prior_timing, case_cofactor_prior_timing_histogram_path(paths))
    write_dataframe(split_summary, model_split_summary_path(paths))
    write_dataframe(eligibility, model_eligibility_summary_path(paths))
    write_dataframe(missingness, missingness_summary_path(paths))
    write_dataframe(model_input, clinical_model_input_path(paths))
    _write_markdown_table(consort, consort_counts_report_path(paths), "CONSORT Counts")
    _write_markdown_table(table1, matched_table1_report_path(paths), "Matched Clinical Table 1")
    _write_markdown_table(split_table1, split_table1_report_path(paths), "Matched Clinical Table 1 by Split")
    _write_markdown_table(critical, critical_illness_summary_report_path(paths), "Sepsis and Renal Injury Timing Summary")
    write_text(
        "# Case Cofactor Prior Timing\n\n"
        "![Case cofactor prior timing histogram](case_cofactor_prior_timing_histogram.svg)\n\n"
        f"{_markdown_table(prior_timing)}\n",
        case_cofactor_prior_timing_report_path(paths),
    )
    _write_markdown_table(split_summary, model_split_summary_report_path(paths), "Model Split Summary")
    _write_markdown_table(eligibility, model_eligibility_summary_report_path(paths), "Model Eligibility Summary")
    _write_markdown_table(missingness, missingness_summary_report_path(paths), "Matched Cohort Missingness Summary")
    report = "\n\n".join(
        [
            "# Clinical Cohort Characterization",
            "## CONSORT Counts",
            _markdown_table(consort),
            "## Matched Clinical Table 1",
            _markdown_table(table1),
            "## Train/Test Split Summary",
            _markdown_table(split_summary),
            "## Model Eligibility Summary",
            _markdown_table(eligibility),
            "## Matched Clinical Table 1 by Split",
            _markdown_table(split_table1),
            "## Sepsis and Renal Injury Timing Summary",
            _markdown_table(critical),
            "## Case Cofactor Prior Timing",
            "![Case cofactor prior timing histogram](case_cofactor_prior_timing_histogram.svg)",
            _markdown_table(prior_timing),
            "## Missingness Summary",
            _markdown_table(missingness),
            "",
        ]
    )
    write_text(report, clinical_characterization_report_path(paths))
    return {
        "consort": consort,
        "table1": table1,
        "split_table1": split_table1,
        "critical_illness": critical,
        "case_cofactor_prior_timing": prior_timing,
        "split_summary": split_summary,
        "eligibility": eligibility,
        "missingness": missingness,
        "clinical_model_input": model_input,
    }


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

    if _has_sex_information(built_wgs) or _has_sex_information(matched_wgs):
        unmatched_sex_left = _sex_category_series(groups[_UNMATCHED_CASES])
        unmatched_sex_right = _sex_category_series(groups[_UNMATCHED_CONTROLS])
        matched_sex_left = _sex_category_series(groups[_MATCHED_CASES])
        matched_sex_right = _sex_category_series(groups[_MATCHED_CONTROLS])
        sex_required = ["female", "male", "other_or_unknown", "missing"]
        sex_values = _ensure_categories(
            _ordered_categories(
                [unmatched_sex_left, unmatched_sex_right, matched_sex_left, matched_sex_right],
                sex_required,
            ),
            sex_required,
        )
        add_row(
            "Demographics",
            "Sex category distribution, overall",
            lambda _frame: _empty_value(),
            unmatched_metrics=(
                _categorical_smd(unmatched_sex_left, unmatched_sex_right, sex_values),
                _categorical_p_value(unmatched_sex_left, unmatched_sex_right, sex_values),
            ),
            matched_metrics=(
                _categorical_smd(matched_sex_left, matched_sex_right, sex_values),
                _categorical_p_value(matched_sex_left, matched_sex_right, sex_values),
            ),
        )
        for sex_category in sex_values:
            add_row(
                "Demographics",
                f"Sex: {_display_category(sex_category)}, n (%)",
                lambda frame, sex_category=sex_category: _format_category_from_series(
                    _sex_category_series(frame),
                    sex_category,
                ),
                unmatched_metrics=(
                    _binary_smd(
                        (unmatched_sex_left == sex_category).astype(float),
                        (unmatched_sex_right == sex_category).astype(float),
                    ),
                    np.nan,
                ),
                matched_metrics=(
                    _binary_smd(
                        (matched_sex_left == sex_category).astype(float),
                        (matched_sex_right == sex_category).astype(float),
                    ),
                    np.nan,
                ),
            )

    unmatched_ancestry_left = _category_series_with_missing(groups[_UNMATCHED_CASES], "ancestry_pred")
    unmatched_ancestry_right = _category_series_with_missing(groups[_UNMATCHED_CONTROLS], "ancestry_pred")
    matched_ancestry_left = _category_series_with_missing(groups[_MATCHED_CASES], "ancestry_pred")
    matched_ancestry_right = _category_series_with_missing(groups[_MATCHED_CONTROLS], "ancestry_pred")
    ancestry_values = _ordered_categories(
        [unmatched_ancestry_left, unmatched_ancestry_right, matched_ancestry_left, matched_ancestry_right],
        ["AFR", "AMR", "EAS", "EUR", "MID", "SAS", "other", "missing"],
    )
    ancestry_values = _ensure_categories(ancestry_values, ["missing"])
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
            lambda frame, ancestry=ancestry: _format_category_from_series(
                _category_series_with_missing(frame, "ancestry_pred"),
                ancestry,
            ),
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
    "add_model_eligibility_flags",
    "add_model_splits",
    "build_case_cofactor_prior_timing_summary",
    "build_clinical_model_input",
    "build_model_eligibility_summary",
    "build_missingness_summary",
    "build_model_split_summary",
    "build_split_table1",
    "characterize_case_control_cohort",
    "case_cofactor_prior_timing_histogram_path",
    "case_cofactor_prior_timing_path",
    "case_cofactor_prior_timing_report_path",
    "clinical_model_input_path",
    "clinical_characterization_report_path",
    "consort_counts_path",
    "consort_counts_report_path",
    "cohort_summary_report_path",
    "cohort_summary_table_path",
    "critical_illness_summary_path",
    "critical_illness_summary_report_path",
    "matched_table1_path",
    "matched_table1_report_path",
    "missingness_summary_path",
    "missingness_summary_report_path",
    "model_eligibility_summary_path",
    "model_eligibility_summary_report_path",
    "model_split_summary_path",
    "model_split_summary_report_path",
    "split_table1_path",
    "split_table1_report_path",
    "summarize_clinical_demographics",
    "write_case_cofactor_prior_timing_histogram",
]
