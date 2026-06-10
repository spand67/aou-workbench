"""Clinical and demographic summary tables for unmatched and matched analyses."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, t, ttest_ind

from .config import ProjectConfig
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


def clinical_characterization_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "clinical_characterization_report.md")


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


def _normalize_category_values(series: pd.Series) -> pd.Series:
    normalized = series.astype(str).str.lower().str.strip()
    normalized = normalized.replace(
        {
            "": "missing",
            "nan": "missing",
            "none": "missing",
            "null": "missing",
            "na": "missing",
            "n/a": "missing",
            "no matching concept": "missing",
            "unknown": "missing",
            "skip": "missing",
            "pmi: skip": "missing",
            "prefer not to answer": "missing",
            "pmi: prefer not to answer": "missing",
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


def _format_category_from_series(series: pd.Series, value: str) -> str:
    if series.empty:
        return _empty_value()
    count = int((series == value).sum())
    return f"{count} ({(100.0 * count / len(series)):.1f}%)"


def _display_category(value: str) -> str:
    return value.replace("_or_", "/").replace("_", " ")


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
    if table.sum() == 0 or (table.sum(axis=0) == 0).any() or (table.sum(axis=1) == 0).any():
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
        sex_categories = _ordered_categories(
            [sex_cases, sex_controls],
            ["female", "male", "other_or_unknown", "missing"],
        )
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
        categories = sorted(matched_df[ancestry].dropna().astype(str).unique().tolist())
        left = _categorical_series(cases, ancestry)
        right = _categorical_series(controls, ancestry)
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
            temp_col = f"__ancestry_{category}"
            matched_df[temp_col] = (matched_df[ancestry].astype(str) == category).astype(int)
            cases[temp_col] = (cases[ancestry].astype(str) == category).astype(int)
            controls[temp_col] = (controls[ancestry].astype(str) == category).astype(int)
            add_binary("Demographics", f"Ancestry: {category}, n (%)", temp_col)

    clinical_columns = [
        "crush_injury",
        "preindex_crush_injury",
        "periindex_crush_injury",
        "sepsis",
        "preindex_sepsis",
        "periindex_sepsis",
        "postindex_sepsis",
        "renal_injury",
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
        "preindex_sepsis",
        "periindex_sepsis",
        "postindex_sepsis",
        "renal_injury",
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


def characterize_case_control_cohort(
    config: ProjectConfig,
    built_df: pd.DataFrame,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
) -> dict[str, pd.DataFrame]:
    consort = build_consort_counts(built_df, matched_df)
    table1 = build_matched_table1(config, matched_df)
    critical = build_critical_illness_summary(config, built_df, matched_df)

    write_dataframe(consort, consort_counts_path(paths))
    write_dataframe(table1, matched_table1_path(paths))
    write_dataframe(critical, critical_illness_summary_path(paths))
    _write_markdown_table(consort, consort_counts_report_path(paths), "CONSORT Counts")
    _write_markdown_table(table1, matched_table1_report_path(paths), "Matched Clinical Table 1")
    _write_markdown_table(critical, critical_illness_summary_report_path(paths), "Sepsis and Renal Injury Timing Summary")
    report = "\n\n".join(
        [
            "# Clinical Cohort Characterization",
            "## CONSORT Counts",
            _markdown_table(consort),
            "## Matched Clinical Table 1",
            _markdown_table(table1),
            "## Sepsis and Renal Injury Timing Summary",
            _markdown_table(critical),
            "",
        ]
    )
    write_text(report, clinical_characterization_report_path(paths))
    return {"consort": consort, "table1": table1, "critical_illness": critical}


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
        sex_values = _ordered_categories(
            [unmatched_sex_left, unmatched_sex_right, matched_sex_left, matched_sex_right],
            ["female", "male", "other_or_unknown", "missing"],
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
    "characterize_case_control_cohort",
    "clinical_characterization_report_path",
    "consort_counts_path",
    "consort_counts_report_path",
    "cohort_summary_report_path",
    "cohort_summary_table_path",
    "critical_illness_summary_path",
    "critical_illness_summary_report_path",
    "matched_table1_path",
    "matched_table1_report_path",
    "summarize_clinical_demographics",
]
