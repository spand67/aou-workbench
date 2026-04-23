"""Clinical and demographic summary tables for unmatched and matched analyses."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import ProjectConfig
from .io_utils import read_table, write_dataframe, write_text
from .paths import ProjectPaths, join_path
from .stage1_prepare import stage1_sample_manifest_path


def cohort_summary_table_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "clinical_demographic_summary.tsv")


def cohort_summary_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "clinical_demographic_summary.md")


def _wgs_present_ids(config: ProjectConfig) -> set[str]:
    stage = config.analysis.stage1
    if stage is None:
        raise RuntimeError("Stage 1 must be configured to derive the WGS-present sample manifest.")
    manifest_path = stage1_sample_manifest_path(stage.variant_table)
    manifest = read_table(manifest_path)
    if "person_id" not in manifest.columns or manifest.empty:
        raise RuntimeError(
            f"Invalid WGS sample manifest at {manifest_path}. Rerun `aou-workbench prepare-stage1` first."
        )
    return set(manifest["person_id"].astype(str))


def _filter_to_wgs_present(frame: pd.DataFrame, present_ids: set[str]) -> pd.DataFrame:
    subset = frame.copy()
    subset["person_id"] = subset["person_id"].astype(str)
    return subset[subset["person_id"].isin(present_ids)].copy()


def _format_n(frame: pd.DataFrame) -> str:
    return str(int(frame["person_id"].astype(str).nunique()))


def _format_mean_sd(frame: pd.DataFrame, column: str) -> str:
    values = pd.to_numeric(frame.get(column), errors="coerce").dropna()
    if values.empty:
        return ""
    return f"{values.mean():.1f} ({values.std(ddof=1):.1f})"


def _format_binary(frame: pd.DataFrame, column: str) -> str:
    values = pd.to_numeric(frame.get(column), errors="coerce").dropna()
    if values.empty:
        return ""
    count = int((values >= 1).sum())
    return f"{count} ({(100.0 * count / len(values)):.1f}%)"


def _format_category(frame: pd.DataFrame, column: str, value: str) -> str:
    series = frame.get(column, pd.Series(dtype=object)).dropna().astype(str)
    if series.empty:
        return ""
    count = int((series == value).sum())
    return f"{count} ({(100.0 * count / len(series)):.1f}%)"


def _group_columns() -> list[str]:
    return [
        "unmatched_rhabdo_wgs",
        "unmatched_non_rhabdo_wgs",
        "matched_rhabdo_wgs",
        "matched_controls_wgs",
    ]


def _render_markdown(summary_df: pd.DataFrame) -> str:
    if summary_df.empty:
        return "# Clinical and Demographic Summary\n\n_No rows available._\n"
    try:
        table = summary_df.to_markdown(index=False)
    except ImportError:
        table = summary_df.to_string(index=False)
    return "\n".join(
        [
            "# Clinical and Demographic Summary",
            "",
            "Columns represent the two WGS-restricted analyses:",
            "- Unmatched OMOP rhabdo vs non-rhabdo",
            "- Matched OMOP rhabdo vs matched controls",
            "",
            table,
            "",
        ]
    )


def summarize_clinical_demographics(
    config: ProjectConfig,
    built_df: pd.DataFrame,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
) -> pd.DataFrame:
    present_ids = _wgs_present_ids(config)
    built_wgs = _filter_to_wgs_present(built_df, present_ids)
    matched_wgs = _filter_to_wgs_present(matched_df, present_ids)

    groups = {
        "unmatched_rhabdo_wgs": built_wgs[built_wgs["rhabdo_case"].fillna(0).astype(int) == 1].copy(),
        "unmatched_non_rhabdo_wgs": built_wgs[built_wgs["rhabdo_case"].fillna(0).astype(int) == 0].copy(),
        "matched_rhabdo_wgs": matched_wgs[matched_wgs[config.analysis.matched_outcome_column].fillna(0).astype(int) == 1].copy(),
        "matched_controls_wgs": matched_wgs[matched_wgs[config.analysis.matched_outcome_column].fillna(0).astype(int) == 0].copy(),
    }

    rows: list[dict[str, str]] = []

    def add_row(section: str, variable: str, formatter) -> None:
        row = {"section": section, "variable": variable}
        for group_name in _group_columns():
            row[group_name] = formatter(groups[group_name])
        rows.append(row)

    add_row("Sample size", "N", _format_n)

    age_column = "age_at_index" if "age_at_index" in built_wgs.columns else "age_raw"
    add_row("Demographics", "Age at index, mean (SD)", lambda frame: _format_mean_sd(frame, age_column))
    add_row("Demographics", "Female sex, n (%)", lambda frame: _format_binary(frame, "is_female"))

    ancestry_values = sorted(
        {
            str(value)
            for value in pd.concat(
                [
                    built_wgs.get("ancestry_pred", pd.Series(dtype=object)),
                    matched_wgs.get("ancestry_pred", pd.Series(dtype=object)),
                ],
                ignore_index=True,
            ).dropna()
            if str(value).strip()
        }
    )
    for ancestry in ancestry_values:
        add_row(
            "Demographics",
            f"Ancestry: {ancestry}, n (%)",
            lambda frame, ancestry=ancestry: _format_category(frame, "ancestry_pred", ancestry),
        )

    clinical_columns = [
        column
        for column in config.phenotype.clinical_cofactor_columns
        if column in built_wgs.columns or column in matched_wgs.columns
    ]
    for column in clinical_columns:
        label = column.replace("_", " ").title() + ", n (%)"
        add_row("Clinical", label, lambda frame, column=column: _format_binary(frame, column))

    summary = pd.DataFrame(rows)
    table_path = cohort_summary_table_path(paths)
    report_path = cohort_summary_report_path(paths)
    write_dataframe(summary, table_path)
    write_text(_render_markdown(summary), report_path)
    return summary


__all__ = [
    "cohort_summary_report_path",
    "cohort_summary_table_path",
    "summarize_clinical_demographics",
]
