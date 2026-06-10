"""Markdown reporting helpers for stage and end-to-end summaries."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable

import pandas as pd

from .io_utils import write_text


def dataframe_markdown(df: pd.DataFrame, *, columns: Iterable[str] | None = None, limit: int | None = 10) -> str:
    if df.empty:
        return "_No rows available._"
    view = df.copy()
    if columns is not None:
        selected = [column for column in columns if column in view.columns]
        if selected:
            view = view[selected]
    if limit is not None:
        view = view.head(limit)
    try:
        return view.to_markdown(index=False)
    except ImportError:
        return view.to_string(index=False)


def write_stage_report(
    *,
    title: str,
    summary_lines: list[str],
    preview_df: pd.DataFrame,
    preview_columns: list[str],
    path: str,
) -> None:
    body = [f"# {title}", ""]
    body.extend(summary_lines)
    body.extend(["", "## Preview", "", dataframe_markdown(preview_df, columns=preview_columns)])
    write_text("\n".join(body), path)


def load_table_if_exists(path: str) -> pd.DataFrame:
    target = Path(path)
    if not target.exists():
        return pd.DataFrame()
    if path.endswith(".tsv") or path.endswith(".tsv.gz") or path.endswith(".bgz"):
        return pd.read_csv(path, sep="\t")
    if path.endswith(".csv") or path.endswith(".csv.gz"):
        return pd.read_csv(path)
    return pd.DataFrame()


def _path_exists(path: str | None) -> bool:
    return bool(path) and not str(path).startswith("gs://") and Path(str(path)).exists()


def _relative_path(target: str, report_path: str) -> str:
    if target.startswith("gs://"):
        return target
    try:
        return os.path.relpath(target, start=Path(report_path).parent)
    except ValueError:
        return target


def _source_line(label: str, source_path: str | None, report_path: str) -> str:
    if not source_path:
        return ""
    if _path_exists(source_path):
        rel = _relative_path(source_path, report_path)
        return f"Source: [{label}]({rel})"
    return f"Source: `{source_path}` not generated yet."


def _figure_block(title: str, figure_path: str | None, report_path: str, caption: str = "") -> list[str]:
    if not figure_path:
        return []
    if not _path_exists(figure_path):
        return [f"- {title}: `{figure_path}` not generated yet."]
    rel = _relative_path(figure_path, report_path)
    lines = [f"### {title}", ""]
    if caption:
        lines.extend([caption, ""])
    lines.append(f"![{title}]({rel})")
    return lines


def _table_section(
    title: str,
    df: pd.DataFrame,
    *,
    report_path: str,
    source_path: str | None = None,
    source_label: str = "table",
    columns: list[str] | None = None,
    limit: int | None = 30,
    note: str = "",
) -> list[str]:
    lines = [f"### {title}", ""]
    if note:
        lines.extend([note, ""])
    source = _source_line(source_label, source_path, report_path)
    if source:
        lines.extend([source, ""])
    lines.append(dataframe_markdown(df, columns=columns, limit=limit))
    return lines


def _filtered_missingness(missingness: pd.DataFrame) -> pd.DataFrame:
    if missingness.empty or "missing_n" not in missingness.columns:
        return missingness
    output = missingness.copy()
    output["missing_n"] = pd.to_numeric(output["missing_n"], errors="coerce").fillna(0)
    output = output[output["missing_n"] > 0].copy()
    if "missing_pct" in output.columns:
        output["_missing_pct_numeric"] = pd.to_numeric(output["missing_pct"], errors="coerce").fillna(0)
        output = output.sort_values(["_missing_pct_numeric", "missing_n"], ascending=False)
        output = output.drop(columns=["_missing_pct_numeric"])
    else:
        output = output.sort_values("missing_n", ascending=False)
    return output


def write_final_report(
    *,
    analysis_name: str,
    output_root: str,
    stage1: pd.DataFrame,
    stage2_genes: pd.DataFrame,
    stage3: pd.DataFrame,
    stage4_hits: pd.DataFrame,
    path: str,
    cohort_tables: dict[str, pd.DataFrame] | None = None,
    clinical_model_tables: dict[str, pd.DataFrame] | None = None,
    preindex_tables: dict[str, pd.DataFrame] | None = None,
    source_paths: dict[str, str] | None = None,
    figure_paths: dict[str, str] | None = None,
) -> None:
    cohort_tables = cohort_tables or {}
    clinical_model_tables = clinical_model_tables or {}
    preindex_tables = preindex_tables or {}
    source_paths = source_paths or {}
    figure_paths = figure_paths or {}

    sections = [
        f"# {analysis_name} Rhabdomyolysis Analysis Summary",
        "",
        f"Run root: `{output_root}`",
        "",
        "This report collects the aggregate tables and figures generated so far. Row-level cohort/model files are linked from command output but are not embedded here.",
        "",
        "## Cohort Characterization",
        "",
    ]
    sections.extend(
        _table_section(
            "CONSORT Counts",
            cohort_tables.get("consort", pd.DataFrame()),
            report_path=path,
            source_path=source_paths.get("consort"),
            source_label="consort_counts.tsv",
            limit=None,
        )
    )
    sections.extend(
        [
            "",
            *_table_section(
                "Train/Test and CV Split",
                cohort_tables.get("split_summary", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("split_summary"),
                source_label="model_split_summary.tsv",
                limit=None,
            ),
            "",
            *_table_section(
                "Model Eligibility",
                cohort_tables.get("eligibility", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("eligibility"),
                source_label="model_eligibility_summary.tsv",
                limit=None,
            ),
            "",
            *_table_section(
                "Matched Clinical Table 1",
                cohort_tables.get("table1", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("table1"),
                source_label="table1_clinical_matched.tsv",
                limit=None,
            ),
            "",
            *_table_section(
                "Matched Clinical Table 1 by Split",
                cohort_tables.get("split_table1", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("split_table1"),
                source_label="table1_clinical_by_split.tsv",
                limit=None,
            ),
            "",
            "## Sepsis and Renal-Injury Timing",
            "",
            *_figure_block(
                "Case Cofactor Prior Timing Histogram",
                figure_paths.get("case_cofactor_prior_timing_histogram"),
                path,
                "Nearest prior sepsis and renal-injury condition records among rhabdomyolysis cases.",
            ),
            "",
            *_table_section(
                "Case Cofactor Prior Timing",
                cohort_tables.get("case_cofactor_prior_timing", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("case_cofactor_prior_timing"),
                source_label="case_cofactor_prior_timing.tsv",
                limit=None,
            ),
            "",
            *_table_section(
                "Critical Illness Timing Summary",
                cohort_tables.get("critical_illness", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("critical_illness"),
                source_label="critical_illness_summary.tsv",
                limit=None,
            ),
            "",
            "## Missingness",
            "",
            *_table_section(
                "Variables With Missing Values",
                _filtered_missingness(cohort_tables.get("missingness", pd.DataFrame())),
                report_path=path,
                source_path=source_paths.get("missingness"),
                source_label="missingness_summary.tsv",
                limit=30,
                note="Shown only for variables with missing values; see the source table for the full matched-cohort missingness audit.",
            ),
            "",
            "## Pre-Index Case Profile",
            "",
            *_table_section(
                "Availability Summary",
                preindex_tables.get("summary", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("preindex_summary"),
                source_label="availability_summary.tsv",
                limit=None,
            ),
            "",
            *_table_section(
                "Biomarker Availability",
                preindex_tables.get("biomarkers", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("preindex_biomarkers"),
                source_label="biomarker_availability.tsv",
                limit=40,
            ),
            "",
            *_table_section(
                "Top Pre-Index Conditions",
                preindex_tables.get("top_conditions", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("preindex_top_conditions"),
                source_label="top_conditions.tsv",
                limit=25,
            ),
            "",
            *_table_section(
                "Top Pre-Index Measurements",
                preindex_tables.get("top_measurements", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("preindex_top_measurements"),
                source_label="top_measurements.tsv",
                limit=25,
            ),
            "",
            "## Clinical-Only Prediction Model",
            "",
            *_figure_block("ROC Curve", figure_paths.get("clinical_model_roc"), path),
            "",
            *_figure_block("Precision-Recall Curve", figure_paths.get("clinical_model_pr"), path),
            "",
            *_figure_block("Calibration Plot", figure_paths.get("clinical_model_calibration"), path),
            "",
            *_table_section(
                "Clinical Model Metrics",
                clinical_model_tables.get("metrics", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("clinical_model_metrics"),
                source_label="metrics.tsv",
                limit=None,
            ),
            "",
            *_table_section(
                "Clinical Model Cross-Validation Metrics",
                clinical_model_tables.get("cv_metrics", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("clinical_model_cv_metrics"),
                source_label="cv_metrics.tsv",
                limit=None,
            ),
            "",
            *_table_section(
                "Clinical Model Coefficients",
                clinical_model_tables.get("coefficients", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("clinical_model_coefficients"),
                source_label="coefficients.tsv",
                limit=40,
            ),
            "",
            *_table_section(
                "Clinical Model Calibration Table",
                clinical_model_tables.get("calibration", pd.DataFrame()),
                report_path=path,
                source_path=source_paths.get("clinical_model_calibration_table"),
                source_label="calibration.tsv",
                limit=None,
            ),
            "",
            "## Genetics",
            "",
            *_figure_block("GWAS Manhattan Plot", figure_paths.get("stage4_manhattan"), path),
            "",
            *_figure_block("GWAS QQ Plot", figure_paths.get("stage4_qq"), path),
            "",
            *_table_section(
                "Stage 1: A Priori Variants",
                stage1,
                report_path=path,
                source_path=source_paths.get("stage1"),
                source_label="prior_variant_results.tsv",
                columns=["label", "gene", "case_carriers", "control_carriers", "fisher_p"],
                limit=30,
            ),
            "",
            *_table_section(
                "Stage 2: P/LP Panel Summary",
                stage2_genes,
                report_path=path,
                source_path=source_paths.get("stage2_genes"),
                source_label="plp_gene_results.tsv",
                columns=["gene", "case_carriers", "control_carriers", "regression_p"],
                limit=30,
            ),
            "",
            *_table_section(
                "Stage 3: Burden Results",
                stage3,
                report_path=path,
                source_path=source_paths.get("stage3"),
                source_label="burden_results.tsv",
                columns=["gene", "mask", "case_carriers", "control_carriers", "regression_p"],
                limit=30,
            ),
            "",
            *_table_section(
                "Stage 4: GWAS Lead Hits",
                stage4_hits,
                report_path=path,
                source_path=source_paths.get("stage4_hits"),
                source_label="gwas_lead_hits.tsv",
                columns=["variant_id", "gene", "chromosome", "position", "regression_p"],
                limit=30,
            ),
        ]
    )
    write_text("\n".join(sections), path)


__all__ = ["load_table_if_exists", "write_final_report", "write_stage_report"]
