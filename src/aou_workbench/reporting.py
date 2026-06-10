"""Markdown reporting helpers for stage and end-to-end summaries."""

from __future__ import annotations

import os
import html
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


def _format_card_value(value: object) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, float):
        if abs(value) < 1 and value != 0:
            return f"{value:.3f}"
        return f"{value:,.2f}".rstrip("0").rstrip(".")
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number.is_integer():
        return f"{int(number):,}"
    if abs(number) < 1 and number != 0:
        return f"{number:.3f}"
    return f"{number:,.2f}".rstrip("0").rstrip(".")


def _consort_value(consort: pd.DataFrame, step: str) -> str:
    if consort.empty or not {"step", "n"}.issubset(consort.columns):
        return ""
    match = consort[consort["step"] == step]
    if match.empty:
        return ""
    return _format_card_value(match["n"].iloc[0])


def _metric_value(metrics: pd.DataFrame, metric: str, evaluation_set: str = "test") -> str:
    if metrics.empty or metric not in metrics.columns or "evaluation_set" not in metrics.columns:
        return ""
    match = metrics[metrics["evaluation_set"].astype(str) == evaluation_set]
    if match.empty:
        return ""
    return _format_card_value(match[metric].iloc[0])


def _card(title: str, value: str, subtitle: str = "") -> str:
    value_html = html.escape(value or "Not generated")
    subtitle_html = f"<span>{html.escape(subtitle)}</span>" if subtitle else ""
    return f'<div class="metric-card"><span>{html.escape(title)}</span><strong>{value_html}</strong>{subtitle_html}</div>'


def _html_table(
    df: pd.DataFrame,
    *,
    columns: list[str] | None = None,
    limit: int | None = 40,
) -> str:
    if df.empty:
        return '<div class="empty">No rows available.</div>'
    view = df.copy()
    if columns is not None:
        selected = [column for column in columns if column in view.columns]
        if selected:
            view = view[selected]
    if limit is not None:
        view = view.head(limit)
    return view.to_html(index=False, escape=True, classes="data-table", border=0)


def _source_link(label: str, source_path: str | None, dashboard_path: str) -> str:
    if not source_path:
        return ""
    if _path_exists(source_path):
        rel = _relative_path(source_path, dashboard_path)
        return f'<a class="source-link" href="{html.escape(rel)}">{html.escape(label)}</a>'
    return f'<span class="missing-source">{html.escape(label)} not generated</span>'


def _dashboard_table(
    title: str,
    df: pd.DataFrame,
    *,
    dashboard_path: str,
    source_path: str | None = None,
    source_label: str = "Source table",
    columns: list[str] | None = None,
    limit: int | None = 40,
    note: str = "",
) -> str:
    source = _source_link(source_label, source_path, dashboard_path)
    note_html = f'<p class="section-note">{html.escape(note)}</p>' if note else ""
    return (
        '<article class="table-card">'
        f'<div class="table-head"><h3>{html.escape(title)}</h3>{source}</div>'
        f"{note_html}"
        '<div class="table-wrap">'
        f"{_html_table(df, columns=columns, limit=limit)}"
        "</div>"
        "</article>"
    )


def _dashboard_figure(title: str, figure_path: str | None, dashboard_path: str, caption: str = "") -> str:
    caption_html = f'<p class="section-note">{html.escape(caption)}</p>' if caption else ""
    if not figure_path or not _path_exists(figure_path):
        return (
            '<article class="figure-card">'
            f"<h3>{html.escape(title)}</h3>"
            f"{caption_html}"
            '<div class="empty">Figure not generated yet.</div>'
            "</article>"
        )
    rel = _relative_path(figure_path, dashboard_path)
    return (
        '<article class="figure-card">'
        f"<h3>{html.escape(title)}</h3>"
        f"{caption_html}"
        f'<img src="{html.escape(rel)}" alt="{html.escape(title)}" />'
        "</article>"
    )


def _section(section_id: str, title: str, body: str) -> str:
    return f'<section id="{html.escape(section_id)}" class="dashboard-section"><h2>{html.escape(title)}</h2>{body}</section>'


def _dashboard_css() -> str:
    return """
:root {
  color-scheme: light;
  --bg: #f6f7f9;
  --panel: #ffffff;
  --panel-soft: #f9fafb;
  --border: #d9dee7;
  --text: #172033;
  --muted: #5d6678;
  --teal: #0f766e;
  --amber: #b45309;
  --blue: #2563eb;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--bg);
  color: var(--text);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
}
.layout {
  display: grid;
  grid-template-columns: 238px 1fr;
  min-height: 100vh;
}
.sidebar {
  position: sticky;
  top: 0;
  height: 100vh;
  overflow: auto;
  padding: 24px 18px;
  background: #101827;
  color: #eef2f7;
}
.sidebar h1 {
  margin: 0 0 8px;
  font-size: 18px;
  line-height: 1.25;
}
.sidebar p {
  margin: 0 0 20px;
  color: #b8c2d2;
  font-size: 12px;
  line-height: 1.45;
}
.sidebar a {
  display: block;
  color: #e5edf7;
  text-decoration: none;
  padding: 8px 10px;
  border-radius: 6px;
  font-size: 14px;
}
.sidebar a:hover { background: #1f2a3d; }
main { padding: 28px; }
.hero {
  background: var(--panel);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 24px;
  margin-bottom: 18px;
}
.hero h1 {
  margin: 0 0 6px;
  font-size: 28px;
  line-height: 1.2;
}
.hero p {
  margin: 0;
  color: var(--muted);
  font-size: 14px;
}
.metrics-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  gap: 12px;
  margin: 14px 0 6px;
}
.metric-card {
  background: var(--panel);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 14px;
  min-height: 94px;
}
.metric-card span {
  display: block;
  color: var(--muted);
  font-size: 12px;
  line-height: 1.35;
}
.metric-card strong {
  display: block;
  margin-top: 8px;
  font-size: 28px;
  line-height: 1.05;
}
.metric-card strong + span { margin-top: 7px; }
.dashboard-section {
  margin: 22px 0;
}
.dashboard-section h2 {
  margin: 0 0 12px;
  font-size: 22px;
}
.cards-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
  gap: 16px;
}
.table-card,
.figure-card {
  background: var(--panel);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 16px;
  overflow: hidden;
}
.table-card h3,
.figure-card h3 {
  margin: 0;
  font-size: 16px;
}
.table-head {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 12px;
  margin-bottom: 10px;
}
.source-link,
.missing-source {
  color: var(--blue);
  font-size: 12px;
  white-space: nowrap;
}
.missing-source { color: var(--muted); }
.section-note {
  margin: 8px 0 12px;
  color: var(--muted);
  font-size: 13px;
  line-height: 1.45;
}
.table-wrap {
  max-height: 520px;
  overflow: auto;
  border: 1px solid var(--border);
  border-radius: 6px;
}
.data-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 12px;
}
.data-table th {
  position: sticky;
  top: 0;
  background: #edf2f7;
  color: #263346;
  text-align: left;
  font-weight: 650;
  z-index: 1;
}
.data-table th,
.data-table td {
  border-bottom: 1px solid #e7ebf1;
  padding: 7px 9px;
  vertical-align: top;
}
.data-table tr:nth-child(even) td { background: var(--panel-soft); }
.figure-card img {
  display: block;
  width: 100%;
  max-height: 620px;
  object-fit: contain;
  border: 1px solid var(--border);
  border-radius: 6px;
  background: #fffdf7;
}
.empty {
  color: var(--muted);
  background: var(--panel-soft);
  border: 1px dashed var(--border);
  border-radius: 6px;
  padding: 18px;
}
@media (max-width: 900px) {
  .layout { grid-template-columns: 1fr; }
  .sidebar {
    position: relative;
    height: auto;
  }
  main { padding: 16px; }
  .cards-grid { grid-template-columns: 1fr; }
}
"""


def write_dashboard_report(
    *,
    analysis_name: str,
    output_root: str,
    path: str,
    cohort_tables: dict[str, pd.DataFrame] | None = None,
    clinical_model_tables: dict[str, pd.DataFrame] | None = None,
    preindex_tables: dict[str, pd.DataFrame] | None = None,
    genetics_tables: dict[str, pd.DataFrame] | None = None,
    source_paths: dict[str, str] | None = None,
    figure_paths: dict[str, str] | None = None,
) -> None:
    cohort_tables = cohort_tables or {}
    clinical_model_tables = clinical_model_tables or {}
    preindex_tables = preindex_tables or {}
    genetics_tables = genetics_tables or {}
    source_paths = source_paths or {}
    figure_paths = figure_paths or {}
    consort = cohort_tables.get("consort", pd.DataFrame())
    metrics = clinical_model_tables.get("metrics", pd.DataFrame())

    cohort_cards = [
        _card("Built cohort rows", _consort_value(consort, "Built cohort rows")),
        _card("Denominator eligible", _consort_value(consort, ">=2 OMOP condition dates")),
        _card("Broad rhabdo cases", _consort_value(consort, "Broad rhabdo cases")),
        _card("Definite rhabdo cases", _consort_value(consort, "Definite rhabdo cases")),
        _card("Eligible controls", _consort_value(consort, "Eligible controls")),
        _card("Matched analytic rows", _consort_value(consort, "Matched analytic rows")),
        _card("Matched cases", _consort_value(consort, "Matched cases")),
        _card("Matched controls", _consort_value(consort, "Matched controls")),
    ]
    model_cards = [
        _card("Test ROC AUC", _metric_value(metrics, "roc_auc")),
        _card("Test average precision", _metric_value(metrics, "average_precision")),
        _card("Sensitivity", _metric_value(metrics, "sensitivity"), "selected training threshold"),
        _card("Specificity", _metric_value(metrics, "specificity"), "selected training threshold"),
        _card("PPV", _metric_value(metrics, "ppv"), "selected training threshold"),
        _card("NPV", _metric_value(metrics, "npv"), "selected training threshold"),
    ]

    sections = [
        _section(
            "overview",
            "Overview",
            '<div class="metrics-grid">' + "".join(cohort_cards) + "</div>"
            + '<div class="metrics-grid">' + "".join(model_cards) + "</div>",
        ),
        _section(
            "cohort",
            "Cohort And Case-Control Tables",
            '<div class="cards-grid">'
            + _dashboard_table(
                "CONSORT Counts",
                consort,
                dashboard_path=path,
                source_path=source_paths.get("consort"),
                source_label="consort_counts.tsv",
                limit=None,
            )
            + _dashboard_table(
                "Train/Test And CV Split",
                cohort_tables.get("split_summary", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("split_summary"),
                source_label="model_split_summary.tsv",
                limit=None,
            )
            + _dashboard_table(
                "Model Eligibility",
                cohort_tables.get("eligibility", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("eligibility"),
                source_label="model_eligibility_summary.tsv",
                limit=None,
            )
            + _dashboard_table(
                "Matched Clinical Table 1",
                cohort_tables.get("table1", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("table1"),
                source_label="table1_clinical_matched.tsv",
                limit=None,
            )
            + _dashboard_table(
                "Matched Clinical Table 1 By Split",
                cohort_tables.get("split_table1", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("split_table1"),
                source_label="table1_clinical_by_split.tsv",
                limit=None,
            )
            + "</div>",
        ),
        _section(
            "timing",
            "Sepsis And Renal-Injury Timing",
            '<div class="cards-grid">'
            + _dashboard_figure(
                "Case Cofactor Prior Timing Histogram",
                figure_paths.get("case_cofactor_prior_timing_histogram"),
                path,
                "Nearest prior sepsis and renal-injury condition records among rhabdomyolysis cases.",
            )
            + _dashboard_table(
                "Case Cofactor Prior Timing",
                cohort_tables.get("case_cofactor_prior_timing", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("case_cofactor_prior_timing"),
                source_label="case_cofactor_prior_timing.tsv",
                limit=None,
            )
            + _dashboard_table(
                "Critical Illness Timing Summary",
                cohort_tables.get("critical_illness", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("critical_illness"),
                source_label="critical_illness_summary.tsv",
                limit=None,
            )
            + "</div>",
        ),
        _section(
            "model",
            "Clinical-Only Model",
            '<div class="cards-grid">'
            + _dashboard_figure("ROC Curve", figure_paths.get("clinical_model_roc"), path)
            + _dashboard_figure("Precision-Recall Curve", figure_paths.get("clinical_model_pr"), path)
            + _dashboard_figure("Calibration Plot", figure_paths.get("clinical_model_calibration"), path)
            + _dashboard_table(
                "Clinical Model Metrics",
                metrics,
                dashboard_path=path,
                source_path=source_paths.get("clinical_model_metrics"),
                source_label="metrics.tsv",
                limit=None,
            )
            + _dashboard_table(
                "Clinical Model Cross-Validation Metrics",
                clinical_model_tables.get("cv_metrics", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("clinical_model_cv_metrics"),
                source_label="cv_metrics.tsv",
                limit=None,
            )
            + _dashboard_table(
                "Clinical Model Coefficients",
                clinical_model_tables.get("coefficients", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("clinical_model_coefficients"),
                source_label="coefficients.tsv",
                limit=40,
            )
            + "</div>",
        ),
        _section(
            "preindex",
            "Pre-Index Case Profile",
            '<div class="cards-grid">'
            + _dashboard_table(
                "Availability Summary",
                preindex_tables.get("summary", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("preindex_summary"),
                source_label="availability_summary.tsv",
                limit=None,
            )
            + _dashboard_table(
                "Biomarker Availability",
                preindex_tables.get("biomarkers", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("preindex_biomarkers"),
                source_label="biomarker_availability.tsv",
                limit=40,
            )
            + _dashboard_table(
                "Top Pre-Index Conditions",
                preindex_tables.get("top_conditions", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("preindex_top_conditions"),
                source_label="top_conditions.tsv",
                limit=25,
            )
            + _dashboard_table(
                "Top Pre-Index Measurements",
                preindex_tables.get("top_measurements", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("preindex_top_measurements"),
                source_label="top_measurements.tsv",
                limit=25,
            )
            + "</div>",
        ),
        _section(
            "missingness",
            "Missingness",
            _dashboard_table(
                "Variables With Missing Values",
                _filtered_missingness(cohort_tables.get("missingness", pd.DataFrame())),
                dashboard_path=path,
                source_path=source_paths.get("missingness"),
                source_label="missingness_summary.tsv",
                limit=40,
                note="Shown only for variables with missing values; see source table for the full audit.",
            ),
        ),
        _section(
            "genetics",
            "Genetics",
            '<div class="cards-grid">'
            + _dashboard_figure("GWAS Manhattan Plot", figure_paths.get("stage4_manhattan"), path)
            + _dashboard_figure("GWAS QQ Plot", figure_paths.get("stage4_qq"), path)
            + _dashboard_table(
                "Stage 1: A Priori Variants",
                genetics_tables.get("stage1", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("stage1"),
                source_label="prior_variant_results.tsv",
                columns=["label", "gene", "case_carriers", "control_carriers", "fisher_p"],
                limit=30,
            )
            + _dashboard_table(
                "Stage 2: P/LP Panel Summary",
                genetics_tables.get("stage2_genes", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("stage2_genes"),
                source_label="plp_gene_results.tsv",
                columns=["gene", "case_carriers", "control_carriers", "regression_p"],
                limit=30,
            )
            + _dashboard_table(
                "Stage 3: Burden Results",
                genetics_tables.get("stage3", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("stage3"),
                source_label="burden_results.tsv",
                columns=["gene", "mask", "case_carriers", "control_carriers", "regression_p"],
                limit=30,
            )
            + _dashboard_table(
                "Stage 4: GWAS Lead Hits",
                genetics_tables.get("stage4_hits", pd.DataFrame()),
                dashboard_path=path,
                source_path=source_paths.get("stage4_hits"),
                source_label="gwas_lead_hits.tsv",
                columns=["variant_id", "gene", "chromosome", "position", "regression_p"],
                limit=30,
            )
            + "</div>",
        ),
    ]
    nav = "".join(
        f'<a href="#{section_id}">{label}</a>'
        for section_id, label in [
            ("overview", "Overview"),
            ("cohort", "Cohort"),
            ("timing", "Timing"),
            ("model", "Model"),
            ("preindex", "Pre-index"),
            ("missingness", "Missingness"),
            ("genetics", "Genetics"),
        ]
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(analysis_name)} Rhabdomyolysis Dashboard</title>
  <style>{_dashboard_css()}</style>
</head>
<body>
  <div class="layout">
    <aside class="sidebar">
      <h1>{html.escape(analysis_name)}</h1>
      <p>Rhabdomyolysis analysis dashboard<br />Run root: {html.escape(output_root)}</p>
      {nav}
    </aside>
    <main>
      <div class="hero">
        <h1>Rhabdomyolysis Analysis Dashboard</h1>
        <p>Formatted aggregate outputs, figures, and model summaries generated from the current run directory.</p>
      </div>
      {''.join(sections)}
    </main>
  </div>
</body>
</html>
"""
    write_text(document, path)


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


__all__ = ["load_table_if_exists", "write_dashboard_report", "write_final_report", "write_stage_report"]
