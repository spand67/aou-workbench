"""Supervisor-facing presentation dashboard for the rhabdomyolysis analysis."""

from __future__ import annotations

import html
import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .clinical_model import clinical_model_metrics_path
from .clinical_prs_model import (
    clinical_prs_model_coefficients_path,
    clinical_prs_model_dir,
    clinical_prs_model_metrics_path,
)
from .cohort import load_clinical_cofactor_events
from .cohort_summary import (
    clinical_model_input_path,
    consort_counts_path,
    model_split_summary_path,
    split_table1_path,
)
from .config import ProjectConfig
from .io_utils import ensure_parent_dir, read_table, slugify, write_dataframe, write_text
from .microarray_plink_gwas import (
    microarray_plink_manhattan_path,
    microarray_plink_qc_path,
    microarray_plink_results_path,
    microarray_plink_variant_qc_summary_path,
)
from .microarray_plink_prs import microarray_prs_metrics_path, microarray_prs_weights_path
from .paths import ProjectPaths, join_path


TIMING_BINS: tuple[tuple[str, int, int], ...] = (
    ("-365 to -91", -365, -91),
    ("-90 to -31", -90, -31),
    ("-30 to -8", -30, -8),
    ("-7 to -1", -7, -1),
    ("index day", 0, 0),
    ("+1 to +7", 1, 7),
    ("+8 to +45", 8, 45),
    ("+46 to +90", 46, 90),
    ("+91 to +365", 91, 365),
)


def presentation_dashboard_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "rhabdo_presentation_dashboard.html")


def presentation_asset_dir(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "presentation_dashboard")


def presentation_consort_svg_path(paths: ProjectPaths) -> str:
    return join_path(presentation_asset_dir(paths), "consort_flow.svg")


def presentation_timing_table_path(paths: ProjectPaths) -> str:
    return join_path(presentation_asset_dir(paths), "sepsis_renal_timing_oneyear.tsv")


def presentation_timing_svg_path(paths: ProjectPaths) -> str:
    return join_path(presentation_asset_dir(paths), "sepsis_renal_timing_oneyear.svg")


def presentation_train_prs_summary_path(paths: ProjectPaths) -> str:
    return join_path(presentation_asset_dir(paths), "train_prs_summary.tsv")


def presentation_train_prs_svg_path(paths: ProjectPaths) -> str:
    return join_path(presentation_asset_dir(paths), "train_prs_distribution.svg")


def _path_exists(path: str | None) -> bool:
    return bool(path) and not str(path).startswith("gs://") and Path(str(path)).exists()


def _relative_path(target: str, dashboard_path: str) -> str:
    try:
        return os.path.relpath(target, start=Path(dashboard_path).parent)
    except ValueError:
        return target


def _format_value(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    if isinstance(value, (int, np.integer)):
        return f"{int(value):,}"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return ""
    if number.is_integer():
        return f"{int(number):,}"
    if abs(number) < 0.001 and number != 0:
        return f"{number:.2e}"
    if abs(number) < 1:
        return f"{number:.3f}"
    return f"{number:,.2f}".rstrip("0").rstrip(".")


def _read_optional_table(path: str) -> pd.DataFrame:
    target = Path(path)
    if not target.exists():
        return pd.DataFrame()
    return read_table(str(target))


def _read_optional_json(path: str) -> dict[str, Any]:
    target = Path(path)
    if not target.exists():
        return {}
    try:
        return json.loads(target.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def _html_table(df: pd.DataFrame, *, columns: list[str] | None = None, limit: int | None = 25) -> str:
    if df.empty:
        return '<div class="placeholder">No rows available.</div>'
    view = df.copy()
    if columns:
        selected = [column for column in columns if column in view.columns]
        if selected:
            view = view[selected]
    if limit is not None:
        view = view.head(limit)
    for column in view.columns:
        if pd.api.types.is_numeric_dtype(view[column]):
            view[column] = view[column].map(_format_value)
    return view.to_html(index=False, escape=True, classes="data-table", border=0)


def _missing_panel(title: str, path: str, command: str = "") -> str:
    command_html = f"<p><code>{html.escape(command)}</code></p>" if command else ""
    return (
        '<div class="placeholder">'
        f"<strong>{html.escape(title)}</strong>"
        f"<p>Expected file not found: <code>{html.escape(path)}</code></p>"
        f"{command_html}"
        "</div>"
    )


def _card(title: str, value: object, subtitle: str = "") -> str:
    subtitle_html = f"<span>{html.escape(subtitle)}</span>" if subtitle else ""
    return (
        '<div class="metric-card">'
        f"<span>{html.escape(title)}</span>"
        f"<strong>{html.escape(_format_value(value) or 'NA')}</strong>"
        f"{subtitle_html}"
        "</div>"
    )


def _section(section_id: str, title: str, body: str) -> str:
    return f'<section id="{html.escape(section_id)}" class="section"><h2>{html.escape(title)}</h2>{body}</section>'


def _source(path: str, dashboard_path: str, label: str = "source") -> str:
    if not _path_exists(path):
        return f'<span class="source missing">{html.escape(label)} missing</span>'
    rel = _relative_path(path, dashboard_path)
    return f'<a class="source" href="{html.escape(rel)}">{html.escape(label)}</a>'


def _figure_card(title: str, figure_path: str, dashboard_path: str, caption: str = "") -> str:
    caption_html = f'<p class="note">{html.escape(caption)}</p>' if caption else ""
    if not _path_exists(figure_path):
        return f'<article class="panel"><div class="panel-head"><h3>{html.escape(title)}</h3></div>{caption_html}{_missing_panel(title, figure_path)}</article>'
    rel = _relative_path(figure_path, dashboard_path)
    return (
        '<article class="panel figure-panel">'
        f'<div class="panel-head"><h3>{html.escape(title)}</h3>{_source(figure_path, dashboard_path, "figure")}</div>'
        f"{caption_html}"
        f'<img src="{html.escape(rel)}" alt="{html.escape(title)}" />'
        "</article>"
    )


def _table_card(
    title: str,
    df: pd.DataFrame,
    *,
    dashboard_path: str,
    source_path: str | None = None,
    columns: list[str] | None = None,
    limit: int | None = 25,
    note: str = "",
) -> str:
    source_html = _source(source_path, dashboard_path, "table") if source_path else ""
    note_html = f'<p class="note">{html.escape(note)}</p>' if note else ""
    table_body = (
        _missing_panel(title, source_path)
        if source_path and not _path_exists(source_path) and df.empty
        else _html_table(df, columns=columns, limit=limit)
    )
    return (
        '<article class="panel">'
        f'<div class="panel-head"><h3>{html.escape(title)}</h3>{source_html}</div>'
        f"{note_html}"
        '<div class="table-wrap">'
        f"{table_body}"
        "</div></article>"
    )


def _consort_n(consort: pd.DataFrame, step: str) -> int | None:
    if consort.empty or not {"step", "n"}.issubset(consort.columns):
        return None
    match = consort[consort["step"].astype(str) == step]
    if match.empty:
        return None
    value = pd.to_numeric(match["n"], errors="coerce").dropna()
    return int(value.iloc[0]) if not value.empty else None


def _split_counts(split_summary: pd.DataFrame, group: str, column: str) -> int | None:
    if split_summary.empty or not {"group", column}.issubset(split_summary.columns):
        return None
    match = split_summary[split_summary["group"].astype(str) == group]
    if match.empty:
        return None
    value = pd.to_numeric(match[column], errors="coerce").dropna()
    return int(value.iloc[0]) if not value.empty else None


def _is_ck_confirmed_primary(config: ProjectConfig) -> bool:
    return str(config.cohort.primary_case_tier).lower() == "definite"


def _dashboard_title(config: ProjectConfig) -> str:
    if _is_ck_confirmed_primary(config):
        return "CK-Confirmed Rhabdomyolysis Susceptibility Dashboard"
    return "Rhabdomyolysis Analysis Dashboard"


def _primary_case_step(config: ProjectConfig) -> str:
    return "Definite rhabdo cases" if _is_ck_confirmed_primary(config) else "Broad rhabdo cases"


def _primary_case_card_label(config: ProjectConfig) -> str:
    return "CK-confirmed cases" if _is_ck_confirmed_primary(config) else "Broad cases"


def _primary_case_flow_label(config: ProjectConfig) -> str:
    if _is_ck_confirmed_primary(config):
        return "CK-confirmed rhabdomyolysis cases"
    return "Broad rhabdomyolysis cases"


def _write_consort_svg(config: ProjectConfig, consort: pd.DataFrame, split_summary: pd.DataFrame, path: str) -> None:
    rows = [
        ("Full built cohort", _consort_n(consort, "Built cohort rows")),
        ("EHR denominator eligible", _consort_n(consort, ">=2 OMOP condition dates")),
        (_primary_case_flow_label(config), _consort_n(consort, _primary_case_step(config))),
        ("Eligible controls", _consort_n(consort, "Eligible controls")),
        ("Matched analytic cohort", _consort_n(consort, "Matched analytic rows")),
        ("Matched cases", _consort_n(consort, "Matched cases")),
        ("Matched controls", _consort_n(consort, "Matched controls")),
        ("Training rows", _split_counts(split_summary, "train", "n_rows")),
        ("Testing rows", _split_counts(split_summary, "test", "n_rows")),
    ]
    excluded = _consort_n(consort, "CK >5000 without rhabdo, excluded")
    width, height = 980, 680
    box_w, box_h = 270, 64
    left_x, right_x = 130, 580
    y0, y_gap = 34, 84
    fill, border, text, muted = "#ffffff", "#cbd5e1", "#111827", "#64748b"
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="40" y="28" font-size="20" font-family="Arial, sans-serif" font-weight="700" fill="#111827">CONSORT Flow</text>',
    ]

    def box(x: int, y: int, label: str, value: int | None, accent: str = "#0f766e") -> None:
        lines.append(f'<rect x="{x}" y="{y}" width="{box_w}" height="{box_h}" rx="8" fill="{fill}" stroke="{border}"/>')
        lines.append(f'<rect x="{x}" y="{y}" width="6" height="{box_h}" rx="3" fill="{accent}"/>')
        lines.append(f'<text x="{x + 18}" y="{y + 25}" font-size="14" font-family="Arial, sans-serif" fill="{muted}">{html.escape(label)}</text>')
        lines.append(f'<text x="{x + 18}" y="{y + 50}" font-size="22" font-family="Arial, sans-serif" font-weight="700" fill="{text}">{html.escape(_format_value(value) if value is not None else "NA")}</text>')

    def arrow(x1: int, y1: int, x2: int, y2: int) -> None:
        lines.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#94a3b8" stroke-width="2" marker-end="url(#arrow)"/>')

    lines.append('<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#94a3b8"/></marker></defs>')
    main_labels = rows[:5]
    for index, (label, value) in enumerate(main_labels):
        y = y0 + index * y_gap
        box(left_x, y, label, value)
        if index:
            arrow(left_x + box_w // 2, y - 20, left_x + box_w // 2, y)
    if excluded is not None:
        y = y0 + 2 * y_gap
        box(right_x, y, "CK-only excluded", excluded, "#b45309")
        arrow(left_x + box_w, y + box_h // 2, right_x, y + box_h // 2)
    y = y0 + 5 * y_gap
    box(left_x, y, rows[5][0], rows[5][1], "#2563eb")
    box(right_x, y, rows[6][0], rows[6][1], "#2563eb")
    arrow(left_x + box_w // 2, y - 20, left_x + box_w // 2, y)
    arrow(left_x + box_w, y + box_h // 2, right_x, y + box_h // 2)
    y = y0 + 6 * y_gap
    box(left_x, y, rows[7][0], rows[7][1], "#7c3aed")
    box(right_x, y, rows[8][0], rows[8][1], "#7c3aed")
    arrow(left_x + box_w // 2, y - 20, left_x + box_w // 2, y)
    arrow(right_x + box_w // 2, y - 20, right_x + box_w // 2, y)
    lines.append("</svg>")
    write_text("".join(lines), path)


def _load_model_input(paths: ProjectPaths) -> pd.DataFrame:
    path = clinical_model_input_path(paths)
    return _read_optional_table(path)


def _cofactor_timing_summary(config: ProjectConfig, paths: ProjectPaths, model_input: pd.DataFrame) -> pd.DataFrame:
    if model_input.empty:
        return pd.DataFrame()
    events = load_clinical_cofactor_events(config)
    if events.empty:
        return pd.DataFrame()
    events = events.copy()
    events["person_id"] = events["person_id"].astype(str)
    events["cofactor"] = events["cofactor"].astype(str)
    events["condition_date"] = pd.to_datetime(events["condition_date"], errors="coerce")
    events = events[events["cofactor"].isin(["sepsis", "renal_injury"])].dropna(subset=["condition_date"]).copy()
    indexed = model_input[["person_id", "analysis_case", "index_date"]].copy()
    indexed["person_id"] = indexed["person_id"].astype(str)
    indexed["analysis_case"] = pd.to_numeric(indexed["analysis_case"], errors="coerce")
    indexed["case_control"] = np.where(indexed["analysis_case"].eq(1), "cases", "controls")
    indexed["index_date"] = pd.to_datetime(indexed["index_date"], errors="coerce")
    indexed = indexed.dropna(subset=["index_date", "analysis_case"]).copy()
    denominators = indexed.groupby("case_control")["person_id"].nunique().to_dict()
    merged = events.merge(indexed[["person_id", "case_control", "index_date"]], on="person_id", how="inner")
    if merged.empty:
        return pd.DataFrame()
    merged["delta_days"] = (merged["condition_date"] - merged["index_date"]).dt.days
    merged = merged[merged["delta_days"].between(-365, 365, inclusive="both")].copy()
    rows = []
    for case_control in ("cases", "controls"):
        denom = int(denominators.get(case_control, 0))
        for cofactor in ("sepsis", "renal_injury"):
            subset = merged[(merged["case_control"] == case_control) & (merged["cofactor"] == cofactor)]
            for label, low, high in TIMING_BINS:
                people = subset[subset["delta_days"].between(low, high, inclusive="both")]["person_id"].nunique()
                rows.append(
                    {
                        "case_control": case_control,
                        "cofactor": cofactor,
                        "timing_bin": label,
                        "bin_start_day": low,
                        "bin_end_day": high,
                        "n_participants": int(people),
                        "denominator": denom,
                        "pct_participants": float(100 * people / denom) if denom else 0.0,
                    }
                )
    return pd.DataFrame(rows)


def _write_timing_svg(summary: pd.DataFrame, path: str) -> None:
    width, height = 1180, 620
    left, top, bottom, right = 70, 70, 82, 30
    panel_w = (width - left - right - 34) / 2
    panel_h = (height - top - bottom - 34) / 2
    colors = {"sepsis": "#0f766e", "renal_injury": "#b45309"}
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="40" y="30" font-size="20" font-family="Arial, sans-serif" font-weight="700" fill="#111827">Sepsis And Renal Injury Timing</text>',
        '<text x="40" y="52" font-size="12" font-family="Arial, sans-serif" fill="#64748b">Participant-level event counts within one year before and after index.</text>',
    ]
    if summary.empty:
        lines.append('<text x="40" y="120" font-size="16" font-family="Arial, sans-serif" fill="#64748b">No timing rows available.</text></svg>')
        write_text("".join(lines), path)
        return
    max_pct = max(float(pd.to_numeric(summary["pct_participants"], errors="coerce").max()), 1.0)
    max_pct = max(max_pct * 1.15, 5.0)
    bins = [label for label, _, _ in TIMING_BINS]
    panels = [("cases", "sepsis"), ("cases", "renal_injury"), ("controls", "sepsis"), ("controls", "renal_injury")]
    for index, (case_control, cofactor) in enumerate(panels):
        col, row = index % 2, index // 2
        x0 = left + col * (panel_w + 34)
        y0 = top + row * (panel_h + 34)
        title = f"{case_control.title()} - {cofactor.replace('_', ' ').title()}"
        lines.append(f'<text x="{x0}" y="{y0 - 12}" font-size="14" font-family="Arial, sans-serif" font-weight="700" fill="#111827">{html.escape(title)}</text>')
        lines.append(f'<line x1="{x0}" y1="{y0 + panel_h}" x2="{x0 + panel_w}" y2="{y0 + panel_h}" stroke="#334155"/>')
        lines.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0 + panel_h}" stroke="#334155"/>')
        chunk = summary[(summary["case_control"] == case_control) & (summary["cofactor"] == cofactor)].copy()
        bar_gap = 5
        bar_w = max((panel_w - bar_gap * (len(bins) - 1)) / len(bins), 1)
        for bin_index, bin_label in enumerate(bins):
            row_match = chunk[chunk["timing_bin"] == bin_label]
            pct = float(row_match["pct_participants"].iloc[0]) if not row_match.empty else 0.0
            count = int(row_match["n_participants"].iloc[0]) if not row_match.empty else 0
            h = (pct / max_pct) * (panel_h - 18)
            x = x0 + bin_index * (bar_w + bar_gap)
            y = y0 + panel_h - h
            lines.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{colors[cofactor]}" opacity="0.86"/>')
            if pct >= max_pct * 0.12:
                lines.append(f'<text x="{x + bar_w / 2:.1f}" y="{y - 3:.1f}" text-anchor="middle" font-size="9" font-family="Arial, sans-serif" fill="#334155">{count}</text>')
            lines.append(f'<text x="{x + bar_w / 2:.1f}" y="{y0 + panel_h + 14}" text-anchor="end" font-size="8" transform="rotate(-35 {x + bar_w / 2:.1f} {y0 + panel_h + 14})" font-family="Arial, sans-serif" fill="#64748b">{html.escape(bin_label)}</text>')
        lines.append(f'<text x="{x0 + 5}" y="{y0 + 12}" font-size="10" font-family="Arial, sans-serif" fill="#64748b">max {max_pct:.1f}%</text>')
    lines.append("</svg>")
    write_text("".join(lines), path)


def _score_column(frame: pd.DataFrame) -> str | None:
    for column in ("SCORE1_SUM", "SCORE1_AVG", "SCORE_SUM", "SCORE_AVG"):
        if column in frame.columns:
            return column
    score_columns = [column for column in frame.columns if "SCORE" in column.upper() or "PRS" in column.upper()]
    return score_columns[-1] if score_columns else None


def _load_train_prs(paths: ProjectPaths, clinical_prs_label: str, model_input: pd.DataFrame) -> pd.DataFrame:
    score_path = Path(clinical_prs_model_dir(paths, clinical_prs_label)) / "train_prs.sscore"
    if not score_path.exists() or model_input.empty:
        return pd.DataFrame()
    raw = pd.read_csv(score_path, sep=r"\s+", dtype=str)
    iid = "IID" if "IID" in raw.columns else "#IID" if "#IID" in raw.columns else None
    score_col = _score_column(raw)
    if iid is None or score_col is None:
        return pd.DataFrame()
    scores = pd.DataFrame(
        {
            "person_id": raw[iid].astype(str),
            "prs_score": pd.to_numeric(raw[score_col], errors="coerce"),
        }
    ).dropna(subset=["prs_score"])
    train = model_input[model_input["analysis_split"].astype(str) == "train"].copy()
    train["person_id"] = train["person_id"].astype(str)
    train["analysis_case"] = pd.to_numeric(train["analysis_case"], errors="coerce")
    return scores.merge(train[["person_id", "analysis_case"]], on="person_id", how="inner").dropna(subset=["analysis_case"])


def _prs_summary(train_prs: pd.DataFrame) -> pd.DataFrame:
    if train_prs.empty:
        return pd.DataFrame()
    rows = []
    for label, group in (("cases", train_prs[train_prs["analysis_case"] == 1]), ("controls", train_prs[train_prs["analysis_case"] == 0])):
        rows.append(
            {
                "group": label,
                "n": int(len(group)),
                "mean": float(group["prs_score"].mean()) if not group.empty else np.nan,
                "sd": float(group["prs_score"].std(ddof=0)) if not group.empty else np.nan,
                "median": float(group["prs_score"].median()) if not group.empty else np.nan,
                "q1": float(group["prs_score"].quantile(0.25)) if not group.empty else np.nan,
                "q3": float(group["prs_score"].quantile(0.75)) if not group.empty else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _write_train_prs_svg(train_prs: pd.DataFrame, path: str) -> None:
    width, height = 920, 360
    left, right, top, bottom = 65, 28, 54, 64
    inner_w, inner_h = width - left - right, height - top - bottom
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="36" y="28" font-size="20" font-family="Arial, sans-serif" font-weight="700" fill="#111827">Training PRS Distribution</text>',
    ]
    if train_prs.empty:
        lines.append('<text x="36" y="88" font-size="14" font-family="Arial, sans-serif" fill="#64748b">No train PRS scores available.</text></svg>')
        write_text("".join(lines), path)
        return
    values = pd.to_numeric(train_prs["prs_score"], errors="coerce").dropna()
    q_low, q_high = float(values.quantile(0.01)), float(values.quantile(0.99))
    if q_low == q_high:
        q_low, q_high = float(values.min()), float(values.max())
    if q_low == q_high:
        q_low -= 1
        q_high += 1
    bins = np.linspace(q_low, q_high, 31)
    max_density = 0.0
    histograms: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for label, group in (("cases", train_prs[train_prs["analysis_case"] == 1]), ("controls", train_prs[train_prs["analysis_case"] == 0])):
        clipped = pd.to_numeric(group["prs_score"], errors="coerce").clip(q_low, q_high).dropna()
        counts, edges = np.histogram(clipped, bins=bins)
        density = counts / max(counts.sum(), 1)
        max_density = max(max_density, float(density.max()) if len(density) else 0.0)
        histograms[label] = (density, edges)
    max_density = max(max_density, 0.01)
    colors = {"cases": "#991b1b", "controls": "#0f766e"}
    lines.append(f'<line x1="{left}" y1="{top + inner_h}" x2="{left + inner_w}" y2="{top + inner_h}" stroke="#334155"/>')
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + inner_h}" stroke="#334155"/>')
    for label, (density, edges) in histograms.items():
        points = []
        for index, value in enumerate(density):
            mid = (edges[index] + edges[index + 1]) / 2
            x = left + ((mid - q_low) / (q_high - q_low)) * inner_w
            y = top + inner_h - (float(value) / max_density) * inner_h
            points.append(f"{x:.1f},{y:.1f}")
        lines.append(f'<polyline points="{" ".join(points)}" fill="none" stroke="{colors[label]}" stroke-width="3" opacity="0.9"/>')
    lines.append(f'<text x="{left}" y="{height - 18}" font-size="11" font-family="Arial, sans-serif" fill="#64748b">PRS score; x-axis clipped to 1st-99th percentile for display</text>')
    lines.append(f'<circle cx="{width - 190}" cy="25" r="5" fill="{colors["cases"]}"/><text x="{width - 178}" y="29" font-size="12" font-family="Arial, sans-serif">cases</text>')
    lines.append(f'<circle cx="{width - 115}" cy="25" r="5" fill="{colors["controls"]}"/><text x="{width - 103}" y="29" font-size="12" font-family="Arial, sans-serif">controls</text>')
    lines.append("</svg>")
    write_text("".join(lines), path)


def _model_comparison_table(paths: ProjectPaths, gwas_label: str, prs_label: str, clinical_prs_label: str) -> pd.DataFrame:
    rows = []
    clinical = _read_optional_table(clinical_model_metrics_path(paths))
    if not clinical.empty:
        test = clinical[clinical["evaluation_set"].astype(str) == "test"].head(1)
        if not test.empty:
            row = test.iloc[0].to_dict()
            rows.append({"model": "Clinical only", **row})
    prs = _read_optional_table(microarray_prs_metrics_path(paths, gwas_label, prs_label))
    if not prs.empty:
        row = prs.head(1).iloc[0].to_dict()
        rows.append(
            {
                "model": "PRS only",
                "evaluation_set": "test",
                "n": row.get("n_participants"),
                "cases": row.get("n_cases"),
                "controls": row.get("n_controls"),
                "roc_auc": row.get("roc_auc"),
                "average_precision": row.get("average_precision"),
                "brier_score": "",
                "log_loss": "",
            }
        )
    clinical_prs = _read_optional_table(clinical_prs_model_metrics_path(paths, clinical_prs_label))
    if not clinical_prs.empty:
        test = clinical_prs[clinical_prs["evaluation_set"].astype(str) == "test"].head(1)
        if not test.empty:
            row = test.iloc[0].to_dict()
            rows.append({"model": "Clinical + PRS", **row})
    if not rows:
        return pd.DataFrame()
    frame = pd.DataFrame(rows)
    columns = ["model", "n", "cases", "controls", "roc_auc", "average_precision", "brier_score", "log_loss", "sensitivity", "specificity", "ppv", "npv"]
    return frame[[column for column in columns if column in frame.columns]]


def _driver_table(paths: ProjectPaths, clinical_prs_label: str) -> pd.DataFrame:
    coefficients = _read_optional_table(clinical_prs_model_coefficients_path(paths, clinical_prs_label))
    if coefficients.empty:
        return coefficients
    coefficients = coefficients.copy()
    coefficients["abs_beta"] = pd.to_numeric(coefficients.get("beta"), errors="coerce").abs()
    coefficients = coefficients[coefficients["feature"].astype(str) != "intercept"].sort_values("abs_beta", ascending=False)
    columns = ["feature", "source_column", "kind", "reference", "beta", "odds_ratio"]
    return coefficients[[column for column in columns if column in coefficients.columns]].head(20)


def _gwas_cards(paths: ProjectPaths, gwas_label: str, prs_label: str) -> str:
    qc = _read_optional_json(microarray_plink_qc_path(paths, gwas_label))
    variant_qc = _read_optional_table(microarray_plink_variant_qc_summary_path(paths, gwas_label))
    weights = _read_optional_table(microarray_prs_weights_path(paths, gwas_label, prs_label))
    sample = qc.get("sample_counts", {}) if isinstance(qc.get("sample_counts"), dict) else {}
    return (
        '<div class="metric-grid">'
        + _card("Data source", "AoU v8 microarray", "PLINK BED/BIM/FAM")
        + _card("Chromosomes", ", ".join(map(str, qc.get("chromosomes_tested", []))) if qc else "NA")
        + _card("Train cases", sample.get("after_microarray_fam_overlap_cases"))
        + _card("Train controls", sample.get("after_microarray_fam_overlap_controls"))
        + _card("Variants tested", qc.get("n_variants_tested") if qc else None)
        + _card("PRS clumped variants", len(weights) if not weights.empty else None)
        + _card("MAF filter", qc.get("min_maf") if qc else None)
        + _card("Call rate filter", qc.get("min_call_rate") if qc else None)
        + "</div>"
        + (
            _html_table(variant_qc, columns=["filter", "rows_before", "rows_after", "rows_removed"], limit=None)
            if not variant_qc.empty
            else ""
        )
    )


def _definition_cards(config: ProjectConfig) -> str:
    definitions = [
        ("Study denominator", "Observation-period data and at least two distinct OMOP condition record dates."),
        ("Broad case", "At least one rhabdomyolysis OMOP condition record; CK not required."),
        ("Definite case", "Broad case plus CK >5000 from 7 days before through 45 days after first rhabdo diagnosis."),
        ("Control", "No rhabdomyolysis diagnosis and no CK >5000 evidence; missing CK allowed."),
        ("Excluded", "CK >5000 without rhabdomyolysis diagnosis is indeterminate and excluded from primary case-control analysis."),
    ]
    if _is_ck_confirmed_primary(config):
        definitions.append(
            (
                "Primary model/GWAS",
                "CK-confirmed cases versus matched eligible controls; crush/major trauma is excluded from primary eligibility, while sepsis and renal injury are characterized but not excluded.",
            )
        )
    else:
        definitions.append(
            (
                "Primary model/GWAS",
                "Broad cases versus matched eligible controls, using primary non-traumatic eligibility.",
            )
        )
    return '<div class="definition-grid">' + "".join(
        f'<article class="definition"><h3>{html.escape(title)}</h3><p>{html.escape(text)}</p></article>'
        for title, text in definitions
    ) + "</div>"


def _split_table(split_table: pd.DataFrame, split: str) -> pd.DataFrame:
    if split_table.empty or "analysis_split" not in split_table.columns:
        return pd.DataFrame()
    return split_table[split_table["analysis_split"].astype(str) == split].copy()


def _css() -> str:
    return """
:root {
  --bg: #f4f6f8;
  --panel: #ffffff;
  --ink: #172033;
  --muted: #667085;
  --border: #d7dde6;
  --accent: #0f766e;
  --accent2: #1d4ed8;
  --warn: #b45309;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--bg);
  color: var(--ink);
  font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  line-height: 1.45;
}
.shell { max-width: 1240px; margin: 0 auto; padding: 30px 24px 60px; }
.hero {
  padding: 28px 0 24px;
  border-bottom: 1px solid var(--border);
  margin-bottom: 22px;
}
.hero h1 { margin: 0 0 8px; font-size: 32px; letter-spacing: 0; }
.hero p { margin: 0; color: var(--muted); }
.section { margin: 26px 0; }
.section h2 { font-size: 21px; margin: 0 0 14px; letter-spacing: 0; }
.grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 16px; }
.metric-grid { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; margin-bottom: 16px; }
.definition-grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 12px; }
.metric-card, .definition, .panel {
  background: var(--panel);
  border: 1px solid var(--border);
  border-radius: 8px;
  box-shadow: 0 1px 2px rgba(16, 24, 40, 0.04);
}
.metric-card { padding: 14px 16px; min-height: 92px; }
.metric-card span { display: block; color: var(--muted); font-size: 12px; }
.metric-card strong { display: block; margin-top: 6px; font-size: 23px; letter-spacing: 0; }
.definition { padding: 16px; }
.definition h3, .panel h3 { margin: 0; font-size: 15px; }
.definition p, .note { color: var(--muted); font-size: 13px; margin: 8px 0 0; }
.panel { padding: 16px; overflow: hidden; }
.panel.full { grid-column: 1 / -1; }
.panel-head { display: flex; align-items: center; justify-content: space-between; gap: 12px; margin-bottom: 12px; }
.source { color: var(--accent2); font-size: 12px; text-decoration: none; white-space: nowrap; }
.source.missing { color: var(--warn); }
.figure-panel img { width: 100%; height: auto; display: block; }
.table-wrap { overflow-x: auto; }
.data-table { width: 100%; border-collapse: collapse; font-size: 12px; }
.data-table th {
  text-align: left;
  background: #eef2f7;
  color: #344054;
  padding: 8px;
  border-bottom: 1px solid var(--border);
  white-space: nowrap;
}
.data-table td {
  padding: 7px 8px;
  border-bottom: 1px solid #edf0f4;
  vertical-align: top;
  white-space: nowrap;
}
.placeholder {
  border: 1px dashed #cbd5e1;
  border-radius: 8px;
  padding: 18px;
  color: var(--muted);
  background: #f8fafc;
  font-size: 13px;
}
code { background: #eef2f7; padding: 1px 4px; border-radius: 4px; }
@media (max-width: 980px) {
  .grid, .definition-grid { grid-template-columns: 1fr; }
  .metric-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }
  .shell { padding: 20px 14px 42px; }
}
"""


def render_presentation_dashboard(
    config: ProjectConfig,
    paths: ProjectPaths,
    *,
    gwas_label: str = "microarray_plink_autosomes_maf05_train_qc",
    prs_label: str = "test-clumped-p001",
    clinical_prs_label: str = "clinical_prs_p001",
    diagnostics_label: str = "prs_p001_diagnostics",
    output_path: str | None = None,
) -> str:
    output = output_path or presentation_dashboard_path(paths)
    ensure_parent_dir(output)
    ensure_parent_dir(join_path(presentation_asset_dir(paths), "placeholder.txt"))

    consort = _read_optional_table(consort_counts_path(paths))
    split_summary = _read_optional_table(model_split_summary_path(paths))
    split_table = _read_optional_table(split_table1_path(paths))
    model_input = _load_model_input(paths)

    _write_consort_svg(config, consort, split_summary, presentation_consort_svg_path(paths))
    timing = _cofactor_timing_summary(config, paths, model_input)
    write_dataframe(timing, presentation_timing_table_path(paths))
    _write_timing_svg(timing, presentation_timing_svg_path(paths))

    train_prs = _load_train_prs(paths, clinical_prs_label, model_input)
    train_summary = _prs_summary(train_prs)
    write_dataframe(train_summary, presentation_train_prs_summary_path(paths))
    _write_train_prs_svg(train_prs, presentation_train_prs_svg_path(paths))

    model_comparison = _model_comparison_table(paths, gwas_label, prs_label, clinical_prs_label)
    drivers = _driver_table(paths, clinical_prs_label)
    manhattan = microarray_plink_manhattan_path(paths, gwas_label)
    gwas_results = microarray_plink_results_path(paths, gwas_label)
    diagnostics_path = join_path(paths.run_root, "clinical", "prs_diagnostics", slugify(diagnostics_label), "overall_metrics.tsv")

    cards = (
        '<div class="metric-grid">'
        + _card("Built cohort", _consort_n(consort, "Built cohort rows"))
        + _card(_primary_case_card_label(config), _consort_n(consort, _primary_case_step(config)))
        + _card("Eligible controls", _consort_n(consort, "Eligible controls"))
        + _card("Matched rows", _consort_n(consort, "Matched analytic rows"))
        + "</div>"
    )

    sections = [
        _section(
            "cohort-flow",
            "Participant Flow",
            cards
            + '<div class="grid">'
            + _figure_card("CONSORT Flow", presentation_consort_svg_path(paths), output)
            + '<article class="panel"><div class="panel-head"><h3>Case And Control Definitions</h3></div>'
            + _definition_cards(config)
            + "</article></div>",
        ),
        _section(
            "table1",
            "Clinical Characteristics",
            '<div class="grid">'
            + _table_card(
                "Training Set Table 1",
                _split_table(split_table, "train"),
                dashboard_path=output,
                source_path=split_table1_path(paths),
                limit=None,
                columns=["section", "variable", "matched_cases", "matched_controls", "effect_measure", "effect_size", "ci_95", "smd", "p_value"],
            )
            + _table_card(
                "Testing Set Table 1",
                _split_table(split_table, "test"),
                dashboard_path=output,
                source_path=split_table1_path(paths),
                limit=None,
                columns=["section", "variable", "matched_cases", "matched_controls", "effect_measure", "effect_size", "ci_95", "smd", "p_value"],
            )
            + "</div>",
        ),
        _section(
            "timing",
            "Sepsis And Renal Injury Timing",
            '<div class="grid">'
            + _figure_card(
                "Participant-Level Timing Within One Year",
                presentation_timing_svg_path(paths),
                output,
                "Each bin counts participants with at least one event in the bin.",
            )
            + _table_card(
                "Timing Bin Counts",
                timing,
                dashboard_path=output,
                source_path=presentation_timing_table_path(paths),
                limit=None,
            )
            + "</div>",
        ),
        _section(
            "gwas",
            "GWAS And PRS Inputs",
            '<div class="grid">'
            + '<article class="panel full"><div class="panel-head"><h3>GWAS Parameters</h3>'
            + _source(gwas_results, output, "GWAS results")
            + "</div>"
            + _gwas_cards(paths, gwas_label, prs_label)
            + "</article>"
            + _figure_card("Manhattan Plot", manhattan, output)
            + _figure_card(
                "Training PRS Distribution",
                presentation_train_prs_svg_path(paths),
                output,
                "Train-split PRS scores from the clumped p<=0.01 weights.",
            )
            + _table_card(
                "Training PRS Summary",
                train_summary,
                dashboard_path=output,
                source_path=presentation_train_prs_summary_path(paths),
                limit=None,
            )
            + "</div>",
        ),
        _section(
            "models",
            "Model Performance And Drivers",
            '<div class="grid">'
            + _table_card(
                "Model Comparison",
                model_comparison,
                dashboard_path=output,
                source_path=diagnostics_path,
                limit=None,
            )
            + _table_card(
                "Clinical + PRS Model Drivers",
                drivers,
                dashboard_path=output,
                source_path=clinical_prs_model_coefficients_path(paths, clinical_prs_label),
                limit=20,
            )
            + "</div>",
        ),
    ]
    nav = "".join(
        f'<a href="#{section_id}">{html.escape(label)}</a>'
        for section_id, label in [
            ("cohort-flow", "Flow"),
            ("table1", "Table 1"),
            ("timing", "Timing"),
            ("gwas", "GWAS/PRS"),
            ("models", "Models"),
        ]
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(config.analysis.analysis_name)} Presentation Dashboard</title>
  <style>{_css()}</style>
</head>
<body>
  <div class="shell">
    <header class="hero">
      <h1>{html.escape(_dashboard_title(config))}</h1>
      <p>Presentation summary for supervisor review. Primary case tier: <code>{html.escape(str(config.cohort.primary_case_tier))}</code>. Run root: <code>{html.escape(paths.run_root)}</code></p>
      <nav class="source" style="margin-top: 12px; display: flex; gap: 14px; flex-wrap: wrap;">{nav}</nav>
    </header>
    {''.join(sections)}
  </div>
</body>
</html>
"""
    write_text(document, output)
    return output


__all__ = [
    "presentation_consort_svg_path",
    "presentation_dashboard_path",
    "presentation_timing_svg_path",
    "presentation_timing_table_path",
    "presentation_train_prs_svg_path",
    "presentation_train_prs_summary_path",
    "render_presentation_dashboard",
]
