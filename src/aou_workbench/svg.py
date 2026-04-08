"""Minimal SVG writers for GWAS visualization without heavyweight plotting dependencies."""

from __future__ import annotations

import math

import pandas as pd

from .io_utils import write_text


def _placeholder_svg(title: str, message: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="900" height="240">'
        f'<rect width="100%" height="100%" fill="#fffaf0" />'
        f'<text x="40" y="80" font-size="28" font-family="Menlo, monospace">{title}</text>'
        f'<text x="40" y="130" font-size="18" font-family="Menlo, monospace">{message}</text>'
        "</svg>"
    )


def _chromosome_rank(value: str) -> int:
    text = str(value).replace("chr", "").upper()
    if text == "X":
        return 23
    if text == "Y":
        return 24
    try:
        return int(text)
    except ValueError:
        return 99


def write_manhattan_svg(df: pd.DataFrame, path: str, *, pvalue_column: str = "regression_p") -> None:
    if df.empty or pvalue_column not in df.columns:
        write_text(_placeholder_svg("Manhattan Plot", "No GWAS results available."), path)
        return
    plot_df = df.copy()
    plot_df = plot_df[pd.to_numeric(plot_df[pvalue_column], errors="coerce").notna()].copy()
    if plot_df.empty:
        write_text(_placeholder_svg("Manhattan Plot", "No valid p-values available."), path)
        return
    plot_df["chromosome_rank"] = plot_df["chromosome"].map(_chromosome_rank)
    plot_df["position"] = pd.to_numeric(plot_df["position"], errors="coerce")
    plot_df["neg_log10_p"] = -plot_df[pvalue_column].astype(float).clip(lower=1e-300).map(math.log10)
    plot_df = plot_df.sort_values(["chromosome_rank", "position"]).reset_index(drop=True)
    max_y = max(float(plot_df["neg_log10_p"].max()), 8.0)
    width = 1100
    height = 460
    left = 70
    right = 30
    top = 30
    bottom = 60
    inner_width = width - left - right
    inner_height = height - top - bottom

    cumulative = 0.0
    x_values = []
    labels = []
    for chrom, chunk in plot_df.groupby("chromosome", sort=False):
        max_pos = float(chunk["position"].max())
        for pos in chunk["position"]:
            x_values.append(cumulative + float(pos))
        labels.append((chrom, cumulative + max_pos / 2))
        cumulative += max_pos + 1_000_000.0
    plot_df["plot_x_raw"] = x_values
    total_width = max(float(plot_df["plot_x_raw"].max()), 1.0)
    plot_df["plot_x"] = left + (plot_df["plot_x_raw"] / total_width) * inner_width
    plot_df["plot_y"] = top + inner_height - (plot_df["neg_log10_p"] / max_y) * inner_height

    points = []
    for index, row in plot_df.iterrows():
        color = "#0f766e" if index % 2 == 0 else "#7c2d12"
        points.append(
            f'<circle cx="{row.plot_x:.2f}" cy="{row.plot_y:.2f}" r="2.2" fill="{color}" opacity="0.8" />'
        )

    threshold = -math.log10(5e-8)
    threshold_y = top + inner_height - (threshold / max_y) * inner_height
    x_axis = f'<line x1="{left}" y1="{top + inner_height}" x2="{width - right}" y2="{top + inner_height}" stroke="#1f2937" />'
    y_axis = f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + inner_height}" stroke="#1f2937" />'
    labels_svg = "".join(
        f'<text x="{left + (value / total_width) * inner_width:.1f}" y="{height - 18}" text-anchor="middle" font-size="12" font-family="Menlo, monospace">{chrom}</text>'
        for chrom, value in labels
    )
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
        '<rect width="100%" height="100%" fill="#fffdf7" />'
        f'<text x="{left}" y="20" font-size="24" font-family="Menlo, monospace">GWAS Manhattan Plot</text>'
        f'{x_axis}{y_axis}'
        f'<line x1="{left}" y1="{threshold_y:.2f}" x2="{width - right}" y2="{threshold_y:.2f}" stroke="#dc2626" stroke-dasharray="6 4" />'
        f'{"".join(points)}'
        f"{labels_svg}"
        f'<text x="18" y="{top + 20}" font-size="12" transform="rotate(-90 18 {top + 20})" font-family="Menlo, monospace">-log10(p)</text>'
        "</svg>"
    )
    write_text(svg, path)


def write_qq_svg(df: pd.DataFrame, path: str, *, pvalue_column: str = "regression_p") -> None:
    if df.empty or pvalue_column not in df.columns:
        write_text(_placeholder_svg("QQ Plot", "No GWAS results available."), path)
        return
    plot_df = df[pd.to_numeric(df[pvalue_column], errors="coerce").notna()].copy()
    if plot_df.empty:
        write_text(_placeholder_svg("QQ Plot", "No valid p-values available."), path)
        return
    plot_df = plot_df.sort_values(pvalue_column).reset_index(drop=True)
    plot_df["observed"] = -plot_df[pvalue_column].astype(float).clip(lower=1e-300).map(math.log10)
    plot_df["expected"] = [
        -math.log10((index + 1) / (len(plot_df) + 1)) for index in range(len(plot_df))
    ]
    max_axis = max(float(plot_df["observed"].max()), float(plot_df["expected"].max()), 6.0)
    width = 520
    height = 520
    left = 70
    right = 30
    top = 30
    bottom = 70
    inner_width = width - left - right
    inner_height = height - top - bottom
    plot_df["plot_x"] = left + (plot_df["expected"] / max_axis) * inner_width
    plot_df["plot_y"] = top + inner_height - (plot_df["observed"] / max_axis) * inner_height
    points = "".join(
        f'<circle cx="{row.plot_x:.2f}" cy="{row.plot_y:.2f}" r="2.2" fill="#0f766e" opacity="0.85" />'
        for row in plot_df.itertuples(index=False)
    )
    diagonal = (
        f'<line x1="{left}" y1="{top + inner_height}" x2="{width - right}" y2="{top}" '
        'stroke="#dc2626" stroke-dasharray="6 4" />'
    )
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
        '<rect width="100%" height="100%" fill="#fffdf7" />'
        f'<text x="{left}" y="20" font-size="24" font-family="Menlo, monospace">GWAS QQ Plot</text>'
        f'<line x1="{left}" y1="{top + inner_height}" x2="{width - right}" y2="{top + inner_height}" stroke="#1f2937" />'
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + inner_height}" stroke="#1f2937" />'
        f"{diagonal}{points}"
        f'<text x="{width / 2:.1f}" y="{height - 18}" text-anchor="middle" font-size="12" font-family="Menlo, monospace">Expected -log10(p)</text>'
        f'<text x="18" y="{height / 2:.1f}" font-size="12" transform="rotate(-90 18 {height / 2:.1f})" font-family="Menlo, monospace">Observed -log10(p)</text>'
        "</svg>"
    )
    write_text(svg, path)


__all__ = ["write_manhattan_svg", "write_qq_svg"]
