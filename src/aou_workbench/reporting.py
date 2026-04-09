"""Markdown reporting helpers for stage and end-to-end summaries."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

from .io_utils import write_text


def dataframe_markdown(df: pd.DataFrame, *, columns: Iterable[str] | None = None, limit: int = 10) -> str:
    if df.empty:
        return "_No rows available._"
    view = df.copy()
    if columns is not None:
        selected = [column for column in columns if column in view.columns]
        if selected:
            view = view[selected]
    return view.head(limit).to_markdown(index=False)


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


def write_final_report(
    *,
    analysis_name: str,
    output_root: str,
    stage1: pd.DataFrame,
    stage2_genes: pd.DataFrame,
    stage3: pd.DataFrame,
    stage4_hits: pd.DataFrame,
    path: str,
) -> None:
    sections = [
        f"# {analysis_name} Rhabdomyolysis Genetics Summary",
        "",
        f"Run root: `{output_root}`",
        "",
        "## Stage 1: A priori variants",
        "",
        dataframe_markdown(stage1, columns=["label", "gene", "case_carriers", "control_carriers", "fisher_p"], limit=8),
        "",
        "## Stage 2: P/LP panel summary",
        "",
        dataframe_markdown(stage2_genes, columns=["gene", "case_carriers", "control_carriers", "regression_p"], limit=8),
        "",
        "## Stage 3: Burden results",
        "",
        dataframe_markdown(stage3, columns=["gene", "mask", "case_carriers", "control_carriers", "regression_p"], limit=8),
        "",
        "## Stage 4: GWAS lead hits",
        "",
        dataframe_markdown(stage4_hits, columns=["variant_id", "gene", "chromosome", "position", "regression_p"], limit=10),
    ]
    write_text("\n".join(sections), path)


__all__ = ["load_table_if_exists", "write_final_report", "write_stage_report"]
