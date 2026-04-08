"""Stage 4: common-variant GWAS on the matched cohort."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd

from .config import ProjectConfig
from .io_utils import read_table, write_dataframe, write_json
from .paths import ProjectPaths
from .reporting import write_stage_report
from .statistics import bh_fdr, genomic_control_lambda, run_binary_logistic_regression
from .svg import write_manhattan_svg, write_qq_svg


def _lead_hit_subset(df: pd.DataFrame, window_bp: int) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    chosen: list[pd.Series] = []
    for row in df.sort_values("regression_p").itertuples(index=False):
        keep = True
        for hit in chosen:
            if str(hit.chromosome) == str(row.chromosome) and abs(int(hit.position) - int(row.position)) <= window_bp:
                keep = False
                break
        if keep:
            chosen.append(pd.Series(row._asdict()))
    return pd.DataFrame(chosen)


def run_stage4_gwas(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    stage = config.analysis.stage4
    if stage is None:
        return pd.DataFrame(), pd.DataFrame()
    raw = read_table(stage.genotype_table).copy()
    raw["person_id"] = raw[stage.person_id_column].astype(str)
    raw["variant_id"] = raw[stage.variant_id_column].astype(str)
    raw["dosage"] = pd.to_numeric(raw[stage.dosage_column], errors="coerce").fillna(0.0)
    raw["chromosome"] = raw[stage.chromosome_column].astype(str)
    raw["position"] = pd.to_numeric(raw[stage.position_column], errors="coerce")
    raw["gene"] = raw.get(stage.gene_column, pd.Series([""] * len(raw), index=raw.index)).astype(str)
    raw = raw[raw["person_id"].isin(set(matched_df["person_id"]))].copy()

    annotation = pd.DataFrame()
    if stage.annotation_table:
        annotation = read_table(stage.annotation_table).copy()
        annotation["variant_id"] = annotation[stage.variant_id_column].astype(str)

    rows: list[dict[str, Any]] = []
    n_samples = matched_df["person_id"].nunique()
    for variant_id, chunk in raw.groupby("variant_id"):
        exposure = chunk.groupby("person_id")["dosage"].max()
        mac = float(exposure.sum())
        meta = chunk.iloc[0]
        if stage.af_column in chunk.columns:
            maf = pd.to_numeric(chunk[stage.af_column], errors="coerce").dropna()
            maf_value = float(maf.iloc[0]) if not maf.empty else mac / max(2 * n_samples, 1)
        else:
            maf_value = mac / max(2 * n_samples, 1)
        if maf_value < stage.min_maf or mac < stage.min_minor_allele_count:
            continue
        regression = run_binary_logistic_regression(
            matched_df,
            exposure,
            outcome_column=config.analysis.matched_outcome_column,
            covariates=stage.covariates,
        )
        rows.append(
            {
                "variant_id": variant_id,
                "chromosome": meta["chromosome"],
                "position": int(meta["position"]),
                "gene": meta["gene"],
                "maf": maf_value,
                "minor_allele_count": mac,
                **regression,
            }
        )
    full = pd.DataFrame(rows)
    if not full.empty:
        full = full.sort_values(["regression_p", "variant_id"]).reset_index(drop=True)
    if not full.empty:
        full["fdr_q"] = bh_fdr(full["regression_p"].values)
        full["minus_log10_p"] = -np.log10(full["regression_p"].clip(lower=1e-300))
        if not annotation.empty:
            full = full.merge(annotation.drop_duplicates("variant_id"), on="variant_id", how="left", suffixes=("", "_annotation"))
    lead_hits = _lead_hit_subset(full, stage.lead_hit_window_bp)
    write_dataframe(full, paths.stage4_full_results_tsv)
    write_dataframe(lead_hits, paths.stage4_lead_hits_tsv)
    write_manhattan_svg(full, paths.stage4_manhattan_svg)
    write_qq_svg(full, paths.stage4_qq_svg)
    write_json(
        {
            "n_variants_tested": int(full.shape[0]),
            "n_lead_hits": int(lead_hits.shape[0]),
            "lambda_gc": genomic_control_lambda(full["regression_p"].values if not full.empty else []),
        },
        paths.stage4_qc_json,
    )
    write_stage_report(
        title="Stage 4: GWAS",
        summary_lines=[
            f"- Variants tested: {full.shape[0]}",
            f"- Lead hits: {lead_hits.shape[0]}",
            f"- MAF threshold: {stage.min_maf}",
        ],
        preview_df=lead_hits if not lead_hits.empty else full,
        preview_columns=["variant_id", "gene", "chromosome", "position", "regression_p", "fdr_q"],
        path=paths.stage4_report_md,
    )
    return full, lead_hits


__all__ = ["run_stage4_gwas"]
