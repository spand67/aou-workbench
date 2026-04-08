"""Stage 2: pathogenic and likely pathogenic panel analysis."""

from __future__ import annotations

from typing import Any

import pandas as pd

from .annotations import annotate_variant_masks
from .config import ProjectConfig
from .io_utils import read_table, write_dataframe, write_json
from .paths import ProjectPaths
from .reporting import write_stage_report
from .statistics import bh_fdr, run_binary_logistic_regression, summarize_binary_exposure


def _prepare_stage2_variants(config: ProjectConfig, matched_df: pd.DataFrame) -> pd.DataFrame:
    stage = config.analysis.stage2
    if stage is None:
        return pd.DataFrame()
    raw = read_table(stage.variant_table).copy()
    raw["person_id"] = raw[stage.person_id_column].astype(str)
    raw["variant_id"] = raw[stage.variant_id_column].astype(str)
    raw["gene"] = raw[stage.gene_column].astype(str)
    raw["dosage"] = pd.to_numeric(raw[stage.dosage_column], errors="coerce").fillna(0.0)
    raw = raw[raw["person_id"].isin(set(matched_df["person_id"]))].copy()
    raw = raw[raw["gene"].isin(set(config.panel.genes_of_interest))].copy()
    return annotate_variant_masks(
        raw,
        clinvar_column=stage.clinvar_column,
        consequence_column=stage.consequence_column,
        revel_column=stage.revel_column,
        af_column=stage.af_column,
        max_af=stage.max_af,
        revel_min=stage.revel_min,
        plof_terms=stage.plof_terms,
        clinvar_plp_terms=stage.clinvar_plp_terms,
    )


def run_stage2_plp_panel(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    stage = config.analysis.stage2
    if stage is None:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    annotated = _prepare_stage2_variants(config, matched_df)
    masked = annotated[annotated["mask_primary"]].copy()

    variant_rows: list[dict[str, Any]] = []
    for variant_id, chunk in masked.groupby("variant_id"):
        exposure = chunk.groupby("person_id")["dosage"].max()
        counts = summarize_binary_exposure(
            matched_df,
            exposure,
            outcome_column=config.analysis.matched_outcome_column,
        )
        regression = run_binary_logistic_regression(
            matched_df,
            exposure,
            outcome_column=config.analysis.matched_outcome_column,
            covariates=stage.covariates,
        )
        meta = chunk.iloc[0]
        variant_rows.append(
            {
                "variant_id": variant_id,
                "gene": meta["gene"],
                "clinvar_significance": meta.get(stage.clinvar_column),
                "consequence": meta.get(stage.consequence_column),
                **counts,
                **regression,
            }
        )

    gene_rows: list[dict[str, Any]] = []
    for gene, chunk in masked.groupby("gene"):
        exposure = chunk.groupby("person_id")["dosage"].max()
        counts = summarize_binary_exposure(
            matched_df,
            exposure,
            outcome_column=config.analysis.matched_outcome_column,
        )
        regression = run_binary_logistic_regression(
            matched_df,
            exposure,
            outcome_column=config.analysis.matched_outcome_column,
            covariates=stage.covariates,
        )
        gene_rows.append({"gene": gene, **counts, **regression})

    person_summary = (
        masked.groupby("person_id", as_index=False)
        .agg(
            n_plp_variants=("variant_id", "nunique"),
            n_plp_genes=("gene", "nunique"),
            max_dosage=("dosage", "max"),
        )
        .sort_values(["n_plp_variants", "n_plp_genes"], ascending=False)
    )

    variant_df = pd.DataFrame(variant_rows)
    gene_df = pd.DataFrame(gene_rows)
    if not variant_df.empty:
        variant_df = variant_df.sort_values(["fisher_p", "regression_p", "variant_id"]).reset_index(drop=True)
    if not gene_df.empty:
        gene_df = gene_df.sort_values(["fisher_p", "regression_p", "gene"]).reset_index(drop=True)
    if not variant_df.empty:
        variant_df["fdr_q"] = bh_fdr(variant_df["fisher_p"].values)
    if not gene_df.empty:
        gene_df["fdr_q"] = bh_fdr(gene_df["fisher_p"].values)

    write_dataframe(variant_df, paths.stage2_variant_tsv)
    write_dataframe(gene_df, paths.stage2_gene_tsv)
    write_dataframe(person_summary, paths.stage2_person_tsv)
    write_json(
        {
            "n_rows_input": int(len(annotated)),
            "n_rows_masked": int(len(masked)),
            "n_panel_genes": int(len(set(config.panel.genes_of_interest))),
            "n_people_with_hits": int(person_summary["person_id"].nunique()) if not person_summary.empty else 0,
        },
        paths.stage2_qc_json,
    )
    write_stage_report(
        title="Stage 2: P/LP panel summary",
        summary_lines=[
            f"- Genes of interest: {len(config.panel.genes_of_interest)}",
            f"- Variants passing the primary mask: {variant_df.shape[0]}",
            f"- Participants with panel hits: {person_summary.shape[0]}",
        ],
        preview_df=gene_df if not gene_df.empty else variant_df,
        preview_columns=["gene", "case_carriers", "control_carriers", "fisher_p", "regression_p"],
        path=paths.stage2_report_md,
    )
    return variant_df, gene_df, person_summary


__all__ = ["run_stage2_plp_panel"]
