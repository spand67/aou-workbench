"""Stage 2: pathogenic and likely pathogenic panel analysis."""

from __future__ import annotations

from typing import Any

import pandas as pd

from .annotations import annotate_variant_masks
from .config import ProjectConfig
from .io_utils import read_table, write_dataframe, write_json
from .paths import ProjectPaths
from .reporting import write_stage_report
from .stage1_prior_variants import _stage1_comparisons
from .stage2_prepare import stage2_sample_manifest_path
from .statistics import bh_fdr, run_binary_logistic_regression, summarize_binary_exposure


def _restrict_to_stage2_samples(sample_df: pd.DataFrame, config: ProjectConfig) -> pd.DataFrame:
    stage = config.analysis.stage2
    if stage is None:
        return sample_df
    manifest_path = stage2_sample_manifest_path(stage.variant_table)
    try:
        present = read_table(manifest_path)
    except Exception as exc:
        raise RuntimeError(
            f"Missing Stage 2 sample manifest at {manifest_path}. Rerun `aou-workbench prepare-stage2` before `run-stage2`."
        ) from exc
    if "person_id" not in present.columns or present.empty:
        raise RuntimeError(
            f"Invalid Stage 2 sample manifest at {manifest_path}. Rerun `aou-workbench prepare-stage2` before `run-stage2`."
        )
    present_ids = set(present["person_id"].astype(str))
    subset = sample_df.copy()
    subset["person_id"] = subset["person_id"].astype(str)
    return subset[subset["person_id"].isin(present_ids)].copy()


def _prepare_stage2_variants(config: ProjectConfig, sample_df: pd.DataFrame) -> pd.DataFrame:
    stage = config.analysis.stage2
    if stage is None:
        return pd.DataFrame()
    sample_df = _restrict_to_stage2_samples(sample_df, config)
    raw = read_table(stage.variant_table).copy()
    raw["person_id"] = raw[stage.person_id_column].astype(str)
    raw["variant_id"] = raw[stage.variant_id_column].astype(str)
    raw["gene"] = raw[stage.gene_column].astype(str)
    raw["dosage"] = pd.to_numeric(raw[stage.dosage_column], errors="coerce").fillna(0.0)
    raw = raw[raw["person_id"].isin(set(sample_df["person_id"]))].copy()
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
    sample_df: pd.DataFrame,
    paths: ProjectPaths,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    stage = config.analysis.stage2
    if stage is None:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    sample_df = _restrict_to_stage2_samples(sample_df, config)
    annotated = _prepare_stage2_variants(config, sample_df)
    masked = annotated[annotated["mask_primary"]].copy()

    variant_rows: list[dict[str, Any]] = []
    gene_rows: list[dict[str, Any]] = []
    comparison_qc: list[dict[str, Any]] = []
    for comparison in _stage1_comparisons(sample_df, config):
        comparison_df = comparison["sample_df"]
        outcome_column = comparison["outcome_column"]
        if outcome_column not in comparison_df.columns:
            continue
        comparison_qc.append(
            {
                "comparison": comparison["comparison"],
                "comparison_label": comparison["comparison_label"],
                "outcome_column": outcome_column,
                "n_cases": int((comparison_df[outcome_column] == 1).sum()),
                "n_controls": int((comparison_df[outcome_column] == 0).sum()),
            }
        )
        for variant_id, chunk in masked.groupby("variant_id"):
            exposure = chunk.groupby("person_id")["dosage"].max()
            counts = summarize_binary_exposure(
                comparison_df,
                exposure,
                outcome_column=outcome_column,
            )
            regression = run_binary_logistic_regression(
                comparison_df,
                exposure,
                outcome_column=outcome_column,
                covariates=stage.covariates,
            )
            meta = chunk.iloc[0]
            variant_rows.append(
                {
                    "comparison": comparison["comparison"],
                    "comparison_label": comparison["comparison_label"],
                    "outcome_column": outcome_column,
                    "variant_id": variant_id,
                    "gene": meta["gene"],
                    "clinvar_significance": meta.get(stage.clinvar_column),
                    "consequence": meta.get(stage.consequence_column),
                    **counts,
                    **regression,
                }
            )

        for gene, chunk in masked.groupby("gene"):
            exposure = chunk.groupby("person_id")["dosage"].max()
            counts = summarize_binary_exposure(
                comparison_df,
                exposure,
                outcome_column=outcome_column,
            )
            regression = run_binary_logistic_regression(
                comparison_df,
                exposure,
                outcome_column=outcome_column,
                covariates=stage.covariates,
            )
            gene_rows.append(
                {
                    "comparison": comparison["comparison"],
                    "comparison_label": comparison["comparison_label"],
                    "outcome_column": outcome_column,
                    "gene": gene,
                    **counts,
                    **regression,
                }
            )

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
        variant_df = variant_df.sort_values(["comparison", "fisher_p", "regression_p", "variant_id"]).reset_index(drop=True)
    if not gene_df.empty:
        gene_df = gene_df.sort_values(["comparison", "fisher_p", "regression_p", "gene"]).reset_index(drop=True)
    if not variant_df.empty:
        variant_df["fdr_q"] = pd.NA
        for comparison_name, comparison_rows in variant_df.groupby("comparison").groups.items():
            indexer = list(comparison_rows)
            variant_df.loc[indexer, "fdr_q"] = bh_fdr(variant_df.loc[indexer, "fisher_p"].values)
    if not gene_df.empty:
        gene_df["fdr_q"] = pd.NA
        for comparison_name, comparison_rows in gene_df.groupby("comparison").groups.items():
            indexer = list(comparison_rows)
            gene_df.loc[indexer, "fdr_q"] = bh_fdr(gene_df.loc[indexer, "fisher_p"].values)

    write_dataframe(variant_df, paths.stage2_variant_tsv)
    write_dataframe(gene_df, paths.stage2_gene_tsv)
    write_dataframe(person_summary, paths.stage2_person_tsv)
    write_json(
        {
            "n_rows_input": int(len(annotated)),
            "n_rows_masked": int(len(masked)),
            "n_panel_genes": int(len(set(config.panel.genes_of_interest))),
            "n_people_with_hits": int(person_summary["person_id"].nunique()) if not person_summary.empty else 0,
            "comparisons": comparison_qc,
        },
        paths.stage2_qc_json,
    )
    summary_lines = [
        f"- Genes of interest: {len(config.panel.genes_of_interest)}",
        f"- Variants passing the primary mask: {variant_df.shape[0]}",
        f"- Participants with panel hits: {person_summary.shape[0]}",
    ]
    for item in comparison_qc:
        summary_lines.append(
            f"- {item['comparison_label']}: {item['n_cases']} cases vs {item['n_controls']} non-rhabdo controls"
        )
    write_stage_report(
        title="Stage 2: P/LP panel summary",
        summary_lines=summary_lines,
        preview_df=gene_df if not gene_df.empty else variant_df,
        preview_columns=["comparison_label", "gene", "case_carriers", "control_carriers", "fisher_p", "regression_p"],
        path=paths.stage2_report_md,
    )
    return variant_df, gene_df, person_summary


__all__ = ["run_stage2_plp_panel", "_restrict_to_stage2_samples"]
