"""Stage 3: gene burden analysis."""

from __future__ import annotations

from typing import Any

import pandas as pd

from .annotations import annotate_variant_masks
from .config import ProjectConfig
from .io_utils import read_table, write_dataframe, write_json
from .paths import ProjectPaths
from .reporting import write_stage_report
from .statistics import bh_fdr, run_binary_logistic_regression, summarize_binary_exposure


def _mask_series(df: pd.DataFrame, mask_name: str) -> pd.Series:
    lookup = {
        "pLoF": df["mask_plof"],
        "ClinVar_PLP": df["mask_clinvar_plp"],
        "pLoF_or_REVEL_0_8": df["mask_plof_or_revel"],
    }
    return lookup.get(mask_name, df["mask_primary"])


def run_stage3_burden(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
) -> pd.DataFrame:
    stage = config.analysis.stage3
    if stage is None:
        return pd.DataFrame()
    raw = read_table(stage.variant_table).copy()
    raw["person_id"] = raw[stage.person_id_column].astype(str)
    raw["variant_id"] = raw[stage.variant_id_column].astype(str)
    raw["gene"] = raw[stage.gene_column].astype(str)
    raw["dosage"] = pd.to_numeric(raw[stage.dosage_column], errors="coerce").fillna(0.0)
    raw = raw[raw["person_id"].isin(set(matched_df["person_id"]))].copy()
    if stage.mode == "targeted" or not stage.run_exome_wide:
        raw = raw[raw["gene"].isin(set(config.panel.burden_target_genes))].copy()
    annotated = annotate_variant_masks(
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

    rows: list[dict[str, Any]] = []
    for mask_name in stage.masks:
        mask = _mask_series(annotated, mask_name)
        subset = annotated[mask].copy()
        for gene, chunk in subset.groupby("gene"):
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
            rows.append(
                {
                    "mask": mask_name,
                    "gene": gene,
                    "n_variants_in_mask": int(chunk["variant_id"].nunique()),
                    **counts,
                    **regression,
                }
            )

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values(["regression_p", "fisher_p", "mask", "gene"]).reset_index(drop=True)
    if not result.empty:
        result["fdr_q"] = bh_fdr(result["regression_p"].values)
    write_dataframe(result, paths.stage3_results_tsv)
    write_json(
        {
            "n_input_rows": int(len(annotated)),
            "n_genes_tested": int(result["gene"].nunique()) if not result.empty else 0,
            "mode": stage.mode,
            "run_exome_wide": bool(stage.run_exome_wide),
        },
        paths.stage3_qc_json,
    )
    write_stage_report(
        title="Stage 3: gene burden analysis",
        summary_lines=[
            f"- Burden masks: {', '.join(stage.masks)}",
            f"- Genes tested: {result['gene'].nunique() if not result.empty else 0}",
            f"- Mode: {stage.mode}",
        ],
        preview_df=result,
        preview_columns=["gene", "mask", "case_carriers", "control_carriers", "regression_p"],
        path=paths.stage3_report_md,
    )
    return result


__all__ = ["run_stage3_burden"]
