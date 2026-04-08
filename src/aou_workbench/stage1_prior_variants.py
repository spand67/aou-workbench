"""Stage 1: a priori candidate variant analysis."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any

import numpy as np
import pandas as pd

from .config import ProjectConfig
from .io_utils import read_table, write_dataframe, write_json
from .paths import ProjectPaths
from .reporting import write_stage_report
from .statistics import bh_fdr, run_binary_logistic_regression, summarize_binary_exposure


def run_stage1_prior_variants(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
) -> pd.DataFrame:
    stage = config.analysis.stage1
    if stage is None:
        return pd.DataFrame()
    raw = read_table(stage.variant_table).copy()
    raw["person_id"] = raw[stage.person_id_column].astype(str)
    raw["variant_id"] = raw[stage.variant_id_column].astype(str)
    raw["dosage"] = pd.to_numeric(raw[stage.dosage_column], errors="coerce").fillna(0.0)
    raw["gene"] = raw.get(stage.gene_column, pd.Series([""] * len(raw), index=raw.index)).astype(str)

    panel_variants = list(config.panel.a_priori_variants)
    if panel_variants:
        requested = pd.DataFrame(
            [
                {
                    "variant_id": variant.variant_id,
                    "label": variant.label,
                    "gene": variant.gene,
                    "rsid": variant.rsid,
                    "source": variant.source,
                    "evidence_tier": variant.evidence_tier,
                    "exact_test_model": variant.exact_test_model,
                }
                for variant in panel_variants
            ]
        )
    else:
        requested = raw[["variant_id", "gene"]].drop_duplicates().copy()
        requested["label"] = requested["variant_id"]
        requested["rsid"] = None
        requested["source"] = None
        requested["evidence_tier"] = None
        requested["exact_test_model"] = "carrier_vs_noncarrier"

    rows: list[dict[str, Any]] = []
    for variant in requested.itertuples(index=False):
        subset = raw[raw["variant_id"] == variant.variant_id]
        exposure = subset.groupby("person_id")["dosage"].max() if not subset.empty else pd.Series(dtype=float)
        counts = summarize_binary_exposure(
            matched_df,
            exposure,
            outcome_column=config.analysis.matched_outcome_column,
            carrier_threshold=stage.carrier_min_dosage,
        )
        regression = run_binary_logistic_regression(
            matched_df,
            exposure,
            outcome_column=config.analysis.matched_outcome_column,
            covariates=stage.covariates,
        )
        rows.append(
            {
                "variant_id": variant.variant_id,
                "label": variant.label,
                "gene": variant.gene,
                "rsid": variant.rsid,
                "source": variant.source,
                "evidence_tier": variant.evidence_tier,
                **counts,
                **regression,
            }
        )
    result = pd.DataFrame(rows).sort_values(["fisher_p", "regression_p", "variant_id"]).reset_index(drop=True)
    if not result.empty:
        result["fdr_q"] = bh_fdr(result["fisher_p"].values)
        result["regression_fdr_q"] = bh_fdr(result["regression_p"].values)
    write_dataframe(result, paths.stage1_results_tsv)
    write_json(
        {
            "n_variants_requested": int(len(requested)),
            "n_variants_with_rows": int(raw["variant_id"].nunique()),
            "n_cases": int((matched_df[config.analysis.matched_outcome_column] == 1).sum()),
            "n_controls": int((matched_df[config.analysis.matched_outcome_column] == 0).sum()),
        },
        paths.stage1_qc_json,
    )
    write_stage_report(
        title="Stage 1: A priori candidate variants",
        summary_lines=[
            f"- Requested variants: {len(requested)}",
            f"- Matched cases: {(matched_df[config.analysis.matched_outcome_column] == 1).sum()}",
            f"- Matched controls: {(matched_df[config.analysis.matched_outcome_column] == 0).sum()}",
        ],
        preview_df=result,
        preview_columns=["label", "gene", "case_carriers", "control_carriers", "fisher_p", "regression_p"],
        path=paths.stage1_report_md,
    )
    return result


__all__ = ["run_stage1_prior_variants"]
