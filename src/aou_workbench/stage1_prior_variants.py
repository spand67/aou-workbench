"""Stage 1: a priori candidate variant analysis."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .config import ProjectConfig
from .io_utils import read_table, write_dataframe, write_json
from .paths import ProjectPaths
from .reporting import write_stage_report
from .stage1_prepare import stage1_sample_manifest_path
from .statistics import bh_fdr, run_binary_logistic_regression, summarize_binary_exposure


def _variant_exposure(
    subset: pd.DataFrame,
    *,
    exact_test_model: str,
) -> tuple[pd.Series, float]:
    dosage = subset.groupby("person_id")["dosage"].max() if not subset.empty else pd.Series(dtype=float)
    if exact_test_model == "hom_alt_vs_rest":
        return (dosage >= 2.0).astype(float), 1.0
    return dosage, 1.0


def _stage1_comparisons(sample_df: pd.DataFrame, config: ProjectConfig) -> list[dict[str, Any]]:
    if {"rhabdo_case", "rhabdo_primary_case"}.issubset(sample_df.columns):
        all_non_rhabdo = sample_df["rhabdo_case"].fillna(0).astype(int) == 0
        comparisons = [
            {
                "comparison": "omop_rhabdo_vs_non_rhabdo",
                "comparison_label": "OMOP rhabdo vs non-rhabdo",
                "outcome_column": "rhabdo_case",
                "sample_df": sample_df[sample_df["rhabdo_case"].fillna(0).astype(int).isin([0, 1])].copy(),
            },
            {
                "comparison": "omop_rhabdo_plus_ck_vs_non_rhabdo",
                "comparison_label": "OMOP rhabdo + CK>=5000 vs non-rhabdo",
                "outcome_column": "rhabdo_primary_case",
                "sample_df": sample_df[(sample_df["rhabdo_primary_case"].fillna(0).astype(int) == 1) | all_non_rhabdo].copy(),
            },
        ]
        if "ancestry_pred" in sample_df.columns:
            ancestry_values = sorted(
                str(value)
                for value in sample_df["ancestry_pred"].dropna().astype(str).unique().tolist()
                if str(value).strip()
            )
            for ancestry in ancestry_values:
                ancestry_df = sample_df[sample_df["ancestry_pred"].astype(str) == ancestry].copy()
                if ancestry_df.empty:
                    continue
                if ancestry_df["rhabdo_case"].fillna(0).astype(int).nunique() < 2:
                    continue
                comparisons.append(
                    {
                        "comparison": f"omop_rhabdo_vs_non_rhabdo_ancestry_{ancestry.lower()}",
                        "comparison_label": f"OMOP rhabdo vs non-rhabdo ({ancestry})",
                        "outcome_column": "rhabdo_case",
                        "sample_df": ancestry_df,
                    }
                )
        return comparisons
    return [
        {
            "comparison": "matched_case_control",
            "comparison_label": "Matched case-control",
            "outcome_column": config.analysis.matched_outcome_column,
            "sample_df": sample_df.copy(),
        }
    ]


def _restrict_to_stage1_wgs_samples(sample_df: pd.DataFrame, config: ProjectConfig) -> pd.DataFrame:
    stage = config.analysis.stage1
    if stage is None:
        return sample_df
    manifest_path = stage1_sample_manifest_path(stage.variant_table)
    try:
        present = read_table(manifest_path)
    except Exception as exc:
        raise RuntimeError(
            f"Missing Stage 1 WGS sample manifest at {manifest_path}. Rerun `aou-workbench prepare-stage1` before `run-stage1`."
        ) from exc
    if "person_id" not in present.columns or present.empty:
        raise RuntimeError(
            f"Invalid Stage 1 WGS sample manifest at {manifest_path}. Rerun `aou-workbench prepare-stage1` before `run-stage1`."
        )
    present_ids = set(present["person_id"].astype(str))
    subset = sample_df.copy()
    subset["person_id"] = subset["person_id"].astype(str)
    return subset[subset["person_id"].isin(present_ids)].copy()


def run_stage1_prior_variants(
    config: ProjectConfig,
    sample_df: pd.DataFrame,
    paths: ProjectPaths,
) -> pd.DataFrame:
    stage = config.analysis.stage1
    if stage is None:
        return pd.DataFrame()
    sample_df = _restrict_to_stage1_wgs_samples(sample_df, config)
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
        for variant in requested.itertuples(index=False):
            subset = raw[raw["variant_id"] == variant.variant_id]
            exposure, carrier_threshold = _variant_exposure(
                subset,
                exact_test_model=variant.exact_test_model,
            )
            counts = summarize_binary_exposure(
                comparison_df,
                exposure,
                outcome_column=outcome_column,
                carrier_threshold=carrier_threshold,
            )
            regression = run_binary_logistic_regression(
                comparison_df,
                exposure,
                outcome_column=outcome_column,
                covariates=stage.covariates,
            )
            rows.append(
                {
                    "comparison": comparison["comparison"],
                    "comparison_label": comparison["comparison_label"],
                    "outcome_column": outcome_column,
                    "variant_id": variant.variant_id,
                    "label": variant.label,
                    "gene": variant.gene,
                    "rsid": variant.rsid,
                    "source": variant.source,
                    "evidence_tier": variant.evidence_tier,
                    "exact_test_model": variant.exact_test_model,
                    **counts,
                    **regression,
                }
            )
    result = pd.DataFrame(rows).sort_values(["comparison", "fisher_p", "regression_p", "variant_id"]).reset_index(drop=True)
    if not result.empty:
        result["fdr_q"] = np.nan
        result["regression_fdr_q"] = np.nan
        for comparison_name, comparison_rows in result.groupby("comparison").groups.items():
            indexer = list(comparison_rows)
            result.loc[indexer, "fdr_q"] = bh_fdr(result.loc[indexer, "fisher_p"].values)
            result.loc[indexer, "regression_fdr_q"] = bh_fdr(result.loc[indexer, "regression_p"].values)
    write_dataframe(result, paths.stage1_results_tsv)
    write_json(
        {
            "n_variants_requested": int(len(requested)),
            "n_variants_with_rows": int(raw["variant_id"].nunique()),
            "comparisons": comparison_qc,
        },
        paths.stage1_qc_json,
    )
    summary_lines = [f"- Requested variants: {len(requested)}"]
    for item in comparison_qc:
        summary_lines.append(
            f"- {item['comparison_label']}: {item['n_cases']} cases vs {item['n_controls']} non-rhabdo controls"
        )
    write_stage_report(
        title="Stage 1: A priori candidate variants",
        summary_lines=summary_lines,
        preview_df=result,
        preview_columns=["comparison_label", "label", "gene", "case_carriers", "control_carriers", "fisher_p", "regression_p"],
        path=paths.stage1_report_md,
    )
    return result


__all__ = [
    "run_stage1_prior_variants",
    "_restrict_to_stage1_wgs_samples",
    "_stage1_comparisons",
    "_variant_exposure",
]
