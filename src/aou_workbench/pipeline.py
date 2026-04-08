"""Top-level orchestration for the AoU rhabdomyolysis workflow."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any

import pandas as pd

from .cohort import build_rhabdo_cohort, cohort_qc_summary
from .config import ProjectConfig
from .io_utils import write_dataframe, write_json, write_text
from .matching import match_case_controls, matching_qc_summary
from .paths import ProjectPaths, build_output_paths
from .phenotype_sql import (
    render_baseline_sql,
    render_case_tier_sql,
    render_clinical_cofactors_sql,
    render_covariate_sql,
)
from .preflight import apply_runtime_defaults, assert_preflight_ok, run_preflight_checks
from .reporting import load_table_if_exists, write_final_report
from .stage1_prior_variants import run_stage1_prior_variants
from .stage2_plp_panel import run_stage2_plp_panel
from .stage3_burden import run_stage3_burden
from .stage4_gwas import run_stage4_gwas


def _write_manifest(config: ProjectConfig, paths: ProjectPaths, extra: dict[str, Any] | None = None) -> None:
    payload = {
        "analysis_name": config.analysis.analysis_name,
        "config_hash": config.config_hash,
        "config": config.to_dict(),
        "paths": paths.as_dict(),
    }
    if extra:
        payload.update(extra)
    write_json(payload, paths.manifest_json)


def _write_rendered_sql(config: ProjectConfig, paths: ProjectPaths) -> None:
    write_text(render_baseline_sql(config), f"{paths.cohort_sql_root}/baseline.sql")
    write_text(render_case_tier_sql(config, config.phenotype.definite), f"{paths.cohort_sql_root}/rhabdo_definite.sql")
    write_text(render_case_tier_sql(config, config.phenotype.probable), f"{paths.cohort_sql_root}/rhabdo_probable.sql")
    if config.phenotype.clinical_cofactors:
        write_text(render_clinical_cofactors_sql(config), f"{paths.cohort_sql_root}/clinical_cofactors.sql")
    write_text(render_covariate_sql(config), f"{paths.cohort_sql_root}/rhabdo_covariates.sql")


def build_cohort_artifacts(config: ProjectConfig) -> tuple[ProjectConfig, ProjectPaths, pd.DataFrame]:
    effective = apply_runtime_defaults(config)
    paths = build_output_paths(effective)
    cohort_df = build_rhabdo_cohort(effective)
    write_dataframe(cohort_df, paths.built_cohort_tsv)
    _write_rendered_sql(effective, paths)
    write_json(cohort_qc_summary(cohort_df), paths.cohort_qc_json)
    _write_manifest(effective, paths, extra={"built_cohort_rows": int(len(cohort_df))})
    return effective, paths, cohort_df


def match_controls_artifacts(config: ProjectConfig, cohort_df: pd.DataFrame | None = None) -> tuple[ProjectConfig, ProjectPaths, pd.DataFrame]:
    effective = apply_runtime_defaults(config)
    paths = build_output_paths(effective)
    if cohort_df is None:
        cohort_df = build_rhabdo_cohort(effective)
    matched_df = match_case_controls(cohort_df, effective)
    write_dataframe(matched_df, paths.matched_cohort_tsv)
    qc_payload = cohort_qc_summary(cohort_df)
    qc_payload.update(matching_qc_summary(matched_df))
    write_json(qc_payload, paths.cohort_qc_json)
    _write_manifest(effective, paths, extra={"matched_rows": int(len(matched_df))})
    return effective, paths, matched_df


def run_all(config: ProjectConfig, *, skip_preflight: bool = False) -> ProjectPaths:
    effective = apply_runtime_defaults(config)
    paths = build_output_paths(effective)
    if not skip_preflight:
        checks = run_preflight_checks(effective)
        assert_preflight_ok(checks)
    cohort_df = build_rhabdo_cohort(effective)
    matched_df = match_case_controls(cohort_df, effective)
    write_dataframe(cohort_df, paths.built_cohort_tsv)
    write_dataframe(matched_df, paths.matched_cohort_tsv)
    _write_rendered_sql(effective, paths)
    qc_payload = cohort_qc_summary(cohort_df)
    qc_payload.update(matching_qc_summary(matched_df))
    write_json(qc_payload, paths.cohort_qc_json)

    stage1_df = pd.DataFrame()
    stage2_variant_df = pd.DataFrame()
    stage2_gene_df = pd.DataFrame()
    stage3_df = pd.DataFrame()
    stage4_lead_hits = pd.DataFrame()
    if effective.analysis.run_stage1:
        stage1_df = run_stage1_prior_variants(effective, matched_df, paths)
    if effective.analysis.run_stage2:
        stage2_variant_df, stage2_gene_df, _ = run_stage2_plp_panel(effective, matched_df, paths)
    if effective.analysis.run_stage3:
        stage3_df = run_stage3_burden(effective, matched_df, paths)
    if effective.analysis.run_stage4:
        _, stage4_lead_hits = run_stage4_gwas(effective, matched_df, paths)

    write_final_report(
        analysis_name=effective.analysis.analysis_name,
        output_root=paths.run_root,
        stage1=stage1_df,
        stage2_genes=stage2_gene_df if not stage2_gene_df.empty else stage2_variant_df,
        stage3=stage3_df,
        stage4_hits=stage4_lead_hits,
        path=paths.final_report_md,
    )
    _write_manifest(
        effective,
        paths,
        extra={
            "built_cohort_rows": int(len(cohort_df)),
            "matched_rows": int(len(matched_df)),
            "stage1_rows": int(stage1_df.shape[0]),
            "stage2_rows": int(stage2_variant_df.shape[0]),
            "stage3_rows": int(stage3_df.shape[0]),
            "stage4_lead_hits": int(stage4_lead_hits.shape[0]),
        },
    )
    return paths


def render_existing_report(config: ProjectConfig) -> ProjectPaths:
    effective = apply_runtime_defaults(config)
    paths = build_output_paths(effective)
    write_final_report(
        analysis_name=effective.analysis.analysis_name,
        output_root=paths.run_root,
        stage1=load_table_if_exists(paths.stage1_results_tsv),
        stage2_genes=load_table_if_exists(paths.stage2_gene_tsv),
        stage3=load_table_if_exists(paths.stage3_results_tsv),
        stage4_hits=load_table_if_exists(paths.stage4_lead_hits_tsv),
        path=paths.final_report_md,
    )
    _write_manifest(effective, paths, extra={"report_only": True})
    return paths


__all__ = [
    "build_cohort_artifacts",
    "match_controls_artifacts",
    "render_existing_report",
    "run_all",
]
