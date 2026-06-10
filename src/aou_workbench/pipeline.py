"""Top-level orchestration for the AoU rhabdomyolysis workflow."""

from __future__ import annotations

from typing import Any

import pandas as pd

from .clinical_model import (
    clinical_model_calibration_path,
    clinical_model_calibration_svg_path,
    clinical_model_coefficients_path,
    clinical_model_cv_metrics_path,
    clinical_model_metrics_path,
    clinical_model_pr_svg_path,
    clinical_model_roc_svg_path,
)
from .cohort import apply_time_anchored_clinical_cofactors, build_rhabdo_cohort, cohort_qc_summary
from .cohort_summary import (
    case_cofactor_prior_timing_histogram_path,
    case_cofactor_prior_timing_path,
    consort_counts_path,
    critical_illness_summary_path,
    matched_table1_path,
    missingness_summary_path,
    model_eligibility_summary_path,
    model_split_summary_path,
    split_table1_path,
)
from .config import ProjectConfig
from .io_utils import write_dataframe, write_json, write_text
from .matching import match_case_controls, matching_qc_summary, matching_universe
from .paths import ProjectPaths, build_output_paths
from .phenotype_sql import (
    render_baseline_sql,
    render_case_tier_sql,
    render_clinical_cofactor_events_sql,
    render_clinical_cofactors_sql,
    render_covariate_sql,
)
from .preflight import apply_runtime_defaults, assert_preflight_ok, run_preflight_checks
from .preindex_profile import (
    preindex_biomarker_path,
    preindex_condition_top_path,
    preindex_measurement_top_path,
    preindex_summary_path,
)
from .reporting import load_table_if_exists, write_dashboard_report, write_final_report
from .sample_restriction import restrict_frame_for_gwas
from .stage1_prepare import prepare_stage1_variant_table
from .stage1_prior_variants import run_stage1_prior_variants
from .stage2_prepare import prepare_stage2_variant_table
from .stage2_plp_panel import run_stage2_plp_panel
from .stage4_prepare import prepare_stage4_acaf_subset
from .stage3_burden import run_stage3_burden
from .stage4_gwas import run_stage4_gwas


def _report_tables(paths: ProjectPaths) -> dict[str, dict[str, pd.DataFrame]]:
    return {
        "cohort_tables": {
            "consort": load_table_if_exists(consort_counts_path(paths)),
            "table1": load_table_if_exists(matched_table1_path(paths)),
            "split_table1": load_table_if_exists(split_table1_path(paths)),
            "split_summary": load_table_if_exists(model_split_summary_path(paths)),
            "eligibility": load_table_if_exists(model_eligibility_summary_path(paths)),
            "critical_illness": load_table_if_exists(critical_illness_summary_path(paths)),
            "case_cofactor_prior_timing": load_table_if_exists(case_cofactor_prior_timing_path(paths)),
            "missingness": load_table_if_exists(missingness_summary_path(paths)),
        },
        "clinical_model_tables": {
            "metrics": load_table_if_exists(clinical_model_metrics_path(paths)),
            "cv_metrics": load_table_if_exists(clinical_model_cv_metrics_path(paths)),
            "coefficients": load_table_if_exists(clinical_model_coefficients_path(paths)),
            "calibration": load_table_if_exists(clinical_model_calibration_path(paths)),
        },
        "preindex_tables": {
            "summary": load_table_if_exists(preindex_summary_path(paths)),
            "biomarkers": load_table_if_exists(preindex_biomarker_path(paths)),
            "top_conditions": load_table_if_exists(preindex_condition_top_path(paths)),
            "top_measurements": load_table_if_exists(preindex_measurement_top_path(paths)),
        },
    }


def _report_source_paths(paths: ProjectPaths) -> dict[str, str]:
    return {
        "consort": consort_counts_path(paths),
        "table1": matched_table1_path(paths),
        "split_table1": split_table1_path(paths),
        "split_summary": model_split_summary_path(paths),
        "eligibility": model_eligibility_summary_path(paths),
        "critical_illness": critical_illness_summary_path(paths),
        "case_cofactor_prior_timing": case_cofactor_prior_timing_path(paths),
        "missingness": missingness_summary_path(paths),
        "preindex_summary": preindex_summary_path(paths),
        "preindex_biomarkers": preindex_biomarker_path(paths),
        "preindex_top_conditions": preindex_condition_top_path(paths),
        "preindex_top_measurements": preindex_measurement_top_path(paths),
        "clinical_model_metrics": clinical_model_metrics_path(paths),
        "clinical_model_cv_metrics": clinical_model_cv_metrics_path(paths),
        "clinical_model_coefficients": clinical_model_coefficients_path(paths),
        "clinical_model_calibration_table": clinical_model_calibration_path(paths),
        "stage1": paths.stage1_results_tsv,
        "stage2_genes": paths.stage2_gene_tsv,
        "stage3": paths.stage3_results_tsv,
        "stage4_hits": paths.stage4_lead_hits_tsv,
    }


def _report_figure_paths(paths: ProjectPaths) -> dict[str, str]:
    return {
        "case_cofactor_prior_timing_histogram": case_cofactor_prior_timing_histogram_path(paths),
        "clinical_model_roc": clinical_model_roc_svg_path(paths),
        "clinical_model_pr": clinical_model_pr_svg_path(paths),
        "clinical_model_calibration": clinical_model_calibration_svg_path(paths),
        "stage4_manhattan": paths.stage4_manhattan_svg,
        "stage4_qq": paths.stage4_qq_svg,
    }


def _write_existing_final_report(effective: ProjectConfig, paths: ProjectPaths) -> None:
    report_tables = _report_tables(paths)
    genetics_tables = {
        "stage1": load_table_if_exists(paths.stage1_results_tsv),
        "stage2_genes": load_table_if_exists(paths.stage2_gene_tsv),
        "stage3": load_table_if_exists(paths.stage3_results_tsv),
        "stage4_hits": load_table_if_exists(paths.stage4_lead_hits_tsv),
    }
    source_paths = _report_source_paths(paths)
    figure_paths = _report_figure_paths(paths)
    write_final_report(
        analysis_name=effective.analysis.analysis_name,
        output_root=paths.run_root,
        stage1=genetics_tables["stage1"],
        stage2_genes=genetics_tables["stage2_genes"],
        stage3=genetics_tables["stage3"],
        stage4_hits=genetics_tables["stage4_hits"],
        path=paths.final_report_md,
        cohort_tables=report_tables["cohort_tables"],
        clinical_model_tables=report_tables["clinical_model_tables"],
        preindex_tables=report_tables["preindex_tables"],
        source_paths=source_paths,
        figure_paths=figure_paths,
    )
    write_dashboard_report(
        analysis_name=effective.analysis.analysis_name,
        output_root=paths.run_root,
        path=paths.final_dashboard_html,
        cohort_tables=report_tables["cohort_tables"],
        clinical_model_tables=report_tables["clinical_model_tables"],
        preindex_tables=report_tables["preindex_tables"],
        genetics_tables=genetics_tables,
        source_paths=source_paths,
        figure_paths=figure_paths,
    )


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
    write_text(render_case_tier_sql(config, config.phenotype.broad), f"{paths.cohort_sql_root}/rhabdo_broad.sql")
    if config.phenotype.clinical_cofactors:
        write_text(render_clinical_cofactors_sql(config), f"{paths.cohort_sql_root}/clinical_cofactors.sql")
        write_text(render_clinical_cofactor_events_sql(config), f"{paths.cohort_sql_root}/clinical_cofactor_events.sql")
    write_text(render_covariate_sql(config), f"{paths.cohort_sql_root}/rhabdo_covariates.sql")


def _restrict_cohort_to_wgs_if_requested(
    config: ProjectConfig,
    cohort_df: pd.DataFrame,
    *,
    require_wgs: bool,
) -> pd.DataFrame:
    if not require_wgs:
        return cohort_df
    return restrict_frame_for_gwas(config, cohort_df, require_wgs=True)


def build_cohort_artifacts(config: ProjectConfig, *, require_wgs: bool = False) -> tuple[ProjectConfig, ProjectPaths, pd.DataFrame]:
    effective = apply_runtime_defaults(config)
    paths = build_output_paths(effective)
    cohort_df = build_rhabdo_cohort(effective)
    cohort_df = _restrict_cohort_to_wgs_if_requested(effective, cohort_df, require_wgs=require_wgs)
    write_dataframe(cohort_df, paths.built_cohort_tsv)
    _write_rendered_sql(effective, paths)
    write_json(cohort_qc_summary(cohort_df), paths.cohort_qc_json)
    _write_manifest(effective, paths, extra={"built_cohort_rows": int(len(cohort_df)), "require_wgs": bool(require_wgs)})
    return effective, paths, cohort_df


def match_controls_artifacts(
    config: ProjectConfig,
    cohort_df: pd.DataFrame | None = None,
    *,
    require_wgs: bool = False,
) -> tuple[ProjectConfig, ProjectPaths, pd.DataFrame]:
    effective = apply_runtime_defaults(config)
    paths = build_output_paths(effective)
    if cohort_df is None:
        cohort_df = build_rhabdo_cohort(effective)
    cohort_df = _restrict_cohort_to_wgs_if_requested(effective, cohort_df, require_wgs=require_wgs)
    matching_universe_df = matching_universe(cohort_df, effective)
    matched_df = apply_time_anchored_clinical_cofactors(effective, match_case_controls(cohort_df, effective))
    write_dataframe(matched_df, paths.matched_cohort_tsv)
    qc_payload = cohort_qc_summary(cohort_df)
    qc_payload["matching_universe_rows"] = int(len(matching_universe_df))
    qc_payload["matching_universe_people"] = int(matching_universe_df["person_id"].astype(str).nunique())
    qc_payload.update(matching_qc_summary(matched_df))
    write_json(qc_payload, paths.cohort_qc_json)
    _write_manifest(effective, paths, extra={"matched_rows": int(len(matched_df)), "require_wgs": bool(require_wgs)})
    return effective, paths, matched_df


def run_all(config: ProjectConfig, *, skip_preflight: bool = False) -> ProjectPaths:
    effective = apply_runtime_defaults(config)
    paths = build_output_paths(effective)
    if not skip_preflight:
        checks = run_preflight_checks(effective)
        assert_preflight_ok(checks)
    cohort_df = build_rhabdo_cohort(effective)
    matching_universe_df = matching_universe(cohort_df, effective)
    matched_df = apply_time_anchored_clinical_cofactors(effective, match_case_controls(cohort_df, effective))
    write_dataframe(cohort_df, paths.built_cohort_tsv)
    write_dataframe(matched_df, paths.matched_cohort_tsv)
    _write_rendered_sql(effective, paths)
    qc_payload = cohort_qc_summary(cohort_df)
    qc_payload["matching_universe_rows"] = int(len(matching_universe_df))
    qc_payload["matching_universe_people"] = int(matching_universe_df["person_id"].astype(str).nunique())
    qc_payload.update(matching_qc_summary(matched_df))
    write_json(qc_payload, paths.cohort_qc_json)

    stage1_df = pd.DataFrame()
    stage2_variant_df = pd.DataFrame()
    stage2_gene_df = pd.DataFrame()
    stage3_df = pd.DataFrame()
    stage4_lead_hits = pd.DataFrame()
    if effective.analysis.run_stage1:
        prepare_stage1_variant_table(effective, cohort_df)
        stage1_df = run_stage1_prior_variants(effective, cohort_df, paths)
    if effective.analysis.run_stage2:
        prepare_stage2_variant_table(effective, cohort_df)
        stage2_variant_df, stage2_gene_df, _ = run_stage2_plp_panel(effective, cohort_df, paths)
    if effective.analysis.run_stage3:
        stage3_df = run_stage3_burden(effective, matched_df, paths)
    if effective.analysis.run_stage4:
        prepare_stage4_acaf_subset(effective, matched_df, paths)
        _, stage4_lead_hits = run_stage4_gwas(effective, matched_df, paths)

    _write_existing_final_report(effective, paths)
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
    _write_existing_final_report(effective, paths)
    _write_manifest(effective, paths, extra={"report_only": True})
    return paths


__all__ = [
    "build_cohort_artifacts",
    "match_controls_artifacts",
    "render_existing_report",
    "run_all",
]
