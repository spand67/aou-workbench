"""CLI entrypoint for the AoU workbench."""

from __future__ import annotations

import argparse
import json
import os
import sys

from .clinical_model import (
    clinical_model_calibration_path,
    clinical_model_calibration_svg_path,
    clinical_model_coefficients_path,
    clinical_model_cv_metrics_path,
    clinical_model_metrics_path,
    clinical_model_pr_svg_path,
    clinical_model_predictions_path,
    clinical_model_report_path,
    clinical_model_roc_svg_path,
    run_clinical_model,
)
from .clinical_prs_model import (
    clinical_prs_model_calibration_path,
    clinical_prs_model_calibration_svg_path,
    clinical_prs_model_coefficients_path,
    clinical_prs_model_default_label,
    clinical_prs_model_metrics_path,
    clinical_prs_model_predictions_path,
    clinical_prs_model_pr_svg_path,
    clinical_prs_model_report_path,
    clinical_prs_model_roc_svg_path,
    run_clinical_prs_model,
)
from .config import load_project_config
from .cohort import apply_time_anchored_clinical_cofactors
from .cohort_summary import (
    case_cofactor_prior_timing_histogram_path,
    case_cofactor_prior_timing_path,
    case_cofactor_prior_timing_report_path,
    characterize_case_control_cohort,
    clinical_model_input_path,
    clinical_characterization_report_path,
    cohort_summary_report_path,
    cohort_summary_table_path,
    consort_counts_report_path,
    consort_counts_path,
    critical_illness_summary_report_path,
    critical_illness_summary_path,
    matched_table1_report_path,
    matched_table1_path,
    model_eligibility_summary_report_path,
    model_eligibility_summary_path,
    missingness_summary_report_path,
    missingness_summary_path,
    model_split_summary_report_path,
    model_split_summary_path,
    split_table1_report_path,
    split_table1_path,
    summarize_clinical_demographics,
)
from .eir import (
    build_eir_cohort_artifacts,
    characterize_eir_artifacts,
    estimate_eir_cohort_artifacts,
    eir_characterization_report_path,
    eir_consort_counts_path,
    eir_missingness_path,
    eir_model_input_path,
    eir_risk_decile_path,
    eir_sparse_report_path,
    eir_sparse_status_path,
    eir_split_summary_path,
    eir_table1_path,
    run_eir_clinical_model,
)
from .gwas_workflow import prepare_terminal_gwas_workspace
from .incident_feasibility import (
    estimate_incident_feasibility_artifacts,
    incident_baseline_history_bins_path,
    incident_case_funnel_path,
    incident_control_funnel_path,
    incident_feasibility_counts_path,
    incident_feasibility_report_path,
    incident_microarray_overlap_counts_path,
    run_incident_feasibility,
)
from .io_utils import read_table
from .microarray_bigsnpr import (
    BIGSNPR_DEFAULT_ALPHAS,
    microarray_bigsnpr_default_label,
    microarray_bigsnpr_metrics_path,
    microarray_bigsnpr_metadata_path,
    microarray_bigsnpr_predictions_path,
    microarray_bigsnpr_qc_path,
    microarray_bigsnpr_report_path,
    microarray_bigsnpr_risk_decile_path,
    microarray_bigsnpr_score_distribution_path,
    microarray_bigsnpr_script_path,
    microarray_bigsnpr_selected_variants_path,
    microarray_bigsnpr_variant_qc_summary_path,
    run_microarray_bigsnpr_model,
)
from .microarray_plink_gwas import (
    microarray_plink_default_label,
    microarray_plink_lead_hits_path,
    microarray_plink_manhattan_path,
    microarray_plink_qc_path,
    microarray_plink_qq_path,
    microarray_plink_report_path,
    microarray_plink_results_path,
    microarray_plink_variant_qc_summary_path,
    run_microarray_plink_gwas,
)
from .microarray_plink_prs import (
    DEFAULT_PRS_THRESHOLDS,
    microarray_prs_case_status_svg_path,
    microarray_prs_default_label,
    microarray_prs_metrics_path,
    microarray_prs_qc_path,
    microarray_prs_range_path,
    microarray_prs_report_path,
    microarray_prs_scores_path,
    microarray_prs_weights_path,
    parse_thresholds,
    run_microarray_plink_prs,
)
from .model_comparison import (
    model_comparison_metrics_path,
    model_comparison_predictions_path,
    model_comparison_report_path,
    run_heldout_model_comparison,
)
from .paths import build_output_paths, project_path
from .pipeline import build_cohort_artifacts, match_controls_artifacts, render_existing_report, run_all
from .preindex_profile import (
    preindex_biomarker_path,
    preindex_condition_top_path,
    preindex_measurement_top_path,
    preindex_report_path,
    preindex_summary_path,
    profile_preindex_case_data,
)
from .presentation_dashboard import render_presentation_dashboard
from .prs_diagnostics import (
    prs_diagnostics_ancestry_path,
    prs_diagnostics_bootstrap_ci_path,
    prs_diagnostics_calibration_path,
    prs_diagnostics_cofactor_path,
    prs_diagnostics_deciles_path,
    prs_diagnostics_definite_path,
    prs_diagnostics_overall_metrics_path,
    prs_diagnostics_qc_path,
    prs_diagnostics_report_path,
    run_prs_diagnostics,
)
from .preflight import apply_runtime_defaults, format_preflight_report, run_preflight_checks
from .regenie import prepare_regenie_inputs
from .stage1_prepare import prepare_stage1_variant_table, prepare_wgs_sample_manifest
from .stage1_prior_variants import run_stage1_prior_variants
from .stage2_prepare import prepare_stage2_variant_table
from .stage2_plp_panel import run_stage2_plp_panel
from .stage4_hail_gwas import (
    hail_pilot_qc_pass_mt_uri,
    hail_pilot_lead_hits_path,
    hail_pilot_manhattan_path,
    hail_pilot_qc_path,
    hail_pilot_default_label,
    hail_pilot_qq_path,
    hail_pilot_report_path,
    hail_pilot_results_ht_uri,
    hail_pilot_results_preview_path,
    hail_pilot_results_path,
    hail_pilot_results_tsv_uri,
    hail_pilot_variant_qc_summary_path,
    hail_stage4_full_results_path,
    hail_stage4_lead_hits_path,
    hail_stage4_manhattan_path,
    hail_stage4_qc_path,
    hail_stage4_qq_path,
    hail_stage4_report_path,
    run_stage4_hail_pilot_gwas,
    run_stage4_hail_gwas,
)
from .stage4_prepare import prepare_stage4_acaf_subset
from .stage3_burden import run_stage3_burden
from .stage4_gwas import run_stage4_gwas


def _add_config_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--workbench-config", default=project_path("configs", "workbench.yaml"))
    parser.add_argument("--phenotype-config", default=project_path("configs", "rhabdo", "phenotype.yaml"))
    parser.add_argument("--cohort-config", default=project_path("configs", "rhabdo", "cohort.yaml"))
    parser.add_argument("--panel-config", default=project_path("configs", "rhabdo", "panel.yaml"))
    parser.add_argument("--analysis-config", default=project_path("configs", "rhabdo", "analysis.yaml"))


def _add_eir_config_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--workbench-config", default=project_path("configs", "workbench.yaml"))
    parser.add_argument("--phenotype-config", default=project_path("configs", "eir", "phenotype.yaml"))
    parser.add_argument("--cohort-config", default=project_path("configs", "eir", "cohort.yaml"))
    parser.add_argument("--panel-config", default=project_path("configs", "eir", "panel.yaml"))
    parser.add_argument("--analysis-config", default=project_path("configs", "eir", "analysis.yaml"))


def _load_config(args: argparse.Namespace):
    return load_project_config(
        workbench_path=args.workbench_config,
        phenotype_path=args.phenotype_config,
        cohort_path=args.cohort_config,
        panel_path=args.panel_config,
        analysis_path=args.analysis_config,
    )


def _load_or_build_matched_artifacts(config, *, require_wgs: bool = False):
    effective = apply_runtime_defaults(config)
    paths = build_output_paths(effective)
    if os.path.exists(paths.matched_cohort_tsv) and _outputs_satisfy_wgs_requirement(paths, require_wgs=require_wgs):
        matched_df = read_table(paths.matched_cohort_tsv)
        matched_df = apply_time_anchored_clinical_cofactors(effective, matched_df)
        return effective, paths, matched_df
    return match_controls_artifacts(config, require_wgs=require_wgs)


def _load_or_build_cohort_artifacts(config, *, require_wgs: bool = False):
    effective = apply_runtime_defaults(config)
    paths = build_output_paths(effective)
    if os.path.exists(paths.built_cohort_tsv) and _outputs_satisfy_wgs_requirement(paths, require_wgs=require_wgs):
        cohort_df = read_table(paths.built_cohort_tsv)
        return effective, paths, cohort_df
    return build_cohort_artifacts(config, require_wgs=require_wgs)


def _manifest_requires_wgs(paths) -> bool:
    try:
        with open(paths.manifest_json, "r", encoding="utf-8") as handle:
            manifest = json.load(handle)
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return False
    return bool(manifest.get("require_wgs"))


def _outputs_satisfy_wgs_requirement(paths, *, require_wgs: bool) -> bool:
    return not require_wgs or _manifest_requires_wgs(paths)


def _clinical_model_input_needs_refresh(paths, *, require_wgs: bool) -> bool:
    if not os.path.exists(clinical_model_input_path(paths)):
        return True
    return require_wgs and not _manifest_requires_wgs(paths)


def _hail_pilot_required_input_columns(config, eligibility_flag: str) -> set[str]:
    return {
        "person_id",
        config.analysis.matched_outcome_column,
        "analysis_split",
        eligibility_flag,
        "match_group_id",
        "case_tier",
        "eligible_control",
        "eligible_ehr_denominator",
        "broad_rhabdo_case",
        "definite_rhabdo_case",
        "high_ck_without_rhabdo",
    }


def _missing_columns(frame, required: set[str]) -> list[str]:
    return sorted(required.difference(frame.columns))


def _load_hail_pilot_matched_input(config, *, eligibility_flag: str):
    effective, paths, cohort_df = _load_or_build_cohort_artifacts(config)
    input_path = clinical_model_input_path(paths)
    required = _hail_pilot_required_input_columns(effective, eligibility_flag)

    needs_refresh = True
    missing = sorted(required)
    matched_df = None
    if os.path.exists(input_path):
        matched_df = read_table(input_path)
        missing = _missing_columns(matched_df, required)
        needs_refresh = bool(missing)

    if needs_refresh:
        missing_text = ", ".join(missing) if missing else "file missing"
        print(f"Refreshing clinical model input for Hail pilot; missing columns: {missing_text}", flush=True)
        _, _, matched_cohort = _load_or_build_matched_artifacts(effective)
        characterize_case_control_cohort(effective, cohort_df, matched_cohort, paths)
        matched_df = read_table(input_path)
        missing = _missing_columns(matched_df, required)
        if missing:
            raise RuntimeError(
                "Existing cohort or matched-cohort outputs are too stale for the Hail pilot. "
                "Still missing required columns after refreshing clinical_model_input.tsv: "
                + ", ".join(missing)
                + ". Rebuild the WGS cohort, rerun matching, and rerun characterization before GWAS."
            )

    return effective, paths, matched_df


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="AoU workbench CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    preflight_parser = subparsers.add_parser("preflight", help="Run runtime and input validation.")
    _add_config_arguments(preflight_parser)

    build_parser = subparsers.add_parser("build-cohort", help="Build the tiered rhabdomyolysis cohort.")
    _add_config_arguments(build_parser)
    build_parser.add_argument("--require-wgs", action="store_true", help="Restrict the saved cohort to participants with CDR WGS availability.")

    match_parser = subparsers.add_parser("match-controls", help="Build the matched case-control cohort.")
    _add_config_arguments(match_parser)
    match_parser.add_argument("--require-wgs", action="store_true", help="Restrict the cohort and matching universe to participants with CDR WGS availability.")

    eir_build_parser = subparsers.add_parser(
        "build-eir-cohort",
        help="Build the incident CK-confirmed, non-traumatic/non-septic EIR-enriched cohort.",
    )
    _add_eir_config_arguments(eir_build_parser)
    eir_build_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Estimate the BigQuery bytes scanned and optional billing cap without executing the cohort query.",
    )
    eir_build_parser.add_argument(
        "--max-tib",
        type=float,
        default=None,
        help="Optional BigQuery maximum bytes billed cap in TiB. The query fails before billing if the estimate is higher.",
    )
    eir_build_parser.add_argument(
        "--write-sql",
        default=None,
        help="Optional local path to write the rendered BigQuery SQL for review.",
    )

    eir_characterize_parser = subparsers.add_parser(
        "characterize-eir-cohort",
        help="Write EIR-enriched CONSORT counts, train/test splits, Table 1, and missingness reports.",
    )
    _add_eir_config_arguments(eir_characterize_parser)

    eir_model_parser = subparsers.add_parser(
        "run-eir-clinical-model",
        help="Train and evaluate the leakage-controlled clinical-only EIR-enriched incident prediction model.",
    )
    _add_eir_config_arguments(eir_model_parser)
    eir_model_parser.add_argument(
        "--l2-penalty",
        type=float,
        default=1.0,
        help="L2 penalty for regularized logistic regression. Default: 1.0.",
    )
    eir_model_parser.add_argument(
        "--run-sparse",
        action="store_true",
        help="Request the optional sparse OMOP model. The v1 command records a clear skipped status if unavailable.",
    )

    incident_feasibility_parser = subparsers.add_parser(
        "incident-rhabdo-feasibility",
        help="Run aggregate feasibility counts for incident non-traumatic rhabdomyolysis prediction.",
    )
    _add_eir_config_arguments(incident_feasibility_parser)
    incident_feasibility_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Estimate BigQuery bytes scanned and optional billing cap without executing the feasibility query.",
    )
    incident_feasibility_parser.add_argument(
        "--max-tib",
        type=float,
        default=None,
        help="Optional BigQuery maximum bytes billed cap in TiB. The query fails before billing if the estimate is higher.",
    )
    incident_feasibility_parser.add_argument(
        "--write-sql",
        default=None,
        help="Optional local path to write the rendered aggregate feasibility SQL for review.",
    )
    incident_feasibility_parser.add_argument(
        "--microarray-fam",
        default=None,
        help="Optional local AoU microarray PLINK .fam path for case/control microarray overlap counts.",
    )
    incident_feasibility_parser.add_argument(
        "--from-cohort-tsv",
        default=None,
        help="Optional existing EIR built_cohort.tsv to summarize locally instead of running BigQuery.",
    )

    wgs_manifest_parser = subparsers.add_parser(
        "prepare-wgs-manifest",
        help="Optionally write the CDR WGS sample IDs for audit/debugging. WGS-restricted cohort runs filter directly in SQL.",
    )
    _add_config_arguments(wgs_manifest_parser)

    prepare_stage1_parser = subparsers.add_parser(
        "prepare-stage1",
        help="Prepare the Stage 1 exact-variant genotype table from AoU genomics callsets.",
    )
    _add_config_arguments(prepare_stage1_parser)

    prepare_stage2_parser = subparsers.add_parser(
        "prepare-stage2",
        help="Prepare the Stage 2 ClinVar genotype table from AoU smaller callsets.",
    )
    _add_config_arguments(prepare_stage2_parser)

    prepare_regenie_parser = subparsers.add_parser(
        "prepare-regenie",
        help="Write matched phenotype/covariate files and a REGENIE ACAF command template.",
    )
    _add_config_arguments(prepare_regenie_parser)

    prepare_gwas_parser = subparsers.add_parser(
        "prepare-gwas",
        help="Generate a terminal-first matched-control GWAS workspace with REGENIE and dsub templates.",
    )
    _add_config_arguments(prepare_gwas_parser)

    prepare_stage4_parser = subparsers.add_parser(
        "prepare-stage4",
        help="Subset the ACAF smaller callset with Hail, export local PLINK files, and run a Hail pilot GWAS.",
    )
    _add_config_arguments(prepare_stage4_parser)
    prepare_stage4_parser.add_argument("--chromosome", default="chr19")

    hail_gwas_parser = subparsers.add_parser(
        "run-hail-gwas",
        help="Run a Hail-native matched-control GWAS directly against the AoU ACAF split MT without local chromosome rewrites.",
    )
    _add_config_arguments(hail_gwas_parser)
    hail_gwas_parser.add_argument(
        "--chromosomes",
        default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22",
        help="Comma-separated chromosome list, e.g. '1,19,22'.",
    )

    hail_pilot_parser = subparsers.add_parser(
        "run-hail-pilot-gwas",
        help="Run a reduced-marker Hail-native GWAS pilot with train-only sample filtering and strict variant QC.",
    )
    _add_config_arguments(hail_pilot_parser)
    hail_pilot_parser.add_argument(
        "--chromosomes",
        default="22",
        help="Comma-separated autosome list. Default: 22.",
    )
    hail_pilot_parser.add_argument(
        "--genotype-source",
        choices=["acaf", "microarray"],
        default="acaf",
        help="Genotype MatrixTable source for the pilot. Default: acaf.",
    )
    hail_pilot_parser.add_argument(
        "--min-maf",
        type=float,
        default=0.05,
        help="Minimum analysis-sample minor allele frequency. Default: 0.05.",
    )
    hail_pilot_parser.add_argument(
        "--min-mac",
        type=int,
        default=20,
        help="Minimum analysis-sample minor allele count. Default: 20.",
    )
    hail_pilot_parser.add_argument(
        "--min-call-rate",
        type=float,
        default=0.98,
        help="Minimum analysis-sample variant call rate. Default: 0.98.",
    )
    hail_pilot_parser.add_argument(
        "--hwe-p-control",
        type=float,
        default=1e-6,
        help="Minimum Hardy-Weinberg equilibrium p-value among controls. Default: 1e-6.",
    )
    hail_pilot_parser.add_argument(
        "--hwe-filter-mode",
        choices=["filter", "report-only"],
        default="filter",
        help="Use control-HWE as a hard filter or report it without filtering. Default: filter.",
    )
    hail_pilot_parser.add_argument(
        "--analysis-split",
        default="train",
        help="Matched cohort analysis_split to use. Default: train.",
    )
    hail_pilot_parser.add_argument(
        "--eligibility-flag",
        default="primary_model_eligible",
        help="Eligibility flag in matched cohort/clinical model input. Default: primary_model_eligible.",
    )
    hail_pilot_parser.add_argument(
        "--label",
        default=None,
        help="Output label under stage4/hail_pilot/. Defaults to <genotype-source>_chr<chromosomes>_maf<min-maf>_train_qc.",
    )
    hail_pilot_parser.add_argument(
        "--target-partitions",
        type=int,
        default=None,
        help="Optional row partition count after interval/sample/biallelic filtering. Default: source-specific.",
    )
    hail_pilot_parser.add_argument(
        "--write-qc-mt",
        action="store_true",
        help="Write the QC-passing genotype MatrixTable to the workspace bucket and run GWAS from that persisted MT.",
    )
    hail_pilot_parser.add_argument(
        "--export-hail-results-tsv",
        action="store_true",
        help="Also export full Hail GWAS results as sharded TSV in the workspace bucket. The Hail Table is always written.",
    )
    hail_pilot_parser.add_argument(
        "--results-preview-n",
        type=int,
        default=100000,
        help="Number of top p-value results to collect locally for preview plots/tables. Default: 100000.",
    )

    microarray_plink_parser = subparsers.add_parser(
        "run-microarray-plink-gwas",
        help="Run the train-only matched-control GWAS pilot using local AoU v8 microarray PLINK BED/BIM/FAM files.",
    )
    _add_config_arguments(microarray_plink_parser)
    microarray_plink_parser.add_argument(
        "--chromosomes",
        default="22",
        help="Comma-separated autosome list. Default: 22.",
    )
    microarray_plink_parser.add_argument(
        "--plink-prefix",
        default=None,
        help="Local PLINK prefix for arrays.bed/bim/fam. Default: ~/plink_microarray/arrays.",
    )
    microarray_plink_parser.add_argument(
        "--copy-plink-to",
        default=None,
        help="Copy AoU microarray PLINK files into this local directory before running. Existing files are reused.",
    )
    microarray_plink_parser.add_argument(
        "--overwrite-plink-copy",
        action="store_true",
        help="Re-copy AoU microarray PLINK files even when local files already exist.",
    )
    microarray_plink_parser.add_argument(
        "--min-maf",
        type=float,
        default=0.05,
        help="Minimum analysis-sample minor allele frequency. Default: 0.05.",
    )
    microarray_plink_parser.add_argument(
        "--min-mac",
        type=int,
        default=20,
        help="Minimum analysis-sample minor allele count. Default: 20.",
    )
    microarray_plink_parser.add_argument(
        "--min-call-rate",
        type=float,
        default=0.98,
        help="Minimum analysis-sample variant call rate. Default: 0.98.",
    )
    microarray_plink_parser.add_argument(
        "--hwe-p-control",
        type=float,
        default=1e-6,
        help="Minimum Hardy-Weinberg equilibrium p-value among controls. Default: 1e-6.",
    )
    microarray_plink_parser.add_argument(
        "--analysis-split",
        default="train",
        help="Matched cohort analysis_split to use. Default: train.",
    )
    microarray_plink_parser.add_argument(
        "--eligibility-flag",
        default="primary_model_eligible",
        help="Eligibility flag in matched cohort/clinical model input. Default: primary_model_eligible.",
    )
    microarray_plink_parser.add_argument(
        "--label",
        default=None,
        help="Output label under stage4/microarray_plink/. Defaults to microarray_plink_chr<chromosomes>_maf<min-maf>_train_qc.",
    )
    microarray_plink_parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Optional PLINK2 thread count.",
    )
    microarray_plink_parser.add_argument(
        "--memory-mb",
        type=int,
        default=None,
        help="Optional PLINK2 memory limit in MiB.",
    )
    microarray_plink_parser.add_argument(
        "--plink2-bin",
        default="plink2",
        help="PLINK2 executable. Default: plink2.",
    )

    microarray_bigsnpr_parser = subparsers.add_parser(
        "run-microarray-bigsnpr-model",
        help="Prepare a QC'd microarray PLINK subset and run a sparse bigsnpr/bigstatsr SNP prediction model.",
    )
    _add_config_arguments(microarray_bigsnpr_parser)
    microarray_bigsnpr_parser.add_argument(
        "--chromosomes",
        default="22",
        help="Comma-separated autosome list. Default: 22.",
    )
    microarray_bigsnpr_parser.add_argument(
        "--plink-prefix",
        default=None,
        help="Local PLINK prefix for arrays.bed/bim/fam. Default: ~/plink_microarray/arrays.",
    )
    microarray_bigsnpr_parser.add_argument(
        "--copy-plink-to",
        default=None,
        help="Copy AoU microarray PLINK files into this local directory before running. Existing files are reused.",
    )
    microarray_bigsnpr_parser.add_argument(
        "--overwrite-plink-copy",
        action="store_true",
        help="Re-copy AoU microarray PLINK files even when local files already exist.",
    )
    microarray_bigsnpr_parser.add_argument("--min-maf", type=float, default=0.05, help="Minimum train-sample MAF. Default: 0.05.")
    microarray_bigsnpr_parser.add_argument("--min-mac", type=int, default=20, help="Minimum train-sample MAC. Default: 20.")
    microarray_bigsnpr_parser.add_argument(
        "--min-call-rate",
        type=float,
        default=0.98,
        help="Minimum train-sample variant call rate. Default: 0.98.",
    )
    microarray_bigsnpr_parser.add_argument(
        "--hwe-p-control",
        type=float,
        default=1e-6,
        help="Minimum HWE p-value among train controls. Default: 1e-6.",
    )
    microarray_bigsnpr_parser.add_argument(
        "--eligibility-flag",
        default="primary_model_eligible",
        help="Eligibility flag in clinical_model_input.tsv. Default: primary_model_eligible.",
    )
    microarray_bigsnpr_parser.add_argument(
        "--label",
        default=None,
        help="Output label under stage4/microarray_bigsnpr/. Defaults to microarray_bigsnpr_chr<chromosomes>_maf<min-maf>_train_qc.",
    )
    microarray_bigsnpr_parser.add_argument("--threads", type=int, default=None, help="Optional PLINK2/R thread count.")
    microarray_bigsnpr_parser.add_argument("--memory-mb", type=int, default=None, help="Optional PLINK2 memory limit in MiB.")
    microarray_bigsnpr_parser.add_argument("--plink2-bin", default="plink2", help="PLINK2 executable. Default: plink2.")
    microarray_bigsnpr_parser.add_argument("--rscript-bin", default="Rscript", help="Rscript executable. Default: Rscript.")
    microarray_bigsnpr_parser.add_argument(
        "--alphas",
        default=",".join(f"{value:g}" for value in BIGSNPR_DEFAULT_ALPHAS),
        help="Comma-separated elastic-net alphas for big_spLogReg. Default: 1,0.5,0.1,0.01.",
    )
    microarray_bigsnpr_parser.add_argument("--folds", type=int, default=5, help="Internal CMSA fold count. Default: 5.")
    microarray_bigsnpr_parser.add_argument("--nlambda", type=int, default=100, help="Number of lambda values. Default: 100.")
    microarray_bigsnpr_parser.add_argument("--dfmax", type=int, default=50000, help="Maximum nonzero coefficients. Default: 50000.")
    microarray_bigsnpr_parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare PLINK subset and R script but do not run Rscript.",
    )
    microarray_bigsnpr_parser.add_argument(
        "--force-plink-subset",
        action="store_true",
        help="Rebuild the QC'd PLINK subset even if it already exists.",
    )

    microarray_prs_parser = subparsers.add_parser(
        "run-microarray-plink-prs",
        help="Build LD-clumped threshold-grid PRS scores from a microarray PLINK GWAS and score the held-out split.",
    )
    _add_config_arguments(microarray_prs_parser)
    microarray_prs_parser.add_argument(
        "--gwas-label",
        required=True,
        help="Microarray PLINK GWAS label to use as the PRS weight source.",
    )
    microarray_prs_parser.add_argument(
        "--plink-prefix",
        default=None,
        help="Local PLINK prefix for arrays.bed/bim/fam. Default: ~/plink_microarray/arrays.",
    )
    microarray_prs_parser.add_argument(
        "--plink2-bin",
        default="plink2",
        help="PLINK2 executable. Default: plink2.",
    )
    microarray_prs_parser.add_argument(
        "--score-split",
        default="test",
        help="Matched cohort analysis_split to score. Default: test.",
    )
    microarray_prs_parser.add_argument(
        "--eligibility-flag",
        default="primary_model_eligible",
        help="Eligibility flag in clinical_model_input.tsv. Default: primary_model_eligible.",
    )
    microarray_prs_parser.add_argument("--clump-r2", type=float, default=0.1, help="PLINK clump r2 threshold. Default: 0.1.")
    microarray_prs_parser.add_argument("--clump-kb", type=int, default=250, help="PLINK clump window in kb. Default: 250.")
    microarray_prs_parser.add_argument("--clump-p1", type=float, default=0.01, help="PLINK clump index p threshold. Default: 0.01.")
    microarray_prs_parser.add_argument("--clump-p2", type=float, default=0.01, help="PLINK clump secondary p threshold. Default: 0.01.")
    microarray_prs_parser.add_argument(
        "--thresholds",
        default=",".join(f"{value:g}" for value in DEFAULT_PRS_THRESHOLDS),
        help="Comma-separated PRS p-value thresholds.",
    )
    microarray_prs_parser.add_argument(
        "--label",
        default=None,
        help="Output label under stage4/microarray_plink/<gwas-label>/prs/. Default: <score-split>_clumped_threshold_grid.",
    )
    microarray_prs_parser.add_argument("--threads", type=int, default=None, help="Optional PLINK2 thread count.")
    microarray_prs_parser.add_argument("--memory-mb", type=int, default=None, help="Optional PLINK2 memory limit in MiB.")

    for name in ("run-stage1", "run-stage2", "run-stage3", "run-stage4", "run-all"):
        stage_parser = subparsers.add_parser(name, help=f"Execute {name}.")
        _add_config_arguments(stage_parser)
        stage_parser.add_argument("--skip-preflight", action="store_true")

    report_parser = subparsers.add_parser("report", help="Rebuild the final markdown report from existing outputs.")
    _add_config_arguments(report_parser)

    presentation_parser = subparsers.add_parser(
        "presentation-dashboard",
        help="Write a supervisor-facing static HTML dashboard from existing aggregate outputs.",
    )
    _add_config_arguments(presentation_parser)
    presentation_parser.add_argument(
        "--gwas-label",
        default="microarray_plink_autosomes_maf05_train_qc",
        help="Microarray PLINK GWAS label to summarize. Default: microarray_plink_autosomes_maf05_train_qc.",
    )
    presentation_parser.add_argument(
        "--prs-label",
        default="test-clumped-p001",
        help="PRS label under the selected microarray GWAS. Default: test-clumped-p001.",
    )
    presentation_parser.add_argument(
        "--clinical-prs-label",
        default="clinical_prs_p001",
        help="Clinical+PRS model label. Default: clinical_prs_p001.",
    )
    presentation_parser.add_argument(
        "--diagnostics-label",
        default="prs_p001_diagnostics",
        help="PRS diagnostics label used as a source note. Default: prs_p001_diagnostics.",
    )
    presentation_parser.add_argument(
        "--output",
        default=None,
        help="Optional dashboard HTML output path. Default: outputs/<run>/rhabdo_presentation_dashboard.html.",
    )

    summary_parser = subparsers.add_parser(
        "summarize-cohort",
        help="Write a clinical/demographic comparison table for the WGS-restricted unmatched and matched analyses.",
    )
    _add_config_arguments(summary_parser)

    characterize_parser = subparsers.add_parser(
        "characterize-cohort",
        help="Write CONSORT counts, matched Table 1, and sepsis/renal injury timing summaries.",
    )
    _add_config_arguments(characterize_parser)
    characterize_parser.add_argument("--require-wgs", action="store_true", help="Build/load WGS-restricted cohort artifacts.")

    clinical_model_parser = subparsers.add_parser(
        "run-clinical-model",
        help="Train and evaluate the clinical-only rhabdomyolysis prediction model.",
    )
    _add_config_arguments(clinical_model_parser)
    clinical_model_parser.add_argument(
        "--eligibility-flag",
        default="primary_model_eligible",
        help="Eligibility flag in clinical_model_input.tsv. Default: primary_model_eligible.",
    )
    clinical_model_parser.add_argument(
        "--l2-penalty",
        type=float,
        default=1.0,
        help="L2 penalty for regularized logistic regression. Default: 1.0.",
    )
    clinical_model_parser.add_argument("--require-wgs", action="store_true", help="Build/load WGS-restricted cohort artifacts.")

    model_comparison_parser = subparsers.add_parser(
        "compare-prs-models",
        help="Compare held-out clinical-only and PRS-only model performance on participants with both predictions.",
    )
    _add_config_arguments(model_comparison_parser)
    model_comparison_parser.add_argument(
        "--gwas-label",
        required=True,
        help="Microarray PLINK GWAS label used as the PRS source.",
    )
    model_comparison_parser.add_argument(
        "--prs-label",
        required=True,
        help="PRS label under stage4/microarray_plink/<gwas-label>/prs/.",
    )
    model_comparison_parser.add_argument(
        "--label",
        default=None,
        help="Output label under clinical/model_comparison/. Default: <gwas-label>_<prs-label>_heldout.",
    )
    model_comparison_parser.add_argument(
        "--eligibility-flag",
        default="primary_model_eligible",
        help="Eligibility flag for rerunning the clinical model if needed. Default: primary_model_eligible.",
    )
    model_comparison_parser.add_argument(
        "--l2-penalty",
        type=float,
        default=1.0,
        help="L2 penalty if the clinical model must be rerun. Default: 1.0.",
    )
    model_comparison_parser.add_argument(
        "--rerun-clinical-model",
        action="store_true",
        help="Force rerunning the clinical model before comparing.",
    )

    clinical_prs_parser = subparsers.add_parser(
        "run-clinical-prs-model",
        help="Train a pragmatic clinical+PRS model on train and evaluate once on held-out test.",
    )
    _add_config_arguments(clinical_prs_parser)
    clinical_prs_parser.add_argument(
        "--gwas-label",
        required=True,
        help="Microarray PLINK GWAS label used as the PRS source.",
    )
    clinical_prs_parser.add_argument(
        "--prs-label",
        required=True,
        help="Existing test PRS label under stage4/microarray_plink/<gwas-label>/prs/.",
    )
    clinical_prs_parser.add_argument(
        "--plink-prefix",
        default=None,
        help="Local PLINK prefix for arrays.bed/bim/fam. Default: ~/plink_microarray/arrays.",
    )
    clinical_prs_parser.add_argument(
        "--plink2-bin",
        default="plink2",
        help="PLINK2 executable. Default: plink2.",
    )
    clinical_prs_parser.add_argument(
        "--threshold-label",
        default="p0_01",
        help="Threshold label to use from prs_scores.tsv. Default: p0_01.",
    )
    clinical_prs_parser.add_argument(
        "--p-threshold",
        type=float,
        default=0.01,
        help="P-value threshold used to filter clumped PRS weights for train scoring. Default: 0.01.",
    )
    clinical_prs_parser.add_argument(
        "--label",
        default=None,
        help="Output label under clinical/clinical_prs_model/. Default: clinical_prs_<prs-label>.",
    )
    clinical_prs_parser.add_argument(
        "--eligibility-flag",
        default="primary_model_eligible",
        help="Eligibility flag in clinical_model_input.tsv. Default: primary_model_eligible.",
    )
    clinical_prs_parser.add_argument(
        "--l2-penalty",
        type=float,
        default=1.0,
        help="L2 penalty for regularized logistic regression. Default: 1.0.",
    )
    clinical_prs_parser.add_argument("--threads", type=int, default=None, help="Optional PLINK2 thread count for train scoring.")
    clinical_prs_parser.add_argument("--memory-mb", type=int, default=None, help="Optional PLINK2 memory limit in MiB.")

    prs_diagnostics_parser = subparsers.add_parser(
        "diagnose-prs",
        help="Generate PRS diagnostics from existing PRS and clinical+PRS outputs without rerunning GWAS or scoring.",
    )
    _add_config_arguments(prs_diagnostics_parser)
    prs_diagnostics_parser.add_argument(
        "--gwas-label",
        required=True,
        help="Microarray PLINK GWAS label used as the PRS source.",
    )
    prs_diagnostics_parser.add_argument(
        "--prs-label",
        default="test-clumped-p001",
        help="Existing PRS label under stage4/microarray_plink/<gwas-label>/prs/. Default: test-clumped-p001.",
    )
    prs_diagnostics_parser.add_argument(
        "--clinical-prs-label",
        default="clinical_prs_p001",
        help="Existing clinical+PRS model label under clinical/clinical_prs_model/. Use empty string to skip clinical+PRS diagnostics.",
    )
    prs_diagnostics_parser.add_argument(
        "--threshold-label",
        default="p0_01",
        help="PRS threshold label to diagnose. Default: p0_01.",
    )
    prs_diagnostics_parser.add_argument(
        "--bootstrap-iterations",
        type=int,
        default=200,
        help="Bootstrap iterations for test-set confidence intervals. Default: 200.",
    )
    prs_diagnostics_parser.add_argument(
        "--seed",
        type=int,
        default=20260611,
        help="Random seed for bootstrap intervals. Default: 20260611.",
    )
    prs_diagnostics_parser.add_argument(
        "--label",
        default="prs_p001_diagnostics",
        help="Output label under clinical/prs_diagnostics/. Default: prs_p001_diagnostics.",
    )

    preindex_parser = subparsers.add_parser(
        "profile-preindex-cases",
        help="Summarize what condition, lab, and biomarker data are available before rhabdo index dates for cases.",
    )
    _add_config_arguments(preindex_parser)
    preindex_parser.add_argument(
        "--case-tier",
        default=None,
        help="Case tier to profile. Defaults to cohort.primary_case_tier.",
    )
    preindex_parser.add_argument(
        "--windows",
        default="365,1095,all",
        help="Comma-separated pre-index windows in days, plus 'all'. Default: 365,1095,all.",
    )
    preindex_parser.add_argument("--top-n", type=int, default=25, help="Top concepts to keep per window.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    config = _load_config(args)

    if args.command == "preflight":
        print(format_preflight_report(run_preflight_checks(config)))
        return 0

    if args.command == "build-cohort":
        _, paths, cohort_df = build_cohort_artifacts(config, require_wgs=args.require_wgs)
        print(f"Built cohort rows: {len(cohort_df)}")
        print(f"Cohort path: {paths.built_cohort_tsv}")
        return 0

    if args.command == "match-controls":
        effective, paths, cohort_df = _load_or_build_cohort_artifacts(config, require_wgs=args.require_wgs)
        _, paths, matched_df = match_controls_artifacts(effective, cohort_df, require_wgs=args.require_wgs)
        print(f"Matched cohort rows: {len(matched_df)}")
        print(f"Matched cohort path: {paths.matched_cohort_tsv}")
        return 0

    if args.command == "build-eir-cohort":
        if args.dry_run:
            _, _, estimate = estimate_eir_cohort_artifacts(
                config,
                max_tib=args.max_tib,
                write_sql_path=args.write_sql,
            )
            total_bytes = int(estimate.get("total_bytes_processed", 0))
            total_tib = float(estimate.get("total_tib_processed", 0.0))
            estimated_cost = float(estimate.get("estimated_query_cost_usd", 0.0))
            max_bytes = estimate.get("maximum_bytes_billed")
            print("EIR cohort BigQuery dry run")
            print(f"Mode: {estimate.get('mode', 'unknown')}")
            print(f"Estimated bytes processed: {total_bytes}")
            print(f"Estimated TiB processed: {total_tib:.4f}")
            print(f"Approx on-demand query cost: ${estimated_cost:.2f}")
            if max_bytes is not None:
                print(f"Maximum bytes billed cap: {int(max_bytes)} ({int(max_bytes) / float(1024**4):.4f} TiB)")
                if estimate.get("would_exceed_maximum_bytes_billed"):
                    print("Cap status: estimated bytes exceed this cap; the real query would be refused before billing.")
                else:
                    print("Cap status: estimated bytes are within this cap.")
            if estimate.get("sql_path"):
                print(f"SQL path: {estimate['sql_path']}")
            if estimate.get("message"):
                print(str(estimate["message"]))
            return 0
        _, paths, cohort_df = build_eir_cohort_artifacts(
            config,
            max_tib=args.max_tib,
            write_sql_path=args.write_sql,
        )
        print(f"EIR cohort rows: {len(cohort_df)}")
        print(f"Primary EIR-enriched cases: {int(cohort_df['eir_primary_case'].sum())}")
        print(f"Eligible controls: {int(cohort_df['eligible_control'].sum())}")
        print(f"Cohort path: {paths.built_cohort_tsv}")
        return 0

    if args.command == "characterize-eir-cohort":
        _, paths, outputs = characterize_eir_artifacts(config)
        print(f"EIR CONSORT rows: {outputs['consort'].shape[0]}")
        print(f"EIR model input rows: {outputs['model_input'].shape[0]}")
        print(f"EIR Table 1 rows: {outputs['table1'].shape[0]}")
        print(f"EIR missingness rows: {outputs['missingness'].shape[0]}")
        print(f"consort: {eir_consort_counts_path(paths)}")
        print(f"split_summary: {eir_split_summary_path(paths)}")
        print(f"table1: {eir_table1_path(paths)}")
        print(f"missingness: {eir_missingness_path(paths)}")
        print(f"clinical_model_input: {eir_model_input_path(paths)}")
        print(f"report: {eir_characterization_report_path(paths)}")
        return 0

    if args.command == "run-eir-clinical-model":
        _, paths, _ = characterize_eir_artifacts(config)
        outputs = run_eir_clinical_model(
            config,
            paths,
            l2_penalty=args.l2_penalty,
            run_sparse=args.run_sparse,
        )
        metrics = outputs["metrics"]
        test_metrics = metrics[metrics["evaluation_set"] == "test"].iloc[0]
        print(f"EIR clinical model test ROC AUC: {test_metrics['roc_auc']:.4f}")
        print(f"EIR clinical model test average precision: {test_metrics['average_precision']:.4f}")
        print(f"metrics: {clinical_model_metrics_path(paths)}")
        print(f"cv_metrics: {clinical_model_cv_metrics_path(paths)}")
        print(f"coefficients: {clinical_model_coefficients_path(paths)}")
        print(f"predictions: {clinical_model_predictions_path(paths)}")
        print(f"calibration: {clinical_model_calibration_path(paths)}")
        print(f"report: {clinical_model_report_path(paths)}")
        print(f"roc: {clinical_model_roc_svg_path(paths)}")
        print(f"precision_recall: {clinical_model_pr_svg_path(paths)}")
        print(f"calibration_plot: {clinical_model_calibration_svg_path(paths)}")
        print(f"risk_deciles: {eir_risk_decile_path(paths)}")
        print(f"sparse_status: {eir_sparse_status_path(paths)}")
        print(f"sparse_report: {eir_sparse_report_path(paths)}")
        return 0

    if args.command == "incident-rhabdo-feasibility":
        if args.dry_run:
            _, _, estimate = estimate_incident_feasibility_artifacts(
                config,
                max_tib=args.max_tib,
                write_sql_path=args.write_sql,
            )
            total_bytes = int(estimate.get("total_bytes_processed", 0))
            total_tib = float(estimate.get("total_tib_processed", 0.0))
            estimated_cost = float(estimate.get("estimated_query_cost_usd", 0.0))
            max_bytes = estimate.get("maximum_bytes_billed")
            print("Incident rhabdo feasibility BigQuery dry run")
            print(f"Mode: {estimate.get('mode', 'unknown')}")
            print(f"Estimated bytes processed: {total_bytes}")
            print(f"Estimated TiB processed: {total_tib:.4f}")
            print(f"Approx on-demand query cost: ${estimated_cost:.2f}")
            if max_bytes is not None:
                print(f"Maximum bytes billed cap: {int(max_bytes)} ({int(max_bytes) / float(1024**4):.4f} TiB)")
                if estimate.get("would_exceed_maximum_bytes_billed"):
                    print("Cap status: estimated bytes exceed this cap; the real query would be refused before billing.")
                else:
                    print("Cap status: estimated bytes are within this cap.")
            if estimate.get("sql_path"):
                print(f"SQL path: {estimate['sql_path']}")
            if estimate.get("message"):
                print(str(estimate["message"]))
            if args.microarray_fam:
                print("Microarray FAM overlap is skipped during dry-run.")
            if args.from_cohort_tsv:
                print("--from-cohort-tsv is ignored during dry-run.")
            return 0
        _, paths, outputs = run_incident_feasibility(
            config,
            max_tib=args.max_tib,
            write_sql_path=args.write_sql,
            microarray_fam=args.microarray_fam,
            from_cohort_tsv=args.from_cohort_tsv,
        )
        print(f"Incident rhabdo feasibility rows: {outputs['feasibility_counts'].shape[0]}")
        print(f"feasibility_counts: {incident_feasibility_counts_path(paths)}")
        print(f"case_funnel: {incident_case_funnel_path(paths)}")
        print(f"control_funnel: {incident_control_funnel_path(paths)}")
        print(f"baseline_history_bins: {incident_baseline_history_bins_path(paths)}")
        if "microarray_overlap_counts" in outputs:
            print(f"microarray_overlap_counts: {incident_microarray_overlap_counts_path(paths)}")
        print(f"report: {incident_feasibility_report_path(paths)}")
        return 0

    if args.command == "prepare-wgs-manifest":
        manifest = prepare_wgs_sample_manifest(config)
        effective = apply_runtime_defaults(config)
        stage = effective.analysis.stage1
        print(f"WGS/ACAF manifest rows: {len(manifest)}")
        if stage is not None:
            from .stage1_prepare import stage1_sample_manifest_path

            print(stage1_sample_manifest_path(stage.variant_table))
        return 0

    if args.command == "prepare-stage1":
        effective, _, cohort_df = _load_or_build_cohort_artifacts(config)
        frame = prepare_stage1_variant_table(effective, cohort_df)
        stage = effective.analysis.stage1
        if stage is None:
            print("Stage 1 is not configured.")
            return 0
        print(f"Prepared Stage 1 rows: {frame.shape[0]}")
        print(stage.variant_table)
        return 0

    if args.command == "prepare-stage2":
        effective, _, cohort_df = _load_or_build_cohort_artifacts(config)
        frame = prepare_stage2_variant_table(effective, cohort_df)
        stage = effective.analysis.stage2
        if stage is None:
            print("Stage 2 is not configured.")
            return 0
        print(f"Prepared Stage 2 rows: {frame.shape[0]}")
        print(stage.variant_table)
        return 0

    if args.command == "prepare-regenie":
        effective, paths, matched_df = _load_or_build_matched_artifacts(config)
        outputs = prepare_regenie_inputs(effective, matched_df, paths)
        keep_df = read_table(outputs["keep"])
        print(
            f"Prepared REGENIE inputs for {len(keep_df)} GWAS-eligible matched samples "
            f"(from {matched_df['person_id'].astype(str).nunique()} matched participants)."
        )
        for name, path in outputs.items():
            print(f"{name}: {path}")
        return 0

    if args.command == "prepare-gwas":
        effective, paths, matched_df = _load_or_build_matched_artifacts(config)
        outputs = prepare_terminal_gwas_workspace(effective, matched_df, paths)
        restricted = read_table(outputs["restricted_manifest"])
        print(
            f"Prepared GWAS workspace for {len(restricted)} GWAS-eligible matched samples "
            f"(from {matched_df['person_id'].astype(str).nunique()} matched participants)."
        )
        for name, path in outputs.items():
            print(f"{name}: {path}")
        return 0

    if args.command == "prepare-stage4":
        effective, paths, matched_df = _load_or_build_matched_artifacts(config)
        outputs = prepare_stage4_acaf_subset(effective, matched_df, paths, chromosome=args.chromosome)
        print(
            f"Prepared Stage 4 ACAF subset on {args.chromosome} for {matched_df['person_id'].astype(str).nunique()} matched samples."
        )
        for name, path in outputs.items():
            print(f"{name}: {path}")
        return 0

    if args.command == "run-hail-gwas":
        effective, paths, matched_df = _load_or_build_matched_artifacts(config)
        chromosomes = [value.strip() for value in args.chromosomes.split(",") if value.strip()]
        full, hits = run_stage4_hail_gwas(effective, matched_df, paths, chromosomes=chromosomes)
        print(f"Hail GWAS variants tested: {full.shape[0]}")
        print(f"Hail GWAS lead hits: {hits.shape[0]}")
        print(f"full_results: {hail_stage4_full_results_path(paths)}")
        print(f"lead_hits: {hail_stage4_lead_hits_path(paths)}")
        print(f"qc: {hail_stage4_qc_path(paths)}")
        print(f"report: {hail_stage4_report_path(paths)}")
        print(f"manhattan: {hail_stage4_manhattan_path(paths)}")
        print(f"qq: {hail_stage4_qq_path(paths)}")
        return 0

    if args.command == "run-hail-pilot-gwas":
        effective, paths, matched_df = _load_hail_pilot_matched_input(
            config,
            eligibility_flag=args.eligibility_flag,
        )
        chromosomes = [value.strip() for value in args.chromosomes.split(",") if value.strip()]
        label = args.label or hail_pilot_default_label(args.genotype_source, chromosomes, min_maf=args.min_maf)
        full, hits = run_stage4_hail_pilot_gwas(
            effective,
            matched_df,
            paths,
            chromosomes=chromosomes,
            min_maf=args.min_maf,
            min_mac=args.min_mac,
            min_call_rate=args.min_call_rate,
            hwe_p_control=args.hwe_p_control,
            hwe_filter_mode=args.hwe_filter_mode,
            analysis_split=args.analysis_split,
            eligibility_flag=args.eligibility_flag,
            label=label,
            genotype_source=args.genotype_source,
            target_partitions=args.target_partitions,
            write_qc_mt=args.write_qc_mt,
            export_hail_results_tsv=args.export_hail_results_tsv,
            results_preview_n=args.results_preview_n,
        )
        print(f"Hail pilot GWAS preview rows: {full.shape[0]}")
        print(f"Hail pilot GWAS lead hits: {hits.shape[0]}")
        print(f"results_preview: {hail_pilot_results_preview_path(paths, label)}")
        print(f"results_legacy_full_if_preview_complete: {hail_pilot_results_path(paths, label)}")
        print(f"hail_results_ht: {hail_pilot_results_ht_uri(effective, label)}")
        if args.export_hail_results_tsv:
            print(f"hail_results_tsv: {hail_pilot_results_tsv_uri(effective, label)}")
        if args.write_qc_mt:
            print(f"qc_pass_mt: {hail_pilot_qc_pass_mt_uri(effective, label)}")
        print(f"lead_hits: {hail_pilot_lead_hits_path(paths, label)}")
        print(f"qc: {hail_pilot_qc_path(paths, label)}")
        print(f"variant_qc_summary: {hail_pilot_variant_qc_summary_path(paths, label)}")
        print(f"report: {hail_pilot_report_path(paths, label)}")
        print(f"manhattan: {hail_pilot_manhattan_path(paths, label)}")
        print(f"qq: {hail_pilot_qq_path(paths, label)}")
        return 0

    if args.command == "run-microarray-plink-gwas":
        effective, paths, matched_df = _load_hail_pilot_matched_input(
            config,
            eligibility_flag=args.eligibility_flag,
        )
        chromosomes = [value.strip() for value in args.chromosomes.split(",") if value.strip()]
        label = args.label or microarray_plink_default_label(chromosomes, min_maf=args.min_maf)
        full, hits = run_microarray_plink_gwas(
            effective,
            matched_df,
            paths,
            chromosomes=chromosomes,
            min_maf=args.min_maf,
            min_mac=args.min_mac,
            min_call_rate=args.min_call_rate,
            hwe_p_control=args.hwe_p_control,
            analysis_split=args.analysis_split,
            eligibility_flag=args.eligibility_flag,
            label=label,
            plink_prefix=args.plink_prefix,
            copy_plink_to=args.copy_plink_to,
            overwrite_plink_copy=args.overwrite_plink_copy,
            threads=args.threads,
            memory_mb=args.memory_mb,
            plink2_bin=args.plink2_bin,
        )
        print(f"Microarray PLINK GWAS variants tested: {full.shape[0]}")
        print(f"Microarray PLINK GWAS lead hits: {hits.shape[0]}")
        print(f"results: {microarray_plink_results_path(paths, label)}")
        print(f"lead_hits: {microarray_plink_lead_hits_path(paths, label)}")
        print(f"qc: {microarray_plink_qc_path(paths, label)}")
        print(f"variant_qc_summary: {microarray_plink_variant_qc_summary_path(paths, label)}")
        print(f"report: {microarray_plink_report_path(paths, label)}")
        print(f"manhattan: {microarray_plink_manhattan_path(paths, label)}")
        print(f"qq: {microarray_plink_qq_path(paths, label)}")
        return 0

    if args.command == "run-microarray-bigsnpr-model":
        effective, paths, matched_df = _load_hail_pilot_matched_input(
            config,
            eligibility_flag=args.eligibility_flag,
        )
        chromosomes = [value.strip() for value in args.chromosomes.split(",") if value.strip()]
        label = args.label or microarray_bigsnpr_default_label(chromosomes, min_maf=args.min_maf)
        outputs = run_microarray_bigsnpr_model(
            effective,
            matched_df,
            paths,
            chromosomes=chromosomes,
            min_maf=args.min_maf,
            min_mac=args.min_mac,
            min_call_rate=args.min_call_rate,
            hwe_p_control=args.hwe_p_control,
            eligibility_flag=args.eligibility_flag,
            label=label,
            plink_prefix=args.plink_prefix,
            copy_plink_to=args.copy_plink_to,
            overwrite_plink_copy=args.overwrite_plink_copy,
            threads=args.threads,
            memory_mb=args.memory_mb,
            plink2_bin=args.plink2_bin,
            rscript_bin=args.rscript_bin,
            alphas=tuple(parse_thresholds(args.alphas)),
            folds=args.folds,
            nlambda=args.nlambda,
            dfmax=args.dfmax,
            prepare_only=args.prepare_only,
            reuse_plink_subset=not args.force_plink_subset,
        )
        variant_qc = outputs["variant_qc"]
        print(f"Microarray bigsnpr variant QC rows: {variant_qc.shape[0]}")
        print(f"sample_metadata: {microarray_bigsnpr_metadata_path(paths, label)}")
        print(f"variant_qc_summary: {microarray_bigsnpr_variant_qc_summary_path(paths, label)}")
        print(f"qc: {microarray_bigsnpr_qc_path(paths, label)}")
        print(f"r_script: {microarray_bigsnpr_script_path(paths, label)}")
        print(f"report: {microarray_bigsnpr_report_path(paths, label)}")
        if args.prepare_only:
            print("Prepared PLINK subset and R script only; rerun without --prepare-only to fit the bigsnpr model.")
        else:
            print(f"metrics: {microarray_bigsnpr_metrics_path(paths, label)}")
            print(f"predictions: {microarray_bigsnpr_predictions_path(paths, label)}")
            print(f"selected_variants: {microarray_bigsnpr_selected_variants_path(paths, label)}")
            print(f"score_distribution: {microarray_bigsnpr_score_distribution_path(paths, label)}")
            print(f"risk_decile: {microarray_bigsnpr_risk_decile_path(paths, label)}")
        return 0

    if args.command == "run-microarray-plink-prs":
        effective, paths, matched_df = _load_hail_pilot_matched_input(
            config,
            eligibility_flag=args.eligibility_flag,
        )
        label = args.label or microarray_prs_default_label(args.score_split)
        outputs = run_microarray_plink_prs(
            effective,
            matched_df,
            paths,
            gwas_label=args.gwas_label,
            plink_prefix=args.plink_prefix,
            plink2_bin=args.plink2_bin,
            score_split=args.score_split,
            eligibility_flag=args.eligibility_flag,
            clump_r2=args.clump_r2,
            clump_kb=args.clump_kb,
            clump_p1=args.clump_p1,
            clump_p2=args.clump_p2,
            thresholds=parse_thresholds(args.thresholds),
            label=label,
            threads=args.threads,
            memory_mb=args.memory_mb,
        )
        metrics = outputs["metrics"]
        print(f"Microarray PLINK PRS thresholds evaluated: {metrics.shape[0]}")
        print(f"weights: {microarray_prs_weights_path(paths, args.gwas_label, label)}")
        print(f"ranges: {microarray_prs_range_path(paths, args.gwas_label, label)}")
        print(f"scores: {microarray_prs_scores_path(paths, args.gwas_label, label)}")
        print(f"metrics: {microarray_prs_metrics_path(paths, args.gwas_label, label)}")
        print(f"qc: {microarray_prs_qc_path(paths, args.gwas_label, label)}")
        print(f"report: {microarray_prs_report_path(paths, args.gwas_label, label)}")
        print(f"case_status_plot: {microarray_prs_case_status_svg_path(paths, args.gwas_label, label)}")
        return 0

    if args.command == "run-all":
        paths = run_all(config, skip_preflight=args.skip_preflight)
        print(f"Run root: {paths.run_root}")
        print(f"Final report: {paths.final_report_md}")
        print(f"Dashboard: {paths.final_dashboard_html}")
        return 0

    if args.command in {"run-stage1", "run-stage2", "run-stage3", "run-stage4"}:
        if args.command == "run-stage1":
            effective, paths, cohort_df = _load_or_build_cohort_artifacts(config)
            frame = run_stage1_prior_variants(effective, cohort_df, paths)
            print(f"Stage 1 rows: {frame.shape[0]}")
            print(paths.stage1_results_tsv)
            return 0
        if args.command == "run-stage2":
            effective, paths, cohort_df = _load_or_build_cohort_artifacts(config)
            _, gene_df, person_df = run_stage2_plp_panel(effective, cohort_df, paths)
            print(f"Stage 2 genes: {gene_df.shape[0]}")
            print(f"Stage 2 people with hits: {person_df.shape[0]}")
            print(paths.stage2_gene_tsv)
            return 0
        effective, paths, matched_df = _load_or_build_matched_artifacts(config)
        if args.command == "run-stage3":
            frame = run_stage3_burden(effective, matched_df, paths)
            print(f"Stage 3 rows: {frame.shape[0]}")
            print(paths.stage3_results_tsv)
            return 0
        frame, hits = run_stage4_gwas(effective, matched_df, paths)
        print(f"Stage 4 variants tested: {frame.shape[0]}")
        print(f"Stage 4 lead hits: {hits.shape[0]}")
        print(paths.stage4_full_results_tsv)
        return 0

    if args.command == "report":
        paths = render_existing_report(config)
        print(f"Final report: {paths.final_report_md}")
        print(f"Dashboard: {paths.final_dashboard_html}")
        return 0

    if args.command == "presentation-dashboard":
        effective = apply_runtime_defaults(config)
        paths = build_output_paths(effective)
        output = render_presentation_dashboard(
            effective,
            paths,
            gwas_label=args.gwas_label,
            prs_label=args.prs_label,
            clinical_prs_label=args.clinical_prs_label,
            diagnostics_label=args.diagnostics_label,
            output_path=args.output,
        )
        print(f"Presentation dashboard: {output}")
        return 0

    if args.command == "summarize-cohort":
        effective, paths, cohort_df = _load_or_build_cohort_artifacts(config)
        _, _, matched_df = _load_or_build_matched_artifacts(config)
        frame = summarize_clinical_demographics(effective, cohort_df, matched_df, paths)
        print(f"Cohort summary rows: {frame.shape[0]}")
        print(cohort_summary_table_path(paths))
        print(cohort_summary_report_path(paths))
        return 0

    if args.command == "characterize-cohort":
        effective, paths, cohort_df = _load_or_build_cohort_artifacts(config, require_wgs=args.require_wgs)
        _, _, matched_df = _load_or_build_matched_artifacts(config, require_wgs=args.require_wgs)
        outputs = characterize_case_control_cohort(effective, cohort_df, matched_df, paths)
        print(f"CONSORT rows: {outputs['consort'].shape[0]}")
        print(f"Table 1 rows: {outputs['table1'].shape[0]}")
        print(f"Split Table 1 rows: {outputs['split_table1'].shape[0]}")
        print(f"Critical illness rows: {outputs['critical_illness'].shape[0]}")
        print(f"Case cofactor prior timing rows: {outputs['case_cofactor_prior_timing'].shape[0]}")
        print(f"Split summary rows: {outputs['split_summary'].shape[0]}")
        print(f"Eligibility rows: {outputs['eligibility'].shape[0]}")
        print(f"Missingness rows: {outputs['missingness'].shape[0]}")
        print(f"consort: {consort_counts_path(paths)}")
        print(f"consort_report: {consort_counts_report_path(paths)}")
        print(f"table1: {matched_table1_path(paths)}")
        print(f"table1_report: {matched_table1_report_path(paths)}")
        print(f"split_table1: {split_table1_path(paths)}")
        print(f"split_table1_report: {split_table1_report_path(paths)}")
        print(f"critical_illness: {critical_illness_summary_path(paths)}")
        print(f"critical_illness_report: {critical_illness_summary_report_path(paths)}")
        print(f"case_cofactor_prior_timing: {case_cofactor_prior_timing_path(paths)}")
        print(f"case_cofactor_prior_timing_report: {case_cofactor_prior_timing_report_path(paths)}")
        print(f"case_cofactor_prior_timing_histogram: {case_cofactor_prior_timing_histogram_path(paths)}")
        print(f"split_summary: {model_split_summary_path(paths)}")
        print(f"split_summary_report: {model_split_summary_report_path(paths)}")
        print(f"eligibility: {model_eligibility_summary_path(paths)}")
        print(f"eligibility_report: {model_eligibility_summary_report_path(paths)}")
        print(f"missingness: {missingness_summary_path(paths)}")
        print(f"missingness_report: {missingness_summary_report_path(paths)}")
        print(f"clinical_model_input: {clinical_model_input_path(paths)}")
        print(f"report: {clinical_characterization_report_path(paths)}")
        report_paths = render_existing_report(config)
        print(f"analysis_report: {report_paths.final_report_md}")
        print(f"analysis_dashboard: {report_paths.final_dashboard_html}")
        return 0

    if args.command == "run-clinical-model":
        effective, paths, cohort_df = _load_or_build_cohort_artifacts(config, require_wgs=args.require_wgs)
        if _clinical_model_input_needs_refresh(paths, require_wgs=args.require_wgs):
            _, _, matched_df = _load_or_build_matched_artifacts(config, require_wgs=args.require_wgs)
            characterize_case_control_cohort(effective, cohort_df, matched_df, paths)
        outputs = run_clinical_model(
            effective,
            paths,
            eligibility_flag=args.eligibility_flag,
            l2_penalty=args.l2_penalty,
        )
        metrics = outputs["metrics"]
        test_metrics = metrics[metrics["evaluation_set"] == "test"].iloc[0]
        print(f"Clinical model test ROC AUC: {test_metrics['roc_auc']:.4f}")
        print(f"Clinical model test average precision: {test_metrics['average_precision']:.4f}")
        print(f"metrics: {clinical_model_metrics_path(paths)}")
        print(f"cv_metrics: {clinical_model_cv_metrics_path(paths)}")
        print(f"coefficients: {clinical_model_coefficients_path(paths)}")
        print(f"predictions: {clinical_model_predictions_path(paths)}")
        print(f"calibration: {clinical_model_calibration_path(paths)}")
        print(f"report: {clinical_model_report_path(paths)}")
        print(f"roc: {clinical_model_roc_svg_path(paths)}")
        print(f"precision_recall: {clinical_model_pr_svg_path(paths)}")
        print(f"calibration_plot: {clinical_model_calibration_svg_path(paths)}")
        report_paths = render_existing_report(config)
        print(f"analysis_report: {report_paths.final_report_md}")
        print(f"analysis_dashboard: {report_paths.final_dashboard_html}")
        return 0

    if args.command == "compare-prs-models":
        effective = apply_runtime_defaults(config)
        paths = build_output_paths(effective)
        outputs = run_heldout_model_comparison(
            effective,
            paths,
            gwas_label=args.gwas_label,
            prs_label=args.prs_label,
            label=args.label,
            rerun_clinical_model=args.rerun_clinical_model,
            eligibility_flag=args.eligibility_flag,
            l2_penalty=args.l2_penalty,
        )
        metrics = outputs["metrics"]
        print(f"Held-out model comparison rows: {metrics.shape[0]}")
        print(f"metrics: {model_comparison_metrics_path(paths, args.label or f'{args.gwas_label}_{args.prs_label}_heldout')}")
        print(f"predictions: {model_comparison_predictions_path(paths, args.label or f'{args.gwas_label}_{args.prs_label}_heldout')}")
        print(f"report: {model_comparison_report_path(paths, args.label or f'{args.gwas_label}_{args.prs_label}_heldout')}")
        return 0

    if args.command == "run-clinical-prs-model":
        effective = apply_runtime_defaults(config)
        paths = build_output_paths(effective)
        label = args.label or clinical_prs_model_default_label(args.prs_label)
        outputs = run_clinical_prs_model(
            effective,
            paths,
            gwas_label=args.gwas_label,
            prs_label=args.prs_label,
            plink_prefix=args.plink_prefix,
            plink2_bin=args.plink2_bin,
            threshold_label=args.threshold_label,
            p_threshold=args.p_threshold,
            label=label,
            eligibility_flag=args.eligibility_flag,
            l2_penalty=args.l2_penalty,
            threads=args.threads,
            memory_mb=args.memory_mb,
        )
        metrics = outputs["metrics"]
        test_metrics = metrics[metrics["evaluation_set"] == "test"].iloc[0]
        print(f"Clinical+PRS model test ROC AUC: {test_metrics['roc_auc']:.4f}")
        print(f"Clinical+PRS model test average precision: {test_metrics['average_precision']:.4f}")
        print(f"metrics: {clinical_prs_model_metrics_path(paths, label)}")
        print(f"coefficients: {clinical_prs_model_coefficients_path(paths, label)}")
        print(f"predictions: {clinical_prs_model_predictions_path(paths, label)}")
        print(f"calibration: {clinical_prs_model_calibration_path(paths, label)}")
        print(f"report: {clinical_prs_model_report_path(paths, label)}")
        print(f"roc: {clinical_prs_model_roc_svg_path(paths, label)}")
        print(f"precision_recall: {clinical_prs_model_pr_svg_path(paths, label)}")
        print(f"calibration_plot: {clinical_prs_model_calibration_svg_path(paths, label)}")
        return 0

    if args.command == "diagnose-prs":
        effective = apply_runtime_defaults(config)
        paths = build_output_paths(effective)
        clinical_prs_label = args.clinical_prs_label or None
        outputs = run_prs_diagnostics(
            effective,
            paths,
            gwas_label=args.gwas_label,
            prs_label=args.prs_label,
            clinical_prs_label=clinical_prs_label,
            threshold_label=args.threshold_label,
            label=args.label,
            bootstrap_iterations=args.bootstrap_iterations,
            seed=args.seed,
        )
        print(f"PRS diagnostics overall rows: {outputs['overall'].shape[0]}")
        print(f"overall: {prs_diagnostics_overall_metrics_path(paths, args.label)}")
        print(f"bootstrap_ci: {prs_diagnostics_bootstrap_ci_path(paths, args.label)}")
        print(f"deciles: {prs_diagnostics_deciles_path(paths, args.label)}")
        print(f"ancestry: {prs_diagnostics_ancestry_path(paths, args.label)}")
        print(f"definite_sensitivity: {prs_diagnostics_definite_path(paths, args.label)}")
        print(f"cofactor_strata: {prs_diagnostics_cofactor_path(paths, args.label)}")
        print(f"calibration: {prs_diagnostics_calibration_path(paths, args.label)}")
        print(f"qc: {prs_diagnostics_qc_path(paths, args.label)}")
        print(f"report: {prs_diagnostics_report_path(paths, args.label)}")
        return 0

    if args.command == "profile-preindex-cases":
        effective, paths, cohort_df = _load_or_build_cohort_artifacts(config)
        windows = [value.strip() for value in args.windows.split(",") if value.strip()]
        outputs = profile_preindex_case_data(
            effective,
            cohort_df,
            paths,
            case_tier=args.case_tier,
            windows=windows,
            top_n=args.top_n,
        )
        print(f"Pre-index case profile rows: {outputs['summary'].shape[0]}")
        print(f"summary: {preindex_summary_path(paths)}")
        print(f"top_conditions: {preindex_condition_top_path(paths)}")
        print(f"top_measurements: {preindex_measurement_top_path(paths)}")
        print(f"biomarkers: {preindex_biomarker_path(paths)}")
        print(f"report: {preindex_report_path(paths)}")
        report_paths = render_existing_report(config)
        print(f"analysis_report: {report_paths.final_report_md}")
        print(f"analysis_dashboard: {report_paths.final_dashboard_html}")
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
