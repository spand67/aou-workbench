"""CLI entrypoint for the AoU workbench."""

from __future__ import annotations

import argparse
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
from .gwas_workflow import prepare_terminal_gwas_workspace
from .io_utils import read_table
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
from .preflight import apply_runtime_defaults, format_preflight_report, run_preflight_checks
from .regenie import prepare_regenie_inputs
from .stage1_prepare import prepare_stage1_variant_table, prepare_wgs_sample_manifest
from .stage1_prior_variants import run_stage1_prior_variants
from .stage2_prepare import prepare_stage2_variant_table
from .stage2_plp_panel import run_stage2_plp_panel
from .stage4_hail_gwas import (
    hail_pilot_lead_hits_path,
    hail_pilot_manhattan_path,
    hail_pilot_qc_path,
    hail_pilot_qq_path,
    hail_pilot_report_path,
    hail_pilot_results_path,
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
    if os.path.exists(paths.matched_cohort_tsv) and not require_wgs:
        matched_df = read_table(paths.matched_cohort_tsv)
        matched_df = apply_time_anchored_clinical_cofactors(effective, matched_df)
        return effective, paths, matched_df
    return match_controls_artifacts(config, require_wgs=require_wgs)


def _load_or_build_cohort_artifacts(config, *, require_wgs: bool = False):
    effective = apply_runtime_defaults(config)
    paths = build_output_paths(effective)
    if os.path.exists(paths.built_cohort_tsv) and not require_wgs:
        cohort_df = read_table(paths.built_cohort_tsv)
        return effective, paths, cohort_df
    return build_cohort_artifacts(config, require_wgs=require_wgs)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="AoU workbench CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    preflight_parser = subparsers.add_parser("preflight", help="Run runtime and input validation.")
    _add_config_arguments(preflight_parser)

    build_parser = subparsers.add_parser("build-cohort", help="Build the tiered rhabdomyolysis cohort.")
    _add_config_arguments(build_parser)
    build_parser.add_argument("--require-wgs", action="store_true", help="Restrict the saved cohort to WGS/ACAF sample IDs.")

    match_parser = subparsers.add_parser("match-controls", help="Build the matched case-control cohort.")
    _add_config_arguments(match_parser)
    match_parser.add_argument("--require-wgs", action="store_true", help="Restrict the cohort and matching universe to WGS/ACAF sample IDs.")

    wgs_manifest_parser = subparsers.add_parser(
        "prepare-wgs-manifest",
        help="Write the WGS/ACAF sample manifest used for WGS-restricted cohort runs.",
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
        default="acaf_chr22_maf05_train_qc",
        help="Output label under stage4/hail_pilot/. Default: acaf_chr22_maf05_train_qc.",
    )

    for name in ("run-stage1", "run-stage2", "run-stage3", "run-stage4", "run-all"):
        stage_parser = subparsers.add_parser(name, help=f"Execute {name}.")
        _add_config_arguments(stage_parser)
        stage_parser.add_argument("--skip-preflight", action="store_true")

    report_parser = subparsers.add_parser("report", help="Rebuild the final markdown report from existing outputs.")
    _add_config_arguments(report_parser)

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
        effective, paths, cohort_df = _load_or_build_cohort_artifacts(config)
        if not os.path.exists(clinical_model_input_path(paths)):
            _, _, matched_df = _load_or_build_matched_artifacts(config)
            characterize_case_control_cohort(effective, cohort_df, matched_df, paths)
        matched_df = read_table(clinical_model_input_path(paths))
        chromosomes = [value.strip() for value in args.chromosomes.split(",") if value.strip()]
        full, hits = run_stage4_hail_pilot_gwas(
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
            label=args.label,
        )
        print(f"Hail pilot GWAS variants tested: {full.shape[0]}")
        print(f"Hail pilot GWAS lead hits: {hits.shape[0]}")
        print(f"results: {hail_pilot_results_path(paths, args.label)}")
        print(f"lead_hits: {hail_pilot_lead_hits_path(paths, args.label)}")
        print(f"qc: {hail_pilot_qc_path(paths, args.label)}")
        print(f"variant_qc_summary: {hail_pilot_variant_qc_summary_path(paths, args.label)}")
        print(f"report: {hail_pilot_report_path(paths, args.label)}")
        print(f"manhattan: {hail_pilot_manhattan_path(paths, args.label)}")
        print(f"qq: {hail_pilot_qq_path(paths, args.label)}")
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
        if args.require_wgs or not os.path.exists(clinical_model_input_path(paths)):
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
