"""CLI entrypoint for the AoU workbench."""

from __future__ import annotations

import argparse
import os
import sys

from .config import load_project_config
from .io_utils import read_table
from .paths import build_output_paths, project_path
from .pipeline import build_cohort_artifacts, match_controls_artifacts, render_existing_report, run_all
from .preflight import apply_runtime_defaults, format_preflight_report, run_preflight_checks
from .regenie import prepare_regenie_inputs
from .stage1_prepare import prepare_stage1_variant_table
from .stage1_prior_variants import run_stage1_prior_variants
from .stage2_prepare import prepare_stage2_variant_table
from .stage2_plp_panel import run_stage2_plp_panel
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


def _load_or_build_matched_artifacts(config):
    effective = apply_runtime_defaults(config)
    paths = build_output_paths(effective)
    if os.path.exists(paths.matched_cohort_tsv):
        matched_df = read_table(paths.matched_cohort_tsv)
        return effective, paths, matched_df
    return match_controls_artifacts(config)


def _load_or_build_cohort_artifacts(config):
    effective = apply_runtime_defaults(config)
    paths = build_output_paths(effective)
    if os.path.exists(paths.built_cohort_tsv):
        cohort_df = read_table(paths.built_cohort_tsv)
        return effective, paths, cohort_df
    return build_cohort_artifacts(config)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="AoU workbench CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    preflight_parser = subparsers.add_parser("preflight", help="Run runtime and input validation.")
    _add_config_arguments(preflight_parser)

    build_parser = subparsers.add_parser("build-cohort", help="Build the tiered rhabdomyolysis cohort.")
    _add_config_arguments(build_parser)

    match_parser = subparsers.add_parser("match-controls", help="Build the matched case-control cohort.")
    _add_config_arguments(match_parser)

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

    prepare_stage4_parser = subparsers.add_parser(
        "prepare-stage4",
        help="Subset the ACAF smaller callset with Hail, export local PLINK files, and run a Hail pilot GWAS.",
    )
    _add_config_arguments(prepare_stage4_parser)
    prepare_stage4_parser.add_argument("--chromosome", default="chr19")

    for name in ("run-stage1", "run-stage2", "run-stage3", "run-stage4", "run-all"):
        stage_parser = subparsers.add_parser(name, help=f"Execute {name}.")
        _add_config_arguments(stage_parser)
        stage_parser.add_argument("--skip-preflight", action="store_true")

    report_parser = subparsers.add_parser("report", help="Rebuild the final markdown report from existing outputs.")
    _add_config_arguments(report_parser)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    config = _load_config(args)

    if args.command == "preflight":
        print(format_preflight_report(run_preflight_checks(config)))
        return 0

    if args.command == "build-cohort":
        _, paths, cohort_df = build_cohort_artifacts(config)
        print(f"Built cohort rows: {len(cohort_df)}")
        print(f"Cohort path: {paths.built_cohort_tsv}")
        return 0

    if args.command == "match-controls":
        _, paths, matched_df = match_controls_artifacts(config)
        print(f"Matched cohort rows: {len(matched_df)}")
        print(f"Matched cohort path: {paths.matched_cohort_tsv}")
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
        print(f"Prepared REGENIE inputs for {matched_df['person_id'].astype(str).nunique()} matched samples.")
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

    if args.command == "run-all":
        paths = run_all(config, skip_preflight=args.skip_preflight)
        print(f"Run root: {paths.run_root}")
        print(f"Final report: {paths.final_report_md}")
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
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
