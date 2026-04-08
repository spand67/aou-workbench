"""Output path planning for the AoU workbench pipeline."""

from __future__ import annotations

from dataclasses import dataclass

from .config import ProjectConfig
from .io_utils import slugify


PROJECT_ROOT = __file__


def project_path(*parts: str) -> str:
    import os
    from pathlib import Path

    root = Path(__file__).resolve().parents[2]
    return str(root.joinpath(*parts))


def join_path(root: str, *parts: str) -> str:
    clean_root = root.rstrip("/")
    clean_parts = [part.strip("/") for part in parts if part]
    if clean_root.startswith("gs://"):
        return "/".join([clean_root, *clean_parts])
    from pathlib import Path

    path = Path(clean_root)
    for part in clean_parts:
        path = path / part
    return str(path)


@dataclass(frozen=True)
class ProjectPaths:
    run_root: str
    manifest_json: str
    cohort_sql_root: str
    built_cohort_tsv: str
    matched_cohort_tsv: str
    cohort_qc_json: str
    stage1_results_tsv: str
    stage1_qc_json: str
    stage1_report_md: str
    stage2_variant_tsv: str
    stage2_gene_tsv: str
    stage2_person_tsv: str
    stage2_qc_json: str
    stage2_report_md: str
    stage3_results_tsv: str
    stage3_qc_json: str
    stage3_report_md: str
    stage4_full_results_tsv: str
    stage4_lead_hits_tsv: str
    stage4_qc_json: str
    stage4_report_md: str
    stage4_manhattan_svg: str
    stage4_qq_svg: str
    final_report_md: str

    def as_dict(self) -> dict[str, str]:
        return self.__dict__.copy()


def build_output_paths(config: ProjectConfig) -> ProjectPaths:
    run_slug = slugify(config.analysis.analysis_name)
    run_root = join_path(config.analysis.output_dir, run_slug)
    return ProjectPaths(
        run_root=run_root,
        manifest_json=join_path(run_root, "run_manifest.json"),
        cohort_sql_root=join_path(run_root, "sql"),
        built_cohort_tsv=join_path(run_root, "cohort", "built_cohort.tsv"),
        matched_cohort_tsv=join_path(run_root, "cohort", "matched_cohort.tsv"),
        cohort_qc_json=join_path(run_root, "cohort", "cohort_qc.json"),
        stage1_results_tsv=join_path(run_root, "stage1", "prior_variant_results.tsv"),
        stage1_qc_json=join_path(run_root, "stage1", "prior_variant_qc.json"),
        stage1_report_md=join_path(run_root, "stage1", "prior_variant_report.md"),
        stage2_variant_tsv=join_path(run_root, "stage2", "plp_variant_results.tsv"),
        stage2_gene_tsv=join_path(run_root, "stage2", "plp_gene_results.tsv"),
        stage2_person_tsv=join_path(run_root, "stage2", "plp_person_results.tsv"),
        stage2_qc_json=join_path(run_root, "stage2", "plp_qc.json"),
        stage2_report_md=join_path(run_root, "stage2", "plp_report.md"),
        stage3_results_tsv=join_path(run_root, "stage3", "burden_results.tsv"),
        stage3_qc_json=join_path(run_root, "stage3", "burden_qc.json"),
        stage3_report_md=join_path(run_root, "stage3", "burden_report.md"),
        stage4_full_results_tsv=join_path(run_root, "stage4", "gwas_full_results.tsv"),
        stage4_lead_hits_tsv=join_path(run_root, "stage4", "gwas_lead_hits.tsv"),
        stage4_qc_json=join_path(run_root, "stage4", "gwas_qc.json"),
        stage4_report_md=join_path(run_root, "stage4", "gwas_report.md"),
        stage4_manhattan_svg=join_path(run_root, "stage4", "manhattan.svg"),
        stage4_qq_svg=join_path(run_root, "stage4", "qq.svg"),
        final_report_md=join_path(run_root, "rhabdo_summary_report.md"),
    )


__all__ = ["PROJECT_ROOT", "ProjectPaths", "build_output_paths", "join_path", "project_path"]
