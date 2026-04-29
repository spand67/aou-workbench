"""Stage 4 Hail-native GWAS on the matched cohort."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import ProjectConfig
from .hail_utils import ensure_hail
from .io_utils import write_dataframe, write_json
from .paths import ProjectPaths, join_path
from .reporting import write_stage_report
from .sample_restriction import restrict_frame_for_gwas
from .statistics import bh_fdr, genomic_control_lambda
from .svg import write_manhattan_svg, write_qq_svg


def hail_stage4_full_results_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "stage4", "hail_gwas_full_results.tsv")


def hail_stage4_lead_hits_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "stage4", "hail_gwas_lead_hits.tsv")


def hail_stage4_qc_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "stage4", "hail_gwas_qc.json")


def hail_stage4_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "stage4", "hail_gwas_report.md")


def hail_stage4_manhattan_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "stage4", "hail_manhattan.svg")


def hail_stage4_qq_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "stage4", "hail_qq.svg")


def _split_acaf_mt_path(path: str) -> str:
    if "/multiMT/" in path:
        return path.replace("/multiMT/", "/splitMT/")
    return path


def _normalize_chromosomes(chromosomes: list[str] | None) -> list[str]:
    if not chromosomes:
        return [str(chrom) for chrom in range(1, 23)]
    normalized: list[str] = []
    seen: set[str] = set()
    for chromosome in chromosomes:
        bare = chromosome.replace("chr", "").strip()
        if not bare:
            continue
        if bare not in seen:
            normalized.append(bare)
            seen.add(bare)
    return normalized


def _hail_sample_frame(
    matched_df: pd.DataFrame,
    config: ProjectConfig,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    restricted = restrict_frame_for_gwas(config, matched_df, require_wgs=True).copy()
    restricted["person_id"] = restricted["person_id"].astype(str)

    stage = config.analysis.stage4
    covariates = [column for column in stage.covariates if column != config.analysis.matched_outcome_column]
    if "ancestry_pred" in config.cohort.exact_match_columns and "ancestry_pred" not in covariates:
        covariates.append("ancestry_pred")

    columns = ["person_id", config.analysis.matched_outcome_column, *[column for column in covariates if column in restricted.columns]]
    sample = restricted[columns].drop_duplicates(subset=["person_id"]).copy()

    numeric_covariates = [column for column in covariates if column != "ancestry_pred" and column in sample.columns]
    for column in [config.analysis.matched_outcome_column, *numeric_covariates]:
        if column in sample.columns:
            sample[column] = pd.to_numeric(sample[column], errors="coerce")

    hail_covariates = list(numeric_covariates)
    if "ancestry_pred" in sample.columns:
        sample["ancestry_pred"] = sample["ancestry_pred"].replace("", pd.NA)
        ancestry_dummies = pd.get_dummies(sample["ancestry_pred"], prefix="ancestry", dtype=float)
        if not ancestry_dummies.empty:
            reference = sorted(ancestry_dummies.columns)[0]
            ancestry_dummies = ancestry_dummies.drop(columns=[reference])
            for column in ancestry_dummies.columns:
                sample[column] = ancestry_dummies[column]
            hail_covariates.extend(ancestry_dummies.columns.tolist())

    required = ["person_id", config.analysis.matched_outcome_column, *hail_covariates]
    complete = sample.dropna(subset=[column for column in required if column in sample.columns]).copy()
    complete["s"] = complete["person_id"].astype(str)
    return complete, hail_covariates, covariates


def _lead_hit_subset(df: pd.DataFrame, window_bp: int) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    chosen: list[pd.Series] = []
    for row in df.sort_values("regression_p").itertuples(index=False):
        keep = True
        for hit in chosen:
            if str(hit.chromosome) == str(row.chromosome) and abs(int(hit.position) - int(row.position)) <= window_bp:
                keep = False
                break
        if keep:
            chosen.append(pd.Series(row._asdict()))
    return pd.DataFrame(chosen)


def _annotation_frame(config: ProjectConfig) -> pd.DataFrame:
    stage = config.analysis.stage4
    if stage is None or not stage.annotation_table:
        return pd.DataFrame()
    path = Path(stage.annotation_table)
    if not path.exists():
        return pd.DataFrame()
    annotation = pd.read_csv(path, sep="\t")
    if stage.variant_id_column in annotation.columns:
        annotation["variant_id"] = annotation[stage.variant_id_column].astype(str)
    return annotation


def _postprocess_results(
    result_df: pd.DataFrame,
    config: ProjectConfig,
    paths: ProjectPaths,
    *,
    matched_samples_requested: int,
    matched_samples_analyzed: int,
    covariates_used: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    stage = config.analysis.stage4
    if stage is None:
        return pd.DataFrame(), pd.DataFrame()

    full = result_df.copy()
    if not full.empty:
        full = full.sort_values(["regression_p", "variant_id"]).reset_index(drop=True)
        full["fdr_q"] = bh_fdr(full["regression_p"].values)
        full["minus_log10_p"] = -np.log10(full["regression_p"].clip(lower=1e-300))
        annotation = _annotation_frame(config)
        if not annotation.empty:
            full = full.merge(annotation.drop_duplicates("variant_id"), on="variant_id", how="left", suffixes=("", "_annotation"))

    lead_hits = _lead_hit_subset(full, stage.lead_hit_window_bp)
    write_dataframe(full, hail_stage4_full_results_path(paths))
    write_dataframe(lead_hits, hail_stage4_lead_hits_path(paths))
    write_manhattan_svg(full, hail_stage4_manhattan_path(paths))
    write_qq_svg(full, hail_stage4_qq_path(paths))
    write_json(
        {
            "chromosomes_tested": sorted(full["chromosome"].astype(str).unique().tolist()) if not full.empty else [],
            "matched_samples_requested": int(matched_samples_requested),
            "matched_samples_analyzed": int(matched_samples_analyzed),
            "n_variants_tested": int(full.shape[0]),
            "n_lead_hits": int(lead_hits.shape[0]),
            "lambda_gc": genomic_control_lambda(full["regression_p"].values if not full.empty else []),
            "covariates_used": covariates_used,
        },
        hail_stage4_qc_path(paths),
    )
    write_stage_report(
        title="Stage 4: Hail GWAS",
        summary_lines=[
            f"- Matched samples requested: {matched_samples_requested}",
            f"- Matched complete-case samples analyzed: {matched_samples_analyzed}",
            f"- Variants tested: {full.shape[0]}",
            f"- Lead hits: {lead_hits.shape[0]}",
            f"- MAF threshold: {stage.min_maf}",
            f"- Covariates: {', '.join(covariates_used) if covariates_used else 'intercept only'}",
        ],
        preview_df=lead_hits if not lead_hits.empty else full,
        preview_columns=["variant_id", "chromosome", "position", "regression_p", "fdr_q", "beta", "odds_ratio"],
        path=hail_stage4_report_path(paths),
    )
    return full, lead_hits


def run_stage4_hail_gwas(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
    *,
    chromosomes: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    stage = config.analysis.stage4
    if stage is None:
        return pd.DataFrame(), pd.DataFrame()

    chromosome_values = _normalize_chromosomes(chromosomes)
    sample_df, hail_covariates, raw_covariates = _hail_sample_frame(matched_df, config)
    if sample_df.empty:
        return _postprocess_results(
            pd.DataFrame(),
            config,
            paths,
            matched_samples_requested=int(matched_df["person_id"].astype(str).nunique()),
            matched_samples_analyzed=0,
            covariates_used=raw_covariates,
        )

    hl = ensure_hail(
        requester_pays_project=config.workbench.requester_pays_project,
        requester_pays_buckets=config.workbench.requester_pays_buckets,
    )
    acaf_mt_path = _split_acaf_mt_path(config.workbench.acaf_mt_path)
    print(f"Reading ACAF MT from {acaf_mt_path}", flush=True)
    mt = hl.read_matrix_table(acaf_mt_path)

    sample_ht = hl.Table.from_pandas(sample_df).key_by("s")
    intervals = [
        hl.parse_locus_interval(f"chr{chromosome}", reference_genome="GRCh38")
        for chromosome in chromosome_values
    ]
    mt = hl.filter_intervals(mt, intervals)
    mt = mt.semi_join_cols(sample_ht)
    mt = mt.annotate_cols(sample=sample_ht[mt.s])
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mt = mt.annotate_rows(call_stats=hl.agg.call_stats(mt.GT, mt.alleles))
    mt = mt.annotate_rows(
        alt_ac=hl.int64(mt.call_stats.AC[1]),
        alt_af=hl.float64(mt.call_stats.AF[1]),
        minor_allele_count=hl.int64(hl.min(mt.call_stats.AC[1], mt.call_stats.AN - mt.call_stats.AC[1])),
        maf=hl.float64(hl.min(mt.call_stats.AF[1], 1.0 - mt.call_stats.AF[1])),
        variant_id=(
            hl.str(mt.locus.contig).replace("chr", "")
            + ":"
            + hl.str(mt.locus.position)
            + ":"
            + mt.alleles[0]
            + ":"
            + mt.alleles[1]
        ),
    )
    mt = mt.filter_rows(
        hl.is_defined(mt.maf)
        & hl.is_defined(mt.minor_allele_count)
        & (mt.maf >= stage.min_maf)
        & (mt.minor_allele_count >= int(stage.min_minor_allele_count))
    )

    print(
        f"Running Hail GWAS on chromosomes {', '.join(chromosome_values)} with {sample_df['s'].nunique()} complete-case samples.",
        flush=True,
    )
    result_ht = hl.logistic_regression_rows(
        test="wald",
        y=hl.float64(mt.sample[config.analysis.matched_outcome_column]),
        x=hl.float64(mt.GT.n_alt_alleles()),
        covariates=[1.0, *[hl.float64(mt.sample[column]) for column in hail_covariates]],
        pass_through=[
            mt.variant_id,
            mt.locus.contig,
            mt.locus.position,
            mt.maf,
            mt.minor_allele_count,
        ],
    )
    result_ht = result_ht.select(
        variant_id=result_ht.variant_id,
        chromosome=result_ht.contig,
        position=result_ht.position,
        maf=result_ht.maf,
        minor_allele_count=result_ht.minor_allele_count,
        beta=result_ht.beta,
        se=result_ht.standard_error,
        z_stat=result_ht.z_stat,
        regression_p=result_ht.p_value,
        odds_ratio=hl.exp(result_ht.beta),
    )
    result_df = result_ht.to_pandas()
    if not result_df.empty:
        result_df["chromosome"] = result_df["chromosome"].astype(str).str.replace("chr", "", regex=False)
        result_df["position"] = pd.to_numeric(result_df["position"], errors="coerce")
        result_df["regression_p"] = pd.to_numeric(result_df["regression_p"], errors="coerce")
        result_df["beta"] = pd.to_numeric(result_df["beta"], errors="coerce")
        result_df["se"] = pd.to_numeric(result_df["se"], errors="coerce")
        result_df["odds_ratio"] = pd.to_numeric(result_df["odds_ratio"], errors="coerce")
        result_df["maf"] = pd.to_numeric(result_df["maf"], errors="coerce")
        result_df["minor_allele_count"] = pd.to_numeric(result_df["minor_allele_count"], errors="coerce")
        result_df["n_samples"] = int(sample_df["s"].nunique())

    return _postprocess_results(
        result_df,
        config,
        paths,
        matched_samples_requested=int(matched_df["person_id"].astype(str).nunique()),
        matched_samples_analyzed=int(sample_df["s"].nunique()),
        covariates_used=raw_covariates,
    )


__all__ = [
    "hail_stage4_full_results_path",
    "hail_stage4_lead_hits_path",
    "hail_stage4_manhattan_path",
    "hail_stage4_qc_path",
    "hail_stage4_qq_path",
    "hail_stage4_report_path",
    "run_stage4_hail_gwas",
]
