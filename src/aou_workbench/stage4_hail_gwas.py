"""Stage 4 Hail-native GWAS on the matched cohort."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import ProjectConfig
from .hail_utils import ensure_hail
from .io_utils import slugify, write_dataframe, write_json
from .paths import ProjectPaths, join_path
from .reporting import write_stage_report
from .sample_restriction import has_wgs_manifest, restrict_frame_for_gwas, restrict_frame_to_ids, wgs_present_ids
from .statistics import bh_fdr, genomic_control_lambda
from .svg import write_manhattan_svg, write_qq_svg


PILOT_COVARIATES = ("age_at_index", "is_female", "pc1", "pc2", "pc3", "pc4", "pc5")
PILOT_CASE_CONTROL_REQUIRED_COLUMNS = ("match_group_id", "case_tier", "eligible_control")
PILOT_GENOTYPE_SOURCES = ("acaf", "microarray")


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


def hail_pilot_output_dir(paths: ProjectPaths, label: str) -> str:
    return join_path(paths.run_root, "stage4", "hail_pilot", slugify(label))


def hail_pilot_results_path(paths: ProjectPaths, label: str) -> str:
    return join_path(hail_pilot_output_dir(paths, label), "gwas_results.tsv")


def hail_pilot_lead_hits_path(paths: ProjectPaths, label: str) -> str:
    return join_path(hail_pilot_output_dir(paths, label), "lead_hits.tsv")


def hail_pilot_qc_path(paths: ProjectPaths, label: str) -> str:
    return join_path(hail_pilot_output_dir(paths, label), "gwas_qc.json")


def hail_pilot_variant_qc_summary_path(paths: ProjectPaths, label: str) -> str:
    return join_path(hail_pilot_output_dir(paths, label), "variant_qc_summary.tsv")


def hail_pilot_manhattan_path(paths: ProjectPaths, label: str) -> str:
    return join_path(hail_pilot_output_dir(paths, label), "manhattan.svg")


def hail_pilot_qq_path(paths: ProjectPaths, label: str) -> str:
    return join_path(hail_pilot_output_dir(paths, label), "qq.svg")


def hail_pilot_report_path(paths: ProjectPaths, label: str) -> str:
    return join_path(hail_pilot_output_dir(paths, label), "report.md")


def _split_acaf_mt_path(path: str) -> str:
    if "/multiMT/" in path:
        return path.replace("/multiMT/", "/splitMT/")
    return path


def hail_pilot_default_label(
    genotype_source: str,
    chromosomes: list[str] | None,
    *,
    min_maf: float = 0.05,
) -> str:
    source = _normalize_genotype_source(genotype_source)
    chromosome_values = _normalize_autosomal_chromosomes(chromosomes or ["22"])
    chromosome_label = "chr" + "-".join(chromosome_values)
    maf_text = f"{min_maf:g}"
    maf_value = maf_text[2:] if maf_text.startswith("0.") else maf_text.replace(".", "")
    maf_label = f"maf{maf_value}"
    return f"{source}_{chromosome_label}_{maf_label}_train_qc"


def _normalize_genotype_source(genotype_source: str) -> str:
    source = genotype_source.strip().lower()
    if source not in PILOT_GENOTYPE_SOURCES:
        raise ValueError(f"Unsupported Hail pilot genotype source `{genotype_source}`; expected one of: {', '.join(PILOT_GENOTYPE_SOURCES)}")
    return source


def _pilot_genotype_mt_path(config: ProjectConfig, genotype_source: str) -> str:
    source = _normalize_genotype_source(genotype_source)
    if source == "microarray":
        return config.workbench.microarray_mt_path
    return _split_acaf_mt_path(config.workbench.acaf_mt_path)


def _pilot_restricts_to_wgs_manifest(genotype_source: str) -> bool:
    return _normalize_genotype_source(genotype_source) == "acaf"


def _key_matrix_table_by_sample_id(mt, hl):
    if "s" in mt.col:
        return mt.key_cols_by("s")
    col_key_fields = list(mt.col_key)
    if not col_key_fields:
        raise RuntimeError("Genotype MatrixTable has no column key and no `s` sample field.")
    sample_field = col_key_fields[0]
    return mt.annotate_cols(s=hl.str(getattr(mt, sample_field))).key_cols_by("s")


def _matrix_table_gt(mt, hl):
    if "GT" in mt.entry:
        return mt.GT
    if "LGT" in mt.entry and "LA" in mt.entry:
        return hl.vds.lgt_to_gt(mt.LGT, mt.LA)
    raise RuntimeError("Genotype MatrixTable must contain GT, or LGT plus LA, for Hail pilot GWAS.")


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


def _normalize_autosomal_chromosomes(chromosomes: list[str] | None) -> list[str]:
    values = _normalize_chromosomes(chromosomes)
    invalid = [value for value in values if not value.isdigit() or int(value) < 1 or int(value) > 22]
    if invalid:
        raise ValueError(f"Hail pilot GWAS only supports autosomes 1-22; got: {', '.join(invalid)}")
    return values


def _is_binary_series(series: pd.Series) -> bool:
    observed = set(series.dropna().astype(float).unique())
    return observed.issubset({0.0, 1.0})


def _stable_hail_covariates(
    sample: pd.DataFrame,
    outcome_column: str,
    candidate_covariates: list[str],
) -> tuple[list[str], list[str]]:
    kept: list[str] = []
    dropped: list[str] = []
    outcome = sample[outcome_column].astype(float)

    for column in candidate_covariates:
        if column not in sample.columns:
            continue
        values = sample[column]
        if values.dropna().nunique() <= 1:
            dropped.append(column)
            continue
        if _is_binary_series(values):
            contingency = (
                pd.crosstab(values.astype(float), outcome)
                .reindex(index=[0.0, 1.0], columns=[0.0, 1.0], fill_value=0)
            )
            if (contingency.loc[1.0, 0.0] == 0) or (contingency.loc[1.0, 1.0] == 0):
                dropped.append(column)
                continue
        kept.append(column)
    return kept, dropped


def _flag_mask(frame: pd.DataFrame, column: str) -> pd.Series:
    values = frame[column]
    numeric = pd.to_numeric(values, errors="coerce").fillna(0)
    string = values.astype("string").str.lower().str.strip()
    return numeric.ne(0) | string.isin({"true", "t", "yes", "y"})


def _pilot_case_control_definition_mask(sample: pd.DataFrame, outcome_column: str) -> pd.Series:
    required = set(PILOT_CASE_CONTROL_REQUIRED_COLUMNS)
    missing = sorted(required.difference(sample.columns))
    if missing:
        raise RuntimeError(
            "GWAS pilot input is missing required case-control columns: "
            + ", ".join(missing)
            + ". Existing outputs are stale; rerun `aou-workbench characterize-cohort`, or rebuild cohort/matching if this persists."
        )

    outcome = pd.to_numeric(sample[outcome_column], errors="coerce")
    tier = sample["case_tier"].astype("string").str.lower()

    if "broad_rhabdo_case" in sample.columns:
        broad_case = _flag_mask(sample, "broad_rhabdo_case")
    else:
        broad_case = tier.isin({"broad", "definite"})

    case_mask = outcome.eq(1) & broad_case
    control_mask = outcome.eq(0) & tier.eq("control") & _flag_mask(sample, "eligible_control")
    return case_mask | control_mask


def _qc_row(name: str, threshold: str, before: int, after: int) -> dict[str, object]:
    before_int = int(before)
    after_int = int(after)
    return {
        "filter": name,
        "threshold": threshold,
        "rows_before": before_int,
        "rows_after": after_int,
        "rows_removed": before_int - after_int,
    }


def _variant_qc_summary_rows(
    *,
    chromosomes: list[str],
    initial_rows: int,
    biallelic_rows: int,
    maf_rows: int,
    mac_rows: int,
    call_rate_rows: int,
    hwe_rows: int,
    min_maf: float,
    min_mac: int,
    min_call_rate: float,
    hwe_p_control: float,
) -> list[dict[str, object]]:
    return [
        _qc_row("interval_and_sample_restriction", ",".join(chromosomes), initial_rows, initial_rows),
        _qc_row("autosomal_biallelic_snp", "biallelic A/C/G/T SNP", initial_rows, biallelic_rows),
        _qc_row("maf", f">= {min_maf}", biallelic_rows, maf_rows),
        _qc_row("minor_allele_count", f">= {int(min_mac)}", maf_rows, mac_rows),
        _qc_row("call_rate", f">= {min_call_rate}", mac_rows, call_rate_rows),
        _qc_row("control_hwe_p", f">= {hwe_p_control}", call_rate_rows, hwe_rows),
    ]


def _hail_sample_frame(
    matched_df: pd.DataFrame,
    config: ProjectConfig,
) -> tuple[pd.DataFrame, list[str], list[str], list[str]]:
    restricted = restrict_frame_for_gwas(config, matched_df, require_wgs=True).copy()
    restricted["person_id"] = restricted["person_id"].astype(str)

    stage = config.analysis.stage4
    covariates = [column for column in stage.covariates if column != config.analysis.matched_outcome_column]
    raw_covariates = list(covariates)
    if "ancestry_pred" in config.cohort.exact_match_columns and "ancestry_pred" not in raw_covariates:
        raw_covariates.append("ancestry_pred")

    columns = ["person_id", config.analysis.matched_outcome_column, *[column for column in raw_covariates if column in restricted.columns]]
    sample = restricted[columns].drop_duplicates(subset=["person_id"]).copy()

    numeric_covariates = [column for column in raw_covariates if column != "ancestry_pred" and column in sample.columns]
    for column in [config.analysis.matched_outcome_column, *numeric_covariates]:
        if column in sample.columns:
            sample[column] = pd.to_numeric(sample[column], errors="coerce")

    if "ancestry_pred" in sample.columns:
        sample["ancestry_pred"] = sample["ancestry_pred"].replace("", pd.NA)

    required = ["person_id", config.analysis.matched_outcome_column, *numeric_covariates]
    complete = sample.dropna(subset=[column for column in required if column in sample.columns]).copy()
    hail_covariates, dropped_covariates = _stable_hail_covariates(
        complete,
        config.analysis.matched_outcome_column,
        numeric_covariates,
    )
    complete = complete.dropna(
        subset=[column for column in [config.analysis.matched_outcome_column, *hail_covariates] if column in complete.columns]
    ).copy()
    complete["s"] = complete["person_id"].astype(str)
    return complete, hail_covariates, raw_covariates, dropped_covariates


def _hail_pilot_sample_frame(
    matched_df: pd.DataFrame,
    config: ProjectConfig,
    *,
    analysis_split: str = "train",
    eligibility_flag: str = "primary_model_eligible",
    covariates: tuple[str, ...] = PILOT_COVARIATES,
    restrict_to_wgs_manifest: bool = True,
) -> tuple[pd.DataFrame, list[str], list[str], list[str], dict[str, int], bool]:
    if "person_id" not in matched_df.columns:
        raise RuntimeError("Matched cohort must contain `person_id` for Hail pilot GWAS.")

    outcome_column = config.analysis.matched_outcome_column
    if outcome_column not in matched_df.columns:
        raise RuntimeError(f"Matched cohort must contain outcome column `{outcome_column}`.")

    sample = matched_df.copy()
    sample["person_id"] = sample["person_id"].astype(str)
    counts: dict[str, int] = {
        "matched_input_participants": int(sample["person_id"].nunique()),
    }

    if analysis_split:
        if "analysis_split" not in sample.columns:
            raise RuntimeError("Matched cohort must contain `analysis_split`; run `aou-workbench characterize-cohort` first.")
        sample = sample[sample["analysis_split"].astype(str) == str(analysis_split)].copy()
    counts["after_analysis_split_participants"] = int(sample["person_id"].nunique())

    if eligibility_flag:
        if eligibility_flag not in sample.columns:
            raise RuntimeError(f"Matched cohort must contain `{eligibility_flag}`; run `aou-workbench characterize-cohort` first.")
        eligibility = pd.to_numeric(sample[eligibility_flag], errors="coerce").fillna(0)
        sample = sample[eligibility == 1].copy()
    counts["after_eligibility_participants"] = int(sample["person_id"].nunique())

    sample = sample[_pilot_case_control_definition_mask(sample, outcome_column)].copy()
    outcome_after_definition = pd.to_numeric(sample[outcome_column], errors="coerce")
    counts["after_case_control_definition_participants"] = int(sample["person_id"].nunique())
    counts["after_case_control_definition_cases"] = int(outcome_after_definition.eq(1).sum())
    counts["after_case_control_definition_controls"] = int(outcome_after_definition.eq(0).sum())

    wgs_manifest_used = False
    if restrict_to_wgs_manifest and has_wgs_manifest(config):
        sample = restrict_frame_to_ids(sample, wgs_present_ids(config, require=True), id_column="person_id")
        wgs_manifest_used = True
    counts["after_wgs_manifest_participants"] = int(sample["person_id"].nunique())

    missing_covariates = [column for column in covariates if column not in sample.columns]
    if missing_covariates:
        raise RuntimeError(f"Matched cohort is missing required Hail pilot covariates: {', '.join(missing_covariates)}")

    raw_covariates = [column for column in covariates if column in sample.columns]
    columns = ["person_id", outcome_column, *raw_covariates]
    complete = sample[columns].drop_duplicates(subset=["person_id"]).copy()
    for column in [outcome_column, *raw_covariates]:
        complete[column] = pd.to_numeric(complete[column], errors="coerce")
    complete = complete.dropna(subset=[outcome_column, *raw_covariates]).copy()

    hail_covariates, dropped_covariates = _stable_hail_covariates(
        complete,
        outcome_column,
        raw_covariates,
    )
    complete = complete.dropna(subset=[outcome_column, *hail_covariates]).copy()
    complete["s"] = complete["person_id"].astype(str)
    counts["after_complete_case_participants"] = int(complete["person_id"].nunique())
    outcome_complete = pd.to_numeric(complete[outcome_column], errors="coerce")
    counts["after_complete_case_cases"] = int(outcome_complete.eq(1).sum())
    counts["after_complete_case_controls"] = int(outcome_complete.eq(0).sum())
    return complete, hail_covariates, raw_covariates, dropped_covariates, counts, wgs_manifest_used


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
    dropped_covariates: list[str],
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
            "dropped_covariates": dropped_covariates,
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
            (
                f"- Dropped unstable Hail covariates: {', '.join(dropped_covariates)}"
                if dropped_covariates
                else "- Dropped unstable Hail covariates: none"
            ),
        ],
        preview_df=lead_hits if not lead_hits.empty else full,
        preview_columns=["variant_id", "chromosome", "position", "regression_p", "fdr_q", "beta", "odds_ratio"],
        path=hail_stage4_report_path(paths),
    )
    return full, lead_hits


def _postprocess_pilot_results(
    result_df: pd.DataFrame,
    variant_qc_summary: pd.DataFrame,
    config: ProjectConfig,
    paths: ProjectPaths,
    *,
    label: str,
    chromosomes: list[str],
    analysis_split: str,
    eligibility_flag: str,
    min_maf: float,
    min_mac: int,
    min_call_rate: float,
    hwe_p_control: float,
    sample_counts: dict[str, int],
    wgs_manifest_used: bool,
    genotype_source: str,
    genotype_mt_path: str,
    covariates_used: list[str],
    dropped_covariates: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    stage = config.analysis.stage4
    if stage is None:
        return pd.DataFrame(), pd.DataFrame()

    full = result_df.copy()
    if not full.empty:
        full = full.sort_values(["regression_p", "variant_id"]).reset_index(drop=True)
        full["fdr_q"] = bh_fdr(full["regression_p"].values)
        full["minus_log10_p"] = -np.log10(full["regression_p"].clip(lower=1e-300))

    lead_hits = _lead_hit_subset(full, stage.lead_hit_window_bp)
    write_dataframe(full, hail_pilot_results_path(paths, label))
    write_dataframe(lead_hits, hail_pilot_lead_hits_path(paths, label))
    write_dataframe(variant_qc_summary, hail_pilot_variant_qc_summary_path(paths, label))
    write_manhattan_svg(full, hail_pilot_manhattan_path(paths, label))
    write_qq_svg(full, hail_pilot_qq_path(paths, label))
    write_json(
        {
            "pilot_label": label,
            "genotype_source": genotype_source,
            "genotype_mt_path": genotype_mt_path,
            "acaf_mt_path": _split_acaf_mt_path(config.workbench.acaf_mt_path),
            "chromosomes_tested": chromosomes,
            "analysis_split": analysis_split,
            "eligibility_flag": eligibility_flag,
            "wgs_manifest_used": wgs_manifest_used,
            "sample_counts": sample_counts,
            "n_variants_tested": int(full.shape[0]),
            "n_lead_hits": int(lead_hits.shape[0]),
            "lambda_gc": genomic_control_lambda(full["regression_p"].values if not full.empty else []),
            "min_maf": float(min_maf),
            "min_mac": int(min_mac),
            "min_call_rate": float(min_call_rate),
            "hwe_p_control": float(hwe_p_control),
            "covariates_used": covariates_used,
            "dropped_covariates": dropped_covariates,
        },
        hail_pilot_qc_path(paths, label),
    )
    write_stage_report(
        title="Stage 4: Hail Pilot GWAS",
        summary_lines=[
            f"- Pilot label: {label}",
            f"- Genotype source: {genotype_source}",
            f"- Genotype MatrixTable: {genotype_mt_path}",
            "- Case-control definition: broad rhabdomyolysis cases versus matched eligible controls",
            f"- Analysis split: {analysis_split}",
            f"- Eligibility flag: {eligibility_flag}",
            f"- Chromosomes: {', '.join(chromosomes)}",
            f"- Complete-case samples analyzed: {sample_counts.get('after_complete_case_participants', 0)}",
            (
                "- Complete-case cases/controls: "
                f"{sample_counts.get('after_complete_case_cases', 0)} / "
                f"{sample_counts.get('after_complete_case_controls', 0)}"
            ),
            f"- MatrixTable samples analyzed: {sample_counts.get('matrix_table_participants', 0)}",
            f"- Variants tested after QC: {full.shape[0]}",
            f"- Lead hits: {lead_hits.shape[0]}",
            f"- Variant QC: MAF >= {min_maf}, MAC >= {min_mac}, call rate >= {min_call_rate}, control HWE p >= {hwe_p_control}",
            f"- Covariates: {', '.join(covariates_used) if covariates_used else 'intercept only'}",
            (
                f"- Dropped unstable covariates: {', '.join(dropped_covariates)}"
                if dropped_covariates
                else "- Dropped unstable covariates: none"
            ),
        ],
        preview_df=lead_hits if not lead_hits.empty else full,
        preview_columns=["variant_id", "chromosome", "position", "regression_p", "fdr_q", "beta", "odds_ratio"],
        path=hail_pilot_report_path(paths, label),
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
    sample_df, hail_covariates, raw_covariates, dropped_covariates = _hail_sample_frame(matched_df, config)
    if sample_df.empty:
        return _postprocess_results(
            pd.DataFrame(),
            config,
            paths,
            matched_samples_requested=int(matched_df["person_id"].astype(str).nunique()),
            matched_samples_analyzed=0,
            covariates_used=raw_covariates,
            dropped_covariates=dropped_covariates,
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
        contig=hl.str(mt.locus.contig),
        position=hl.int32(mt.locus.position),
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
            mt.contig,
            mt.position,
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
        covariates_used=hail_covariates,
        dropped_covariates=dropped_covariates,
    )


def run_stage4_hail_pilot_gwas(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
    *,
    chromosomes: list[str] | None = None,
    min_maf: float = 0.05,
    min_mac: int = 20,
    min_call_rate: float = 0.98,
    hwe_p_control: float = 1e-6,
    analysis_split: str = "train",
    eligibility_flag: str = "primary_model_eligible",
    label: str | None = None,
    genotype_source: str = "acaf",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    stage = config.analysis.stage4
    if stage is None:
        return pd.DataFrame(), pd.DataFrame()

    source = _normalize_genotype_source(genotype_source)
    chromosome_values = _normalize_autosomal_chromosomes(chromosomes or ["22"])
    label = label or hail_pilot_default_label(source, chromosome_values, min_maf=min_maf)
    genotype_mt_path = _pilot_genotype_mt_path(config, source)
    restrict_to_wgs_manifest = _pilot_restricts_to_wgs_manifest(source)
    (
        sample_df,
        hail_covariates,
        raw_covariates,
        dropped_covariates,
        sample_counts,
        wgs_manifest_used,
    ) = _hail_pilot_sample_frame(
        matched_df,
        config,
        analysis_split=analysis_split,
        eligibility_flag=eligibility_flag,
        restrict_to_wgs_manifest=restrict_to_wgs_manifest,
    )
    empty_qc = pd.DataFrame(
        columns=["filter", "threshold", "rows_before", "rows_after", "rows_removed"]
    )
    if sample_df.empty:
        return _postprocess_pilot_results(
            pd.DataFrame(),
            empty_qc,
            config,
            paths,
            label=label,
            chromosomes=chromosome_values,
            analysis_split=analysis_split,
            eligibility_flag=eligibility_flag,
            min_maf=min_maf,
            min_mac=min_mac,
            min_call_rate=min_call_rate,
            hwe_p_control=hwe_p_control,
            sample_counts=sample_counts,
            wgs_manifest_used=wgs_manifest_used,
            genotype_source=source,
            genotype_mt_path=genotype_mt_path,
            covariates_used=raw_covariates,
            dropped_covariates=dropped_covariates,
        )

    hl = ensure_hail(
        requester_pays_project=config.workbench.requester_pays_project,
        requester_pays_buckets=config.workbench.requester_pays_buckets,
    )
    print(f"Reading {source} MT from {genotype_mt_path}", flush=True)
    mt = hl.read_matrix_table(genotype_mt_path)

    sample_ht = hl.Table.from_pandas(sample_df).key_by("s")
    intervals = [
        hl.parse_locus_interval(f"chr{chromosome}", reference_genome="GRCh38")
        for chromosome in chromosome_values
    ]
    mt = hl.filter_intervals(mt, intervals)
    mt = _key_matrix_table_by_sample_id(mt, hl)
    mt = mt.semi_join_cols(sample_ht)
    mt = mt.annotate_cols(sample=sample_ht[mt.s])
    print("Counting MatrixTable sample overlap after chromosome and sample filters.", flush=True)
    matrix_table_participants = int(mt.count_cols())
    sample_counts["matrix_table_participants"] = matrix_table_participants
    print(f"MatrixTable sample overlap: {matrix_table_participants}", flush=True)
    if matrix_table_participants == 0:
        return _postprocess_pilot_results(
            pd.DataFrame(),
            empty_qc,
            config,
            paths,
            label=label,
            chromosomes=chromosome_values,
            analysis_split=analysis_split,
            eligibility_flag=eligibility_flag,
            min_maf=min_maf,
            min_mac=min_mac,
            min_call_rate=min_call_rate,
            hwe_p_control=hwe_p_control,
            sample_counts=sample_counts,
            wgs_manifest_used=wgs_manifest_used,
            genotype_source=source,
            genotype_mt_path=genotype_mt_path,
            covariates_used=raw_covariates,
            dropped_covariates=dropped_covariates,
        )

    bases = hl.literal({"A", "C", "G", "T"})
    biallelic_snp = (
        (hl.len(mt.alleles) == 2)
        & (hl.len(mt.alleles[0]) == 1)
        & (hl.len(mt.alleles[1]) == 1)
        & bases.contains(mt.alleles[0])
        & bases.contains(mt.alleles[1])
    )
    print("Counting variant rows and autosomal biallelic SNP rows.", flush=True)
    row_counts = mt.aggregate_rows(
        hl.struct(
            initial=hl.agg.count(),
            biallelic=hl.agg.count_where(biallelic_snp),
        )
    )
    initial_rows = int(row_counts.initial)
    biallelic_rows = int(row_counts.biallelic)
    mt = mt.filter_rows(biallelic_snp)
    gt_expr = _matrix_table_gt(mt, hl)
    mt = mt.select_entries(_GT=gt_expr)

    outcome = config.analysis.matched_outcome_column
    mt = mt.annotate_rows(call_stats=hl.agg.call_stats(mt._GT, mt.alleles))
    mt = mt.annotate_rows(
        n_called=hl.int64(hl.agg.count_where(hl.is_defined(mt._GT))),
        n_total=hl.int64(hl.agg.count()),
        hwe_control=hl.agg.filter(
            hl.float64(mt.sample[outcome]) == 0.0,
            hl.agg.hardy_weinberg_test(mt._GT),
        ),
    )
    mt = mt.annotate_rows(
        contig=hl.str(mt.locus.contig),
        position=hl.int32(mt.locus.position),
        alt_ac=hl.int64(mt.call_stats.AC[1]),
        alt_af=hl.float64(mt.call_stats.AF[1]),
        minor_allele_count=hl.int64(hl.min(mt.call_stats.AC[1], mt.call_stats.AN - mt.call_stats.AC[1])),
        maf=hl.float64(hl.min(mt.call_stats.AF[1], 1.0 - mt.call_stats.AF[1])),
        call_rate=hl.float64(mt.n_called) / hl.float64(mt.n_total),
        hwe_control_p=mt.hwe_control.p_value,
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

    maf_filter = hl.is_defined(mt.maf) & (mt.maf >= float(min_maf))
    mac_filter = maf_filter & hl.is_defined(mt.minor_allele_count) & (mt.minor_allele_count >= int(min_mac))
    call_rate_filter = mac_filter & hl.is_defined(mt.call_rate) & (mt.call_rate >= float(min_call_rate))
    final_qc_filter = call_rate_filter & hl.is_defined(mt.hwe_control_p) & (mt.hwe_control_p >= float(hwe_p_control))
    print("Computing cumulative variant QC counts.", flush=True)
    qc_counts = mt.aggregate_rows(
        hl.struct(
            maf=hl.agg.count_where(maf_filter),
            mac=hl.agg.count_where(mac_filter),
            call_rate=hl.agg.count_where(call_rate_filter),
            hwe=hl.agg.count_where(final_qc_filter),
        )
    )
    qc_rows = _variant_qc_summary_rows(
        chromosomes=chromosome_values,
        initial_rows=initial_rows,
        biallelic_rows=biallelic_rows,
        maf_rows=int(qc_counts.maf),
        mac_rows=int(qc_counts.mac),
        call_rate_rows=int(qc_counts.call_rate),
        hwe_rows=int(qc_counts.hwe),
        min_maf=min_maf,
        min_mac=min_mac,
        min_call_rate=min_call_rate,
        hwe_p_control=hwe_p_control,
    )
    mt = mt.filter_rows(final_qc_filter)
    mt = mt.annotate_rows(
        hwe_case=hl.agg.filter(
            hl.float64(mt.sample[outcome]) == 1.0,
            hl.agg.hardy_weinberg_test(mt._GT),
        )
    )
    mt = mt.annotate_rows(
        hwe_case_p=mt.hwe_case.p_value,
    )

    print(
        "Running Hail pilot GWAS on "
        f"{source} chromosomes {', '.join(chromosome_values)} with {matrix_table_participants} MatrixTable samples.",
        flush=True,
    )
    result_ht = hl.logistic_regression_rows(
        test="wald",
        y=hl.float64(mt.sample[outcome]),
        x=hl.float64(mt._GT.n_alt_alleles()),
        covariates=[1.0, *[hl.float64(mt.sample[column]) for column in hail_covariates]],
        pass_through=[
            mt.variant_id,
            mt.contig,
            mt.position,
            mt.maf,
            mt.minor_allele_count,
            mt.call_rate,
            mt.hwe_control_p,
            mt.hwe_case_p,
        ],
    )
    result_ht = result_ht.select(
        variant_id=result_ht.variant_id,
        chromosome=result_ht.contig,
        position=result_ht.position,
        maf=result_ht.maf,
        minor_allele_count=result_ht.minor_allele_count,
        call_rate=result_ht.call_rate,
        hwe_control_p=result_ht.hwe_control_p,
        hwe_case_p=result_ht.hwe_case_p,
        beta=result_ht.beta,
        se=result_ht.standard_error,
        z_stat=result_ht.z_stat,
        regression_p=result_ht.p_value,
        odds_ratio=hl.exp(result_ht.beta),
    )
    result_df = result_ht.to_pandas()
    if not result_df.empty:
        result_df["chromosome"] = result_df["chromosome"].astype(str).str.replace("chr", "", regex=False)
        for column in (
            "position",
            "regression_p",
            "beta",
            "se",
            "odds_ratio",
            "maf",
            "minor_allele_count",
            "call_rate",
            "hwe_control_p",
            "hwe_case_p",
        ):
            result_df[column] = pd.to_numeric(result_df[column], errors="coerce")
        result_df["n_samples"] = matrix_table_participants

    variant_qc_summary = pd.DataFrame(qc_rows)
    return _postprocess_pilot_results(
        result_df,
        variant_qc_summary,
        config,
        paths,
        label=label,
        chromosomes=chromosome_values,
        analysis_split=analysis_split,
        eligibility_flag=eligibility_flag,
        min_maf=min_maf,
        min_mac=min_mac,
        min_call_rate=min_call_rate,
        hwe_p_control=hwe_p_control,
        sample_counts=sample_counts,
        wgs_manifest_used=wgs_manifest_used,
        genotype_source=source,
        genotype_mt_path=genotype_mt_path,
        covariates_used=hail_covariates,
        dropped_covariates=dropped_covariates,
    )


__all__ = [
    "PILOT_COVARIATES",
    "hail_stage4_full_results_path",
    "hail_stage4_lead_hits_path",
    "hail_stage4_manhattan_path",
    "hail_stage4_qc_path",
    "hail_stage4_qq_path",
    "hail_stage4_report_path",
    "hail_pilot_lead_hits_path",
    "hail_pilot_manhattan_path",
    "hail_pilot_output_dir",
    "hail_pilot_qc_path",
    "hail_pilot_default_label",
    "hail_pilot_qq_path",
    "hail_pilot_report_path",
    "hail_pilot_results_path",
    "hail_pilot_variant_qc_summary_path",
    "run_stage4_hail_gwas",
    "run_stage4_hail_pilot_gwas",
]
