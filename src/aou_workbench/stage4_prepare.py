"""Stage 4 preparation from the AoU ACAF smaller callset."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import ProjectConfig
from .hail_utils import ensure_hail
from .io_utils import ensure_parent_dir, write_dataframe, write_json
from .paths import ProjectPaths, join_path
from .preflight import apply_runtime_defaults


def stage4_plink_prefix(paths: ProjectPaths, chromosome: str) -> str:
    chrom = chromosome.replace("chr", "")
    return join_path(paths.run_root, "stage4", "plink", f"acaf_chr{chrom}")


def stage4_hail_result_path(paths: ProjectPaths, chromosome: str) -> str:
    chrom = chromosome.replace("chr", "")
    return join_path(paths.run_root, "stage4", f"hail_chr{chrom}_gwas.tsv")


def _normalized_chromosome_values(chromosome: str) -> set[str]:
    chrom = chromosome.strip()
    bare = chrom.replace("chr", "")
    return {f"chr{bare}"}


def _matched_sample_frame(matched_df: pd.DataFrame, config: ProjectConfig) -> pd.DataFrame:
    covariates = [column for column in config.analysis.stage4.covariates if column in matched_df.columns]
    columns = ["person_id", config.analysis.matched_outcome_column, *covariates]
    sample = matched_df[columns].drop_duplicates(subset=["person_id"]).copy()
    sample["person_id"] = sample["person_id"].astype(str)
    return sample


def prepare_stage4_acaf_subset(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
    *,
    chromosome: str = "chr19",
) -> dict[str, str]:
    config = apply_runtime_defaults(config)
    stage = config.analysis.stage4
    if stage is None:
        return {}

    sample_df = _matched_sample_frame(matched_df, config)
    sample_ids = sample_df["person_id"].tolist()
    chromosome_values = _normalized_chromosome_values(chromosome)
    genotype_path = stage.genotype_table
    annotation_path = stage.annotation_table or join_path(paths.run_root, "stage4", "gwas_annotations.tsv")
    plink_prefix = stage4_plink_prefix(paths, chromosome)
    pilot_result_path = stage4_hail_result_path(paths, chromosome)

    ensure_parent_dir(genotype_path)
    ensure_parent_dir(annotation_path)
    ensure_parent_dir(f"{plink_prefix}.bed")
    ensure_parent_dir(pilot_result_path)

    hl = ensure_hail(
        requester_pays_project=config.workbench.requester_pays_project,
        requester_pays_buckets=config.workbench.requester_pays_buckets,
    )
    print(f"Reading ACAF MT from {config.workbench.acaf_mt_path}", flush=True)
    mt = hl.read_matrix_table(config.workbench.acaf_mt_path)

    sample_set = hl.literal(set(sample_ids))
    chrom_set = hl.literal(chromosome_values)
    mt = mt.filter_cols(sample_set.contains(mt.s))
    mt = mt.filter_rows(chrom_set.contains(mt.locus.contig))
    print(f"Filtering ACAF MT to matched samples and {chromosome} before splitting multiallelics.", flush=True)
    mt = hl.split_multi_hts(mt)

    sample_ht = hl.Table.from_pandas(sample_df).key_by("person_id")
    mt = mt.annotate_cols(sample=sample_ht[mt.s])
    mt = mt.filter_cols(hl.is_defined(mt.sample))
    mt = mt.annotate_rows(
        variant_id=(
            hl.str(mt.locus.contig).replace("chr", "")
            + "-"
            + hl.str(mt.locus.position)
            + "-"
            + mt.alleles[0]
            + "-"
            + mt.alleles[1]
        ),
        maf=hl.if_else(
            hl.len(mt.info.AF) > 0,
            hl.float64(mt.info.AF[0]),
            hl.missing(hl.tfloat64),
        ),
    )
    mt = mt.filter_rows((hl.is_missing(mt.maf)) | (mt.maf >= stage.min_maf))

    row_count = mt.count_rows()
    col_count = mt.count_cols()
    print(
        f"Stage 4 ACAF subset contains {row_count} variant rows across {col_count} matched samples on {chromosome}.",
        flush=True,
    )

    # Export a local PLINK BED set for REGENIE or PLINK-based downstream analyses.
    hl.export_plink(
        mt,
        plink_prefix,
        ind_id=mt.s,
        fam_id=mt.s,
        is_female=(hl.float64(mt.sample.is_female) == 1.0),
        pheno=(hl.float64(mt.sample[config.analysis.matched_outcome_column]) == 1.0),
        varid=mt.variant_id,
    )

    entries = mt.filter_entries(hl.is_defined(mt.GT) & mt.GT.is_non_ref()).entries().key_by()
    entries = entries.select(
        person_id=entries.s,
        variant_id=entries.variant_id,
        chromosome=entries.locus.contig,
        position=entries.locus.position,
        gene=hl.str(""),
        dosage=hl.float64(entries.GT.n_alt_alleles()),
        maf=entries.maf,
    )
    entry_rows = entries.collect()
    genotype_df = pd.DataFrame(
        [
            {
                "person_id": row.person_id,
                "variant_id": row.variant_id,
                "chromosome": row.chromosome,
                "position": int(row.position),
                "gene": row.gene,
                "dosage": float(row.dosage),
                "maf": float(row.maf) if row.maf is not None else None,
            }
            for row in entry_rows
        ]
    )
    annotation_df = (
        genotype_df[["variant_id", "gene"]]
        .drop_duplicates()
        .rename(columns={"gene": "nearest_gene"})
        .reset_index(drop=True)
    )
    write_dataframe(genotype_df, genotype_path)
    write_dataframe(annotation_df, annotation_path)

    result_ht = hl.logistic_regression_rows(
        test="wald",
        y=hl.float64(mt.sample[config.analysis.matched_outcome_column]),
        x=mt.GT.n_alt_alleles(),
        covariates=[
            1.0,
            *[
                hl.float64(mt.sample[column])
                for column in stage.covariates
                if column in sample_df.columns and column != config.analysis.matched_outcome_column
            ],
        ],
        pass_through=[
            mt.variant_id,
            mt.locus.contig,
            mt.locus.position,
            mt.maf,
        ],
    )
    result_ht = result_ht.select(
        variant_id=result_ht.variant_id,
        chromosome=result_ht.contig,
        position=result_ht.position,
        maf=result_ht.maf,
        beta=result_ht.beta,
        standard_error=result_ht.standard_error,
        z_stat=result_ht.z_stat,
        p_value=result_ht.p_value,
    )
    result_ht.export(pilot_result_path)

    meta_path = str(Path(genotype_path).with_suffix(Path(genotype_path).suffix + ".meta.json"))
    write_json(
        {
            "chromosome": chromosome,
            "matched_samples_requested": int(len(sample_ids)),
            "matched_samples_in_subset": int(col_count),
            "variant_rows_in_subset": int(row_count),
            "nonref_entries_written": int(len(genotype_df)),
            "plink_prefix": plink_prefix,
            "pilot_result_path": pilot_result_path,
            "acaf_mt_path": config.workbench.acaf_mt_path,
        },
        meta_path,
    )
    return {
        "genotypes": genotype_path,
        "annotations": annotation_path,
        "plink_prefix": plink_prefix,
        "hail_gwas": pilot_result_path,
    }


__all__ = [
    "prepare_stage4_acaf_subset",
    "stage4_hail_result_path",
    "stage4_plink_prefix",
]
