"""Helpers to prepare REGENIE inputs from the matched cohort."""

from __future__ import annotations

import pandas as pd

from .config import ProjectConfig
from .io_utils import ensure_parent_dir, write_dataframe, write_text
from .paths import ProjectPaths, join_path
from .sample_restriction import restrict_frame_for_gwas


def regenie_output_dir(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "stage4", "regenie")


def regenie_matched_samples_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "matched_samples.tsv")


def regenie_gwas_samples_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "matched_gwas_samples.tsv")


def regenie_keep_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "keep.tsv")


def regenie_keep_iid_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "keep_iid.txt")


def regenie_pheno_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "phenotypes.tsv")


def regenie_covar_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "covariates.tsv")


def regenie_numeric_covar_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "covariates_numeric.tsv")


def regenie_keep_bed_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "keep_bed.tsv")


def regenie_pheno_bed_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "phenotypes_bed.tsv")


def regenie_covar_bed_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "covariates_bed.tsv")


def regenie_keep_bed_complete_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "keep_bed_complete.tsv")


def regenie_pheno_bed_complete_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "phenotypes_bed_complete.tsv")


def regenie_covar_bed_complete_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "covariates_bed_complete.tsv")


def regenie_commands_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "run_regenie_acaf.sh")


def _restricted_matched_df(config: ProjectConfig, matched_df: pd.DataFrame) -> pd.DataFrame:
    return restrict_frame_for_gwas(config, matched_df, require_wgs=True).copy()


def _regenie_sample_frame(matched_df: pd.DataFrame) -> pd.DataFrame:
    sample = matched_df.copy()
    sample["IID"] = sample["person_id"].astype(str)
    sample["FID"] = sample["IID"]
    return sample


def _sample_manifest(sample_df: pd.DataFrame) -> pd.DataFrame:
    return sample_df[["FID", "IID"]].drop_duplicates().sort_values(["FID", "IID"]).reset_index(drop=True)


def _keep_file(sample_df: pd.DataFrame) -> pd.DataFrame:
    return _sample_manifest(sample_df)


def _keep_iid_file(sample_df: pd.DataFrame) -> pd.DataFrame:
    keep = _keep_file(sample_df)
    return keep[["IID"]].copy()


def _phenotype_file(sample_df: pd.DataFrame, outcome_column: str) -> pd.DataFrame:
    pheno = sample_df[["FID", "IID", outcome_column]].drop_duplicates(subset=["FID", "IID"]).copy()
    pheno[outcome_column] = pd.to_numeric(pheno[outcome_column], errors="coerce")
    return pheno.rename(columns={outcome_column: "rhabdo_case"})


def _covariate_columns(config: ProjectConfig) -> list[str]:
    stage = config.analysis.stage4
    if stage is None:
        return []
    covariates = [column for column in stage.covariates if column != config.analysis.matched_outcome_column]
    if "ancestry_pred" in config.cohort.exact_match_columns and "ancestry_pred" not in covariates:
        covariates.append("ancestry_pred")
    return covariates


def _covariate_file(sample_df: pd.DataFrame, covariates: list[str]) -> pd.DataFrame:
    columns = ["FID", "IID", *[column for column in covariates if column in sample_df.columns]]
    covar = sample_df[columns].drop_duplicates(subset=["FID", "IID"]).copy()
    return covar


def _numeric_covariate_file(sample_df: pd.DataFrame, covariates: list[str]) -> pd.DataFrame:
    numeric_covars = [column for column in covariates if column != "ancestry_pred" and column in sample_df.columns]
    columns = ["FID", "IID", *numeric_covars]
    covar = sample_df[columns].drop_duplicates(subset=["FID", "IID"]).copy()
    return covar


def _bed_id_frame(df: pd.DataFrame) -> pd.DataFrame:
    bed = df.copy()
    bed["FID"] = "0"
    bed["IID"] = bed["IID"].astype(str)
    columns = [column for column in df.columns if column not in {"FID", "IID"}]
    return bed[["FID", "IID", *columns]]


def _complete_case_frames(
    keep_bed: pd.DataFrame,
    pheno_bed: pd.DataFrame,
    covar_bed: pd.DataFrame,
    covariate_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = keep_bed.merge(pheno_bed, on=["FID", "IID"], how="inner")
    df = df.merge(covar_bed, on=["FID", "IID"], how="inner")

    numeric_columns = ["rhabdo_case", *[column for column in covariate_columns if column != "ancestry_pred"]]
    for column in numeric_columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    if "ancestry_pred" in df.columns:
        df["ancestry_pred"] = df["ancestry_pred"].replace("", pd.NA)

    needed = ["FID", "IID", "rhabdo_case", *[column for column in covariate_columns if column in df.columns]]
    complete = df.dropna(subset=needed).copy()

    keep_complete = complete[["FID", "IID"]].drop_duplicates().reset_index(drop=True)
    pheno_complete = complete[["FID", "IID", "rhabdo_case"]].drop_duplicates().reset_index(drop=True)
    covar_columns = ["FID", "IID", *[column for column in covariate_columns if column in complete.columns]]
    covar_complete = complete[covar_columns].drop_duplicates().reset_index(drop=True)
    return keep_complete, pheno_complete, covar_complete


def _command_template(paths: ProjectPaths, config: ProjectConfig) -> str:
    keep_path = regenie_keep_bed_complete_path(paths)
    pheno_path = regenie_pheno_bed_complete_path(paths)
    covar_path = regenie_covar_bed_complete_path(paths)
    covariate_columns = _covariate_columns(config)
    cat_covars = [column for column in covariate_columns if column == "ancestry_pred"]
    covar_list = ",".join(covariate_columns) if covariate_columns else ""
    cat_covar_list = ",".join(cat_covars)
    covar_line = f'  --covarColList "{covar_list}" \\\n' if covar_list else ""
    cat_covar_line = f'  --catCovarList "{cat_covar_list}" \\\n' if cat_covar_list else ""
    return f"""#!/usr/bin/env bash
set -euo pipefail

REGENIE_BIN="${{REGENIE_BIN:-regenie}}"
STEP1_BED_PREFIX="${{STEP1_BED_PREFIX:-/path/to/acaf_step1_pruned}}"
STEP2_BED_PREFIX_TEMPLATE="${{STEP2_BED_PREFIX_TEMPLATE:-/path/to/chr{{CHR}}_matched_biallelic_uid}}"
OUTDIR="{regenie_output_dir(paths)}"

mkdir -p "$OUTDIR"

"$REGENIE_BIN" \\
  --step 1 \\
  --bed "$STEP1_BED_PREFIX" \\
  --keep "{keep_path}" \\
  --covarFile "{covar_path}" \\
  --phenoFile "{pheno_path}" \\
  --phenoCol rhabdo_case \\
{covar_line}{cat_covar_line}  --bt \\
  --bsize 1000 \\
  --lowmem \\
  --lowmem-prefix "$OUTDIR/tmp" \\
  --out "$OUTDIR/step1"

for CHR in {{1..22}}; do
  STEP2_BED_PREFIX="${{STEP2_BED_PREFIX_TEMPLATE/\{{CHR\}}/${{CHR}}}}"
  "$REGENIE_BIN" \\
    --step 2 \\
    --bed "${{STEP2_BED_PREFIX}}" \\
    --chr "${{CHR}}" \\
    --keep "{keep_path}" \\
    --covarFile "{covar_path}" \\
    --phenoFile "{pheno_path}" \\
    --phenoCol rhabdo_case \\
{covar_line}{cat_covar_line}    --bt \\
    --firth --approx \\
    --pThresh 0.01 \\
    --minMAC {int(config.analysis.stage4.min_minor_allele_count)} \\
    --pred "$OUTDIR/step1_pred.list" \\
    --bsize 400 \\
    --out "$OUTDIR/step2_chr${{CHR}}"
done
"""


def prepare_regenie_inputs(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
) -> dict[str, str]:
    stage = config.analysis.stage4
    if stage is None:
        return {}

    restricted = _restricted_matched_df(config, matched_df)
    sample_df = _regenie_sample_frame(restricted)
    full_matched_manifest = _sample_manifest(_regenie_sample_frame(matched_df))
    keep = _keep_file(sample_df)
    keep_iid = _keep_iid_file(sample_df)
    pheno = _phenotype_file(sample_df, config.analysis.matched_outcome_column)
    covariate_columns = _covariate_columns(config)
    covar = _covariate_file(sample_df, covariate_columns)
    numeric_covar = _numeric_covariate_file(sample_df, covariate_columns)
    keep_bed = _bed_id_frame(keep)
    pheno_bed = _bed_id_frame(pheno)
    covar_bed = _bed_id_frame(covar)
    keep_bed_complete, pheno_bed_complete, covar_bed_complete = _complete_case_frames(
        keep_bed,
        pheno_bed,
        covar_bed,
        covariate_columns,
    )

    matched_samples_path = regenie_matched_samples_path(paths)
    gwas_samples_path = regenie_gwas_samples_path(paths)
    keep_path = regenie_keep_path(paths)
    keep_iid_path = regenie_keep_iid_path(paths)
    pheno_path = regenie_pheno_path(paths)
    covar_path = regenie_covar_path(paths)
    numeric_covar_path = regenie_numeric_covar_path(paths)
    keep_bed_path = regenie_keep_bed_path(paths)
    pheno_bed_path = regenie_pheno_bed_path(paths)
    covar_bed_path = regenie_covar_bed_path(paths)
    keep_bed_complete_path = regenie_keep_bed_complete_path(paths)
    pheno_bed_complete_path = regenie_pheno_bed_complete_path(paths)
    covar_bed_complete_path = regenie_covar_bed_complete_path(paths)
    commands_path = regenie_commands_path(paths)

    ensure_parent_dir(keep_path)
    write_dataframe(full_matched_manifest, matched_samples_path)
    write_dataframe(keep, gwas_samples_path)
    write_dataframe(keep, keep_path)
    write_dataframe(keep_iid, keep_iid_path)
    write_dataframe(pheno, pheno_path)
    write_dataframe(covar, covar_path)
    write_dataframe(numeric_covar, numeric_covar_path)
    write_dataframe(keep_bed, keep_bed_path)
    write_dataframe(pheno_bed, pheno_bed_path)
    write_dataframe(covar_bed, covar_bed_path)
    write_dataframe(keep_bed_complete, keep_bed_complete_path)
    write_dataframe(pheno_bed_complete, pheno_bed_complete_path)
    write_dataframe(covar_bed_complete, covar_bed_complete_path)
    write_text(_command_template(paths, config), commands_path)

    return {
        "matched_samples": matched_samples_path,
        "gwas_samples": gwas_samples_path,
        "keep": keep_path,
        "keep_iid": keep_iid_path,
        "phenotypes": pheno_path,
        "covariates": covar_path,
        "covariates_numeric": numeric_covar_path,
        "keep_bed": keep_bed_path,
        "phenotypes_bed": pheno_bed_path,
        "covariates_bed": covar_bed_path,
        "keep_bed_complete": keep_bed_complete_path,
        "phenotypes_bed_complete": pheno_bed_complete_path,
        "covariates_bed_complete": covar_bed_complete_path,
        "commands": commands_path,
    }


__all__ = [
    "prepare_regenie_inputs",
    "regenie_commands_path",
    "regenie_covar_path",
    "regenie_covar_bed_complete_path",
    "regenie_covar_bed_path",
    "regenie_gwas_samples_path",
    "regenie_keep_path",
    "regenie_keep_bed_complete_path",
    "regenie_keep_bed_path",
    "regenie_keep_iid_path",
    "regenie_matched_samples_path",
    "regenie_numeric_covar_path",
    "regenie_output_dir",
    "regenie_pheno_bed_complete_path",
    "regenie_pheno_bed_path",
    "regenie_pheno_path",
]
