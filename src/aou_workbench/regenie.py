"""Helpers to prepare REGENIE inputs from the matched cohort."""

from __future__ import annotations

import pandas as pd

from .config import ProjectConfig
from .io_utils import ensure_parent_dir, write_dataframe, write_text
from .paths import ProjectPaths, join_path


def regenie_output_dir(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "stage4", "regenie")


def regenie_keep_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "keep.tsv")


def regenie_pheno_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "phenotypes.tsv")


def regenie_covar_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "covariates.tsv")


def regenie_commands_path(paths: ProjectPaths) -> str:
    return join_path(regenie_output_dir(paths), "run_regenie_acaf.sh")


def _regenie_sample_frame(matched_df: pd.DataFrame) -> pd.DataFrame:
    sample = matched_df.copy()
    sample["IID"] = sample["person_id"].astype(str)
    sample["FID"] = sample["IID"]
    return sample


def _keep_file(sample_df: pd.DataFrame) -> pd.DataFrame:
    return sample_df[["FID", "IID"]].drop_duplicates().sort_values(["FID", "IID"]).reset_index(drop=True)


def _phenotype_file(sample_df: pd.DataFrame, outcome_column: str) -> pd.DataFrame:
    pheno = sample_df[["FID", "IID", outcome_column]].drop_duplicates(subset=["FID", "IID"]).copy()
    pheno[outcome_column] = pd.to_numeric(pheno[outcome_column], errors="coerce")
    return pheno.rename(columns={outcome_column: "rhabdo_case"})


def _covariate_file(sample_df: pd.DataFrame, covariates: tuple[str, ...]) -> pd.DataFrame:
    columns = ["FID", "IID", *[column for column in covariates if column in sample_df.columns]]
    if "ancestry_pred" in sample_df.columns and "ancestry_pred" not in columns:
        columns.append("ancestry_pred")
    covar = sample_df[columns].drop_duplicates(subset=["FID", "IID"]).copy()
    return covar


def _command_template(paths: ProjectPaths, config: ProjectConfig) -> str:
    keep_path = regenie_keep_path(paths)
    pheno_path = regenie_pheno_path(paths)
    covar_path = regenie_covar_path(paths)
    covariate_columns = [column for column in config.analysis.stage4.covariates if column != config.analysis.matched_outcome_column]
    if "ancestry_pred" in config.cohort.exact_match_columns and "ancestry_pred" not in covariate_columns:
        covariate_columns.append("ancestry_pred")
    cat_covars = []
    if "ancestry_pred" in covariate_columns:
        cat_covars.append("ancestry_pred")
    covar_list = ",".join(covariate_columns) if covariate_columns else ""
    cat_covar_list = ",".join(cat_covars)
    covar_line = f'  --covarColList "{covar_list}" \\\n' if covar_list else ""
    cat_covar_line = f'  --catCovarList "{cat_covar_list}" \\\n' if cat_covar_list else ""
    return f"""#!/usr/bin/env bash
set -euo pipefail

# Edit these to point at locally accessible or mounted ACAF PGEN prefixes.
# AoU documents the source data under:
# gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/pgen/
STEP1_PGEN_PREFIX="/path/to/acaf_step1_pruned_prefix"
STEP2_PGEN_PREFIX_TEMPLATE="/path/to/acaf_chr{{CHR}}"
OUTDIR="{regenie_output_dir(paths)}"

mkdir -p "$OUTDIR"

regenie \\
  --step 1 \\
  --pgen "$STEP1_PGEN_PREFIX" \\
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
  STEP2_PGEN_PREFIX="${{STEP2_PGEN_PREFIX_TEMPLATE/\{{CHR\}}/${{CHR}}}}"
  regenie \\
    --step 2 \\
    --pgen "${{STEP2_PGEN_PREFIX}}" \\
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

    sample_df = _regenie_sample_frame(matched_df)
    keep = _keep_file(sample_df)
    pheno = _phenotype_file(sample_df, config.analysis.matched_outcome_column)
    covar = _covariate_file(sample_df, stage.covariates)

    keep_path = regenie_keep_path(paths)
    pheno_path = regenie_pheno_path(paths)
    covar_path = regenie_covar_path(paths)
    commands_path = regenie_commands_path(paths)

    ensure_parent_dir(keep_path)
    write_dataframe(keep, keep_path)
    write_dataframe(pheno, pheno_path)
    write_dataframe(covar, covar_path)
    write_text(_command_template(paths, config), commands_path)

    return {
        "keep": keep_path,
        "phenotypes": pheno_path,
        "covariates": covar_path,
        "commands": commands_path,
    }


__all__ = [
    "prepare_regenie_inputs",
    "regenie_commands_path",
    "regenie_covar_path",
    "regenie_keep_path",
    "regenie_output_dir",
    "regenie_pheno_path",
]
