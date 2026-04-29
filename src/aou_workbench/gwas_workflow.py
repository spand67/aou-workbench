"""Terminal-first GWAS workspace generation for matched AoU analyses."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import ProjectConfig
from .io_utils import ensure_parent_dir, write_dataframe, write_text
from .paths import ProjectPaths, join_path
from .regenie import (
    prepare_regenie_inputs,
    regenie_commands_path,
    regenie_covar_path,
    regenie_covar_bed_complete_path,
    regenie_gwas_samples_path,
    regenie_keep_path,
    regenie_keep_bed_complete_path,
    regenie_keep_iid_path,
    regenie_matched_samples_path,
    regenie_output_dir,
    regenie_pheno_bed_complete_path,
    regenie_pheno_path,
)
from .sample_restriction import restrict_frame_for_gwas


def gwas_workspace_dir(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "stage4", "gwas_terminal")


def gwas_readme_path(paths: ProjectPaths) -> str:
    return join_path(gwas_workspace_dir(paths), "README.md")


def gwas_rewrite_script_path(paths: ProjectPaths) -> str:
    return join_path(gwas_workspace_dir(paths), "01_rewrite_acaf_to_bed.sh")


def gwas_prune_script_path(paths: ProjectPaths) -> str:
    return join_path(gwas_workspace_dir(paths), "02_build_step1_pruned_dataset.sh")


def gwas_step1_script_path(paths: ProjectPaths) -> str:
    return join_path(gwas_workspace_dir(paths), "03_run_regenie_step1.sh")


def gwas_step2_script_path(paths: ProjectPaths) -> str:
    return join_path(gwas_workspace_dir(paths), "04_run_regenie_step2.sh")


def gwas_merge_script_path(paths: ProjectPaths) -> str:
    return join_path(gwas_workspace_dir(paths), "05_merge_regenie_results.sh")


def gwas_sync_script_path(paths: ProjectPaths) -> str:
    return join_path(gwas_workspace_dir(paths), "06_sync_summary_stats.sh")


def gwas_dsub_step1_path(paths: ProjectPaths) -> str:
    return join_path(gwas_workspace_dir(paths), "submit_step1_dsub.sh")


def gwas_dsub_step2_path(paths: ProjectPaths) -> str:
    return join_path(gwas_workspace_dir(paths), "submit_step2_dsub.sh")


def _absolute(path: str) -> str:
    return str(Path(path).resolve()) if not path.startswith("gs://") else path


def _matched_sample_manifests(config: ProjectConfig, matched_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    full = matched_df.copy()
    full["IID"] = full["person_id"].astype(str)
    full["FID"] = full["IID"]
    full_manifest = full[["FID", "IID"]].drop_duplicates().sort_values(["FID", "IID"]).reset_index(drop=True)
    restricted = restrict_frame_for_gwas(config, matched_df, require_wgs=True).copy()
    restricted["IID"] = restricted["person_id"].astype(str)
    restricted["FID"] = restricted["IID"]
    restricted_manifest = (
        restricted[["FID", "IID"]].drop_duplicates().sort_values(["FID", "IID"]).reset_index(drop=True)
    )
    return full_manifest, restricted_manifest


def _covariate_columns(config: ProjectConfig) -> list[str]:
    stage = config.analysis.stage4
    if stage is None:
        return []
    covariates = [column for column in stage.covariates if column != config.analysis.matched_outcome_column]
    if "ancestry_pred" in config.cohort.exact_match_columns and "ancestry_pred" not in covariates:
        covariates.append("ancestry_pred")
    return covariates


def _bucket_output_root(config: ProjectConfig) -> str:
    bucket = config.workbench.workspace_bucket or "gs://YOUR_WORKSPACE_BUCKET"
    return join_path(bucket, "gwas", config.analysis.analysis_slug)


def _rewrite_script(config: ProjectConfig, paths: ProjectPaths) -> str:
    keep_iid_path = _absolute(regenie_keep_iid_path(paths))
    acaf_prefix = (
        config.workbench.acaf_mt_path.replace("/splitMT/hail.mt", "/pgen").replace("/multiMT/hail.mt", "/pgen")
    )
    return f"""#!/usr/bin/env bash
set -euo pipefail

PLINK2_BIN="${{PLINK2_BIN:-plink2}}"
GSUTIL_BIN="${{GSUTIL_BIN:-gsutil}}"
REQUESTER_PAYS_PROJECT="${{REQUESTER_PAYS_PROJECT:-${{GCS_REQUESTER_PAYS_PROJECT:-${{GOOGLE_PROJECT:-}}}}}}"
WORKDIR="${{WORKDIR:-$PWD/gwas_work}}"
DOWNLOAD_DIR="$WORKDIR/pgen"
BED_DIR="$WORKDIR/bed"
KEEP_FILE="{keep_iid_path}"
ACAF_PGEN_ROOT="{acaf_prefix}"

mkdir -p "$DOWNLOAD_DIR" "$BED_DIR"

if [ "$#" -gt 0 ]; then
  CHROMS=("$@")
elif [ -n "${{CHR:-}}" ]; then
  CHROMS=("$CHR")
else
  CHROMS=($(seq 1 22))
fi

for CHR in "${{CHROMS[@]}}"; do
  PREFIX="$DOWNLOAD_DIR/acaf_threshold.chr${{CHR}}"
  "$GSUTIL_BIN" -u "$REQUESTER_PAYS_PROJECT" cp \
    "$ACAF_PGEN_ROOT/acaf_threshold.chr${{CHR}}.pgen" \
    "$ACAF_PGEN_ROOT/acaf_threshold.chr${{CHR}}.pvar" \
    "$ACAF_PGEN_ROOT/acaf_threshold.chr${{CHR}}.psam" \
    "$DOWNLOAD_DIR/"
  "$PLINK2_BIN" \
    --pfile "$PREFIX" \
    --keep "$KEEP_FILE" \
    --maf {config.analysis.stage4.min_maf} \
    --mac {int(config.analysis.stage4.min_minor_allele_count)} \
    --max-alleles 2 \
    --threads "${{PLINK_THREADS:-8}}" \
    --memory "${{PLINK_MEMORY_MB:-24000}}" \
    --make-bed \
    --out "$BED_DIR/chr${{CHR}}_matched_biallelic"
  "$PLINK2_BIN" \
    --bfile "$BED_DIR/chr${{CHR}}_matched_biallelic" \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 1000 \
    --make-bed \
    --out "$BED_DIR/chr${{CHR}}_matched_biallelic_uid"
  rm -f "$PREFIX".pgen "$PREFIX".pvar "$PREFIX".psam
  rm -f \
    "$BED_DIR/chr${{CHR}}_matched_biallelic.bed" \
    "$BED_DIR/chr${{CHR}}_matched_biallelic.bim" \
    "$BED_DIR/chr${{CHR}}_matched_biallelic.fam" \
    "$BED_DIR/chr${{CHR}}_matched_biallelic.log"
done
"""


def _prune_script(paths: ProjectPaths) -> str:
    return f"""#!/usr/bin/env bash
set -euo pipefail

PLINK2_BIN="${{PLINK2_BIN:-plink2}}"
WORKDIR="${{WORKDIR:-$PWD/gwas_work}}"
BED_DIR="$WORKDIR/bed"
STEP1_DIR="$WORKDIR/step1"

mkdir -p "$STEP1_DIR"

if [ "$#" -gt 0 ]; then
  CHROMS=("$@")
elif [ -n "${{CHR:-}}" ]; then
  CHROMS=("$CHR")
else
  CHROMS=($(seq 1 22))
fi

for CHR in "${{CHROMS[@]}}"; do
  "$PLINK2_BIN" \
    --bfile "$BED_DIR/chr${{CHR}}_matched_biallelic_uid" \
    --indep-pairwise 200 50 0.2 \
    --out "$STEP1_DIR/prune_chr${{CHR}}"
  "$PLINK2_BIN" \
    --bfile "$BED_DIR/chr${{CHR}}_matched_biallelic_uid" \
    --extract "$STEP1_DIR/prune_chr${{CHR}}.prune.in" \
    --make-bed \
    --out "$STEP1_DIR/chr${{CHR}}_pruned"
done

if [ "${{#CHROMS[@]}}" -eq 1 ]; then
  CHR="${{CHROMS[0]}}"
  cp "$STEP1_DIR/chr${{CHR}}_pruned.bed" "$STEP1_DIR/acaf_step1_pruned.bed"
  cp "$STEP1_DIR/chr${{CHR}}_pruned.bim" "$STEP1_DIR/acaf_step1_pruned.bim"
  cp "$STEP1_DIR/chr${{CHR}}_pruned.fam" "$STEP1_DIR/acaf_step1_pruned.fam"
else
  ANCHOR_CHR="${{CHROMS[0]}}"
  : > "$STEP1_DIR/pruned_merge_list.txt"
  for CHR in "${{CHROMS[@]:1}}"; do
    echo "$STEP1_DIR/chr${{CHR}}_pruned" >> "$STEP1_DIR/pruned_merge_list.txt"
  done

  "$PLINK2_BIN" \
    --bfile "$STEP1_DIR/chr${{ANCHOR_CHR}}_pruned" \
    --pmerge-list "$STEP1_DIR/pruned_merge_list.txt" bfile \
    --make-bed \
    --out "$STEP1_DIR/acaf_step1_pruned"
fi

for CHR in "${{CHROMS[@]}}"; do
  rm -f \
    "$STEP1_DIR/chr${{CHR}}_pruned.bed" \
    "$STEP1_DIR/chr${{CHR}}_pruned.bim" \
    "$STEP1_DIR/chr${{CHR}}_pruned.fam" \
    "$STEP1_DIR/prune_chr${{CHR}}.prune.in" \
    "$STEP1_DIR/prune_chr${{CHR}}.prune.out"
done
"""


def _step1_script(config: ProjectConfig, paths: ProjectPaths) -> str:
    covariates = _covariate_columns(config)
    covar_list = ",".join(covariates)
    cat_covars = "ancestry_pred" if "ancestry_pred" in covariates else ""
    cat_line = f'  --catCovarList "{cat_covars}" \\\n' if cat_covars else ""
    return f"""#!/usr/bin/env bash
set -euo pipefail

REGENIE_BIN="${{REGENIE_BIN:-regenie}}"
WORKDIR="${{WORKDIR:-$PWD/gwas_work}}"
STEP1_PREFIX="${{STEP1_PREFIX:-$WORKDIR/step1/acaf_step1_pruned}}"
OUTDIR="$WORKDIR/regenie"

mkdir -p "$OUTDIR"

"$REGENIE_BIN" \\
  --step 1 \\
  --bed "$STEP1_PREFIX" \\
  --keep "{_absolute(regenie_keep_bed_complete_path(paths))}" \\
  --covarFile "{_absolute(regenie_covar_bed_complete_path(paths))}" \\
  --phenoFile "{_absolute(regenie_pheno_bed_complete_path(paths))}" \\
  --phenoCol rhabdo_case \\
  --covarColList "{covar_list}" \\
{cat_line}  --bt \\
  --bsize 1000 \\
  --lowmem \\
  --lowmem-prefix "$OUTDIR/tmp" \\
  --out "$OUTDIR/step1"
"""


def _step2_script(config: ProjectConfig, paths: ProjectPaths) -> str:
    covariates = _covariate_columns(config)
    covar_list = ",".join(covariates)
    cat_covars = "ancestry_pred" if "ancestry_pred" in covariates else ""
    cat_line = f'    --catCovarList "{cat_covars}" \\\n' if cat_covars else ""
    return f"""#!/usr/bin/env bash
set -euo pipefail

REGENIE_BIN="${{REGENIE_BIN:-regenie}}"
WORKDIR="${{WORKDIR:-$PWD/gwas_work}}"
BED_DIR="$WORKDIR/bed"
OUTDIR="$WORKDIR/regenie"

mkdir -p "$OUTDIR"

DELETE_CHROM_AFTER_STEP2="${{DELETE_CHROM_AFTER_STEP2:-1}}"

if [ "$#" -gt 0 ]; then
  CHROMS=("$@")
else
  CHROMS=($(seq 1 22))
fi

for CHR in "${{CHROMS[@]}}"; do
  "$REGENIE_BIN" \\
    --step 2 \\
    --bed "$BED_DIR/chr${{CHR}}_matched_biallelic_uid" \\
    --chr "${{CHR}}" \\
    --keep "{_absolute(regenie_keep_bed_complete_path(paths))}" \\
    --covarFile "{_absolute(regenie_covar_bed_complete_path(paths))}" \\
    --phenoFile "{_absolute(regenie_pheno_bed_complete_path(paths))}" \\
    --phenoCol rhabdo_case \\
    --covarColList "{covar_list}" \\
{cat_line}    --bt \\
    --firth --approx \\
    --pThresh 0.01 \\
    --minMAC {int(config.analysis.stage4.min_minor_allele_count)} \\
    --pred "$OUTDIR/step1_pred.list" \\
    --bsize 400 \\
    --out "$OUTDIR/step2_chr${{CHR}}"
  if [ "$DELETE_CHROM_AFTER_STEP2" = "1" ]; then
    rm -f \
      "$BED_DIR/chr${{CHR}}_matched_biallelic_uid.bed" \
      "$BED_DIR/chr${{CHR}}_matched_biallelic_uid.bim" \
      "$BED_DIR/chr${{CHR}}_matched_biallelic_uid.fam" \
      "$BED_DIR/chr${{CHR}}_matched_biallelic_uid.log"
  fi
done
"""


def _merge_script(paths: ProjectPaths) -> str:
    return f"""#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${{WORKDIR:-$PWD/gwas_work}}"
OUTDIR="$WORKDIR/regenie"
MERGED="$OUTDIR/step2_all_chr_rhabdo_case.regenie"

first_file=$(ls "$OUTDIR"/step2_chr*_rhabdo_case.regenie | head -n 1)
head -n 1 "$first_file" > "$MERGED"
for file in "$OUTDIR"/step2_chr*_rhabdo_case.regenie; do
  tail -n +2 "$file" >> "$MERGED"
done

echo "$MERGED"
"""


def _sync_script(config: ProjectConfig) -> str:
    bucket_root = _bucket_output_root(config)
    return f"""#!/usr/bin/env bash
set -euo pipefail

GSUTIL_BIN="${{GSUTIL_BIN:-gsutil}}"
WORKDIR="${{WORKDIR:-$PWD/gwas_work}}"
OUTDIR="$WORKDIR/regenie"
DEST="{bucket_root}/summary_stats"

"$GSUTIL_BIN" -m cp "$OUTDIR"/step2_chr*_rhabdo_case.regenie "$DEST"/
if [ -f "$OUTDIR/step2_all_chr_rhabdo_case.regenie" ]; then
  "$GSUTIL_BIN" cp "$OUTDIR/step2_all_chr_rhabdo_case.regenie" "$DEST/step2_all_chr_rhabdo_case.regenie"
fi
"""


def _dsub_step1_script(config: ProjectConfig, paths: ProjectPaths) -> str:
    bucket_root = _bucket_output_root(config)
    return f"""#!/usr/bin/env bash
set -euo pipefail

DSUB_BIN="${{DSUB_BIN:-dsub}}"
IMAGE="${{IMAGE:-gcr.io/YOUR_PROJECT/regenie-worker:latest}}"
PROJECT="${{PROJECT:-${{GOOGLE_PROJECT:-${{GCS_REQUESTER_PAYS_PROJECT:-YOUR_PROJECT}}}}}}"

"$DSUB_BIN" \\
  --provider google-cls-v2 \\
  --regions us-central1 \\
  --project "$PROJECT" \\
  --image "$IMAGE" \\
  --logging "{bucket_root}/dsub_logs/step1" \\
  --input KEEP_FILE="{_absolute(regenie_keep_bed_complete_path(paths))}" \\
  --input PHENO_FILE="{_absolute(regenie_pheno_bed_complete_path(paths))}" \\
  --input COVAR_FILE="{_absolute(regenie_covar_bed_complete_path(paths))}" \\
  --script "{_absolute(gwas_step1_script_path(paths))}" \\
  --output-recursive OUTDIR="{bucket_root}/step1"
"""


def _dsub_step2_script(config: ProjectConfig, paths: ProjectPaths) -> str:
    bucket_root = _bucket_output_root(config)
    return f"""#!/usr/bin/env bash
set -euo pipefail

DSUB_BIN="${{DSUB_BIN:-dsub}}"
IMAGE="${{IMAGE:-gcr.io/YOUR_PROJECT/regenie-worker:latest}}"
PROJECT="${{PROJECT:-${{GOOGLE_PROJECT:-${{GCS_REQUESTER_PAYS_PROJECT:-YOUR_PROJECT}}}}}}"

for CHR in $(seq 1 22); do
  "$DSUB_BIN" \\
    --provider google-cls-v2 \\
    --regions us-central1 \\
    --project "$PROJECT" \\
    --image "$IMAGE" \\
    --logging "{bucket_root}/dsub_logs/step2/chr${{CHR}}" \\
    --input KEEP_FILE="{_absolute(regenie_keep_bed_complete_path(paths))}" \\
    --input PHENO_FILE="{_absolute(regenie_pheno_bed_complete_path(paths))}" \\
    --input COVAR_FILE="{_absolute(regenie_covar_bed_complete_path(paths))}" \\
    --input STEP1_PRED="{bucket_root}/step1/step1_pred.list" \\
    --script "{_absolute(gwas_step2_script_path(paths))}" \\
    --env CHR="${{CHR}}" \\
    --output-recursive OUTDIR="{bucket_root}/step2/chr${{CHR}}"
done
"""


def _readme(config: ProjectConfig, paths: ProjectPaths) -> str:
    return f"""# Terminal-First Matched GWAS Workspace

This workspace is designed for a fresh Verily JupyterLab Spark cluster terminal.

## Included sample restrictions

- WGS-present participants from the Stage 1 sample manifest
- Max-unrelated subset {'enabled' if config.workbench.max_unrelated_path else 'not configured'}

## Key input files

- Matched sample manifest: `{regenie_matched_samples_path(paths)}`
- GWAS analysis sample manifest: `{regenie_gwas_samples_path(paths)}`
- REGENIE keep file: `{regenie_keep_path(paths)}`
- AoU IID-only keep file for PLINK2: `{regenie_keep_iid_path(paths)}`
- REGENIE phenotype file: `{regenie_pheno_path(paths)}`
- REGENIE covariate file: `{regenie_covar_path(paths)}`
- Complete-case BED-side keep file: `{regenie_keep_bed_complete_path(paths)}`
- Complete-case BED-side phenotype file: `{regenie_pheno_bed_complete_path(paths)}`
- Complete-case BED-side covariate file: `{regenie_covar_bed_complete_path(paths)}`

## Run order

1. Rewrite ACAF chromosomes to matched biallelic BED files:
   - `{gwas_rewrite_script_path(paths)}`
2. Build the genome-wide pruned Step 1 BED set:
   - `{gwas_prune_script_path(paths)}`
3. Run REGENIE Step 1:
   - `{gwas_step1_script_path(paths)}`
4. Run REGENIE Step 2 across chromosomes:
   - `{gwas_step2_script_path(paths)}`
5. Merge chromosome summary statistics:
   - `{gwas_merge_script_path(paths)}`
6. Copy retained summary statistics to the workspace bucket:
   - `{gwas_sync_script_path(paths)}`

## Notes

- Genotype handling is chromosome-scoped and ephemeral by design.
- The local terminal workflow uses PLINK2 to rewrite AoU ACAF PGEN inputs into matched, common, biallelic BED files with unique variant IDs before REGENIE.
- The rewrite step deletes localized `pgen/pvar/psam` files after each chromosome rewrite, and Step 2 deletes each chromosome BED after association testing completes unless `DELETE_CHROM_AFTER_STEP2=0`.
- The dsub templates target `us-central1` and write outputs under `{_bucket_output_root(config)}`.
- The legacy single-file REGENIE template remains at `{regenie_commands_path(paths)}`.
"""


def prepare_terminal_gwas_workspace(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
) -> dict[str, str]:
    outputs = prepare_regenie_inputs(config, matched_df, paths)
    workspace_dir = gwas_workspace_dir(paths)
    ensure_parent_dir(gwas_readme_path(paths))

    full_manifest, restricted_manifest = _matched_sample_manifests(config, matched_df)
    matched_manifest_path = join_path(workspace_dir, "matched_samples.tsv")
    restricted_manifest_path = join_path(workspace_dir, "matched_gwas_samples.tsv")
    write_dataframe(full_manifest, matched_manifest_path)
    write_dataframe(restricted_manifest, restricted_manifest_path)

    scripts = {
        gwas_rewrite_script_path(paths): _rewrite_script(config, paths),
        gwas_prune_script_path(paths): _prune_script(paths),
        gwas_step1_script_path(paths): _step1_script(config, paths),
        gwas_step2_script_path(paths): _step2_script(config, paths),
        gwas_merge_script_path(paths): _merge_script(paths),
        gwas_sync_script_path(paths): _sync_script(config),
        gwas_dsub_step1_path(paths): _dsub_step1_script(config, paths),
        gwas_dsub_step2_path(paths): _dsub_step2_script(config, paths),
        gwas_readme_path(paths): _readme(config, paths),
    }
    for path, contents in scripts.items():
        write_text(contents, path)

    outputs.update(
        {
            "gwas_workspace": workspace_dir,
            "matched_manifest": matched_manifest_path,
            "restricted_manifest": restricted_manifest_path,
            "rewrite_script": gwas_rewrite_script_path(paths),
            "prune_script": gwas_prune_script_path(paths),
            "step1_script": gwas_step1_script_path(paths),
            "step2_script": gwas_step2_script_path(paths),
            "merge_script": gwas_merge_script_path(paths),
            "sync_script": gwas_sync_script_path(paths),
            "dsub_step1": gwas_dsub_step1_path(paths),
            "dsub_step2": gwas_dsub_step2_path(paths),
            "readme": gwas_readme_path(paths),
        }
    )
    return outputs


__all__ = ["gwas_workspace_dir", "prepare_terminal_gwas_workspace"]
