#!/usr/bin/env bash

set -eu

# Broad rhabdomyolysis WGS-restricted clinical characterization plus a chr21
# ACAF/Hail GWAS smoke test in an All of Us/Verily Dataproc app.
#
# Deliberately avoids `set -o pipefail`; run_logged captures command status from
# PIPESTATUS so tee logging is readable without hiding the real command exit.

GIT_URL="${GIT_URL:-https://github.com/spand67/aou-workbench.git}"
REPO_DIR="${REPO_DIR:-$HOME/genomic-analyses/rhabdo-v1}"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"

COHORT_CONFIG="${COHORT_CONFIG:-configs/rhabdo/cohort.yaml}"
ANALYSIS_CONFIG="${ANALYSIS_CONFIG:-configs/rhabdo/analysis.yaml}"
GENOTYPE_SOURCE="${GENOTYPE_SOURCE:-acaf}"
CHROMOSOMES="${CHROMOSOMES:-21}"
GWAS_LABEL="${GWAS_LABEL:-acaf_chr21_broad_train_maf05_qc}"
GWAS_SLUG="${GWAS_SLUG:-${GWAS_LABEL//_/-}}"

MIN_MAF="${MIN_MAF:-0.05}"
MIN_MAC="${MIN_MAC:-20}"
MIN_CALL_RATE="${MIN_CALL_RATE:-0.98}"
HWE_P_CONTROL="${HWE_P_CONTROL:-1e-6}"
ANALYSIS_SPLIT="${ANALYSIS_SPLIT:-train}"
ELIGIBILITY_FLAG="${ELIGIBILITY_FLAG:-primary_model_eligible}"

clone_or_update_repo() {
  if [ ! -d "$REPO_DIR/.git" ]; then
    mkdir -p "$(dirname "$REPO_DIR")"
    git clone "$GIT_URL" "$REPO_DIR"
  fi
  cd "$REPO_DIR"
  git pull --ff-only origin main
}

install_package() {
  export PATH="$HOME/.local/bin:$PATH"
  export PYTHONPATH="$PWD/src:${PYTHONPATH:-}"
  python -m pip install --user --no-deps --force-reinstall -e .
}

run_logged() {
  local name="$1"
  shift
  mkdir -p outputs/logs
  local log_path="outputs/logs/${name}_${RUN_ID}.log"
  echo ""
  echo "===== $name ====="
  echo "Running: $*"
  echo "Running: $*" > "$log_path"
  "$@" 2>&1 | tee -a "$log_path"
  local status="${PIPESTATUS[0]}"
  if [ "$status" -ne 0 ]; then
    echo "Step failed: $name (exit $status). Log: $log_path" >&2
    exit "$status"
  fi
}

print_file_if_exists() {
  local title="$1"
  local path="$2"
  echo ""
  echo "===== $title ====="
  if [ -f "$path" ]; then
    cat "$path"
  else
    echo "Missing: $path"
  fi
}

print_head_if_exists() {
  local title="$1"
  local path="$2"
  local lines="${3:-40}"
  echo ""
  echo "===== $title ====="
  if [ -f "$path" ]; then
    head -n "$lines" "$path"
  else
    echo "Missing: $path"
  fi
}

main() {
  clone_or_update_repo
  install_package
  mkdir -p outputs/logs

  echo ""
  echo "===== run_config ====="
  echo "RUN_ID=$RUN_ID"
  echo "REPO_DIR=$REPO_DIR"
  echo "COHORT_CONFIG=$COHORT_CONFIG"
  echo "ANALYSIS_CONFIG=$ANALYSIS_CONFIG"
  echo "GENOTYPE_SOURCE=$GENOTYPE_SOURCE"
  echo "CHROMOSOMES=$CHROMOSOMES"
  echo "GWAS_LABEL=$GWAS_LABEL"
  echo "GWAS_SLUG=$GWAS_SLUG"
  echo "ANALYSIS_SPLIT=$ANALYSIS_SPLIT"
  echo "ELIGIBILITY_FLAG=$ELIGIBILITY_FLAG"

  run_logged broad_wgs_build \
    aou-workbench build-cohort \
      --require-wgs \
      --cohort-config "$COHORT_CONFIG" \
      --analysis-config "$ANALYSIS_CONFIG"

  run_logged broad_wgs_match \
    aou-workbench match-controls \
      --require-wgs \
      --cohort-config "$COHORT_CONFIG" \
      --analysis-config "$ANALYSIS_CONFIG"

  run_logged broad_wgs_characterize \
    aou-workbench characterize-cohort \
      --require-wgs \
      --cohort-config "$COHORT_CONFIG" \
      --analysis-config "$ANALYSIS_CONFIG"

  {
    print_file_if_exists "CONSORT" outputs/rhabdo-end-to-end/cohort/consort_counts.tsv
    print_file_if_exists "SPLITS" outputs/rhabdo-end-to-end/cohort/model_split_summary.tsv
    print_head_if_exists "TABLE 1 PREVIEW" outputs/rhabdo-end-to-end/cohort/table1_clinical_by_split.tsv 40
  } 2>&1 | tee "outputs/logs/broad_wgs_clinical_review_${RUN_ID}.log"

  run_logged hail_gwas_chr21_broad \
    aou-workbench run-hail-pilot-gwas \
      --cohort-config "$COHORT_CONFIG" \
      --analysis-config "$ANALYSIS_CONFIG" \
      --genotype-source "$GENOTYPE_SOURCE" \
      --chromosomes "$CHROMOSOMES" \
      --min-maf "$MIN_MAF" \
      --min-mac "$MIN_MAC" \
      --min-call-rate "$MIN_CALL_RATE" \
      --hwe-p-control "$HWE_P_CONTROL" \
      --analysis-split "$ANALYSIS_SPLIT" \
      --eligibility-flag "$ELIGIBILITY_FLAG" \
      --label "$GWAS_LABEL"

  {
    print_file_if_exists "CHR21 VARIANT QC" "outputs/rhabdo-end-to-end/stage4/hail_pilot/${GWAS_SLUG}/variant_qc_summary.tsv"
    print_file_if_exists "CHR21 GWAS QC" "outputs/rhabdo-end-to-end/stage4/hail_pilot/${GWAS_SLUG}/gwas_qc.json"
    print_head_if_exists "CHR21 LEAD HITS" "outputs/rhabdo-end-to-end/stage4/hail_pilot/${GWAS_SLUG}/lead_hits.tsv" 25
  } 2>&1 | tee "outputs/logs/hail_gwas_chr21_broad_review_${RUN_ID}.log"

  echo ""
  echo "Done."
  echo "Clinical dashboard: outputs/rhabdo-end-to-end/rhabdo_dashboard.html"
  echo "GWAS report: outputs/rhabdo-end-to-end/stage4/hail_pilot/${GWAS_SLUG}/report.md"
  echo "GWAS Manhattan: outputs/rhabdo-end-to-end/stage4/hail_pilot/${GWAS_SLUG}/manhattan.svg"
  echo "GWAS QQ: outputs/rhabdo-end-to-end/stage4/hail_pilot/${GWAS_SLUG}/qq.svg"
}

main "$@"
