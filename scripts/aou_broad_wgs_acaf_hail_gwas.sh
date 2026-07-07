#!/usr/bin/env bash

set -eu

# Broad rhabdomyolysis WGS-restricted clinical characterization plus an ACAF/Hail
# GWAS smoke test on chr5. Set RUN_AUTOSOMES=1 after reviewing chr5 to launch the
# full autosomal run. Deliberately avoids `set -o pipefail`; run_logged captures
# command status from PIPESTATUS so tee logging stays readable while failures
# still stop the script.

GIT_URL="${GIT_URL:-https://github.com/spand67/aou-workbench.git}"
REPO_DIR="${REPO_DIR:-$HOME/genomic-analyses/rhabdo-v1}"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"

RUN_CLINICAL="${RUN_CLINICAL:-1}"
RUN_SMOKE="${RUN_SMOKE:-1}"
RUN_AUTOSOMES="${RUN_AUTOSOMES:-0}"
AUTOSOMES_BACKGROUND="${AUTOSOMES_BACKGROUND:-1}"

COHORT_CONFIG="${COHORT_CONFIG:-configs/rhabdo/cohort.yaml}"
ANALYSIS_CONFIG="${ANALYSIS_CONFIG:-configs/rhabdo/analysis.yaml}"
GENOTYPE_SOURCE="${GENOTYPE_SOURCE:-acaf}"
SMOKE_CHROMOSOMES="${SMOKE_CHROMOSOMES:-5}"
AUTOSOME_CHROMOSOMES="${AUTOSOME_CHROMOSOMES:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"
SMOKE_LABEL="${SMOKE_LABEL:-acaf_chr5_broad_train_maf01_hwe_report}"
AUTOSOME_LABEL="${AUTOSOME_LABEL:-acaf_autosomes_broad_train_maf01_hwe_report}"

MIN_MAF="${MIN_MAF:-0.01}"
MIN_MAC="${MIN_MAC:-20}"
MIN_CALL_RATE="${MIN_CALL_RATE:-0.98}"
HWE_P_CONTROL="${HWE_P_CONTROL:-1e-6}"
HWE_FILTER_MODE="${HWE_FILTER_MODE:-report-only}"
ANALYSIS_SPLIT="${ANALYSIS_SPLIT:-train}"
ELIGIBILITY_FLAG="${ELIGIBILITY_FLAG:-primary_model_eligible}"
SMOKE_PREVIEW_N="${SMOKE_PREVIEW_N:-250000}"
AUTOSOME_PREVIEW_N="${AUTOSOME_PREVIEW_N:-0}"

slugify_label() {
  echo "$1" | tr '[:upper:]_' '[:lower:]-' | sed -E 's/[^a-z0-9-]+/-/g; s/^-+//; s/-+$//; s/-+/-/g'
}

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

review_clinical_outputs() {
  {
    print_file_if_exists "CONSORT" outputs/rhabdo-end-to-end/cohort/consort_counts.tsv
    print_file_if_exists "SPLITS" outputs/rhabdo-end-to-end/cohort/model_split_summary.tsv
    print_head_if_exists "TABLE 1 PREVIEW" outputs/rhabdo-end-to-end/cohort/table1_clinical_by_split.tsv 40
  } 2>&1 | tee "outputs/logs/broad_wgs_clinical_review_${RUN_ID}.log"
}

review_gwas_outputs() {
  local title="$1"
  local label="$2"
  local slug
  slug="$(slugify_label "$label")"
  local gwas_dir="outputs/rhabdo-end-to-end/stage4/hail_pilot/${slug}"
  {
    print_file_if_exists "${title} VARIANT QC" "${gwas_dir}/variant_qc_summary.tsv"
    print_file_if_exists "${title} GWAS QC" "${gwas_dir}/gwas_qc.json"
    print_head_if_exists "${title} LEAD HITS" "${gwas_dir}/lead_hits.tsv" 30
    print_head_if_exists "${title} RESULTS PREVIEW" "${gwas_dir}/gwas_results_preview.tsv" 20
    echo ""
    echo "===== ${title} PATHS ====="
    echo "Local GWAS dir: ${gwas_dir}"
    echo "Report: ${gwas_dir}/report.md"
    echo "Manhattan: ${gwas_dir}/manhattan.svg"
    echo "QQ: ${gwas_dir}/qq.svg"
  } 2>&1 | tee "outputs/logs/hail_gwas_${slug}_review_${RUN_ID}.log"
}

run_hail_gwas_command() {
  local chromosomes="$1"
  local label="$2"
  local preview_n="$3"
  aou-workbench run-hail-pilot-gwas \
    --cohort-config "$COHORT_CONFIG" \
    --analysis-config "$ANALYSIS_CONFIG" \
    --genotype-source "$GENOTYPE_SOURCE" \
    --chromosomes "$chromosomes" \
    --min-maf "$MIN_MAF" \
    --min-mac "$MIN_MAC" \
    --min-call-rate "$MIN_CALL_RATE" \
    --hwe-p-control "$HWE_P_CONTROL" \
    --hwe-filter-mode "$HWE_FILTER_MODE" \
    --analysis-split "$ANALYSIS_SPLIT" \
    --eligibility-flag "$ELIGIBILITY_FLAG" \
    --skip-variant-row-counts \
    --results-preview-n "$preview_n" \
    --label "$label"
}

write_autosomes_runner() {
  mkdir -p outputs/logs
  local runner="outputs/logs/hail_gwas_autosomes_broad_${RUN_ID}.sh"
  cat > "$runner" <<EOF
#!/usr/bin/env bash
set -eu
cd "$REPO_DIR"
export PATH="\$HOME/.local/bin:\$PATH"
export PYTHONPATH="\$PWD/src:\${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

aou-workbench run-hail-pilot-gwas \\
  --cohort-config "$COHORT_CONFIG" \\
  --analysis-config "$ANALYSIS_CONFIG" \\
  --genotype-source "$GENOTYPE_SOURCE" \\
  --chromosomes "$AUTOSOME_CHROMOSOMES" \\
  --min-maf "$MIN_MAF" \\
  --min-mac "$MIN_MAC" \\
  --min-call-rate "$MIN_CALL_RATE" \\
  --hwe-p-control "$HWE_P_CONTROL" \\
  --hwe-filter-mode "$HWE_FILTER_MODE" \\
  --analysis-split "$ANALYSIS_SPLIT" \\
  --eligibility-flag "$ELIGIBILITY_FLAG" \\
  --skip-variant-row-counts \\
  --results-preview-n "$AUTOSOME_PREVIEW_N" \\
  --label "$AUTOSOME_LABEL"
EOF
  chmod +x "$runner"
  echo "$runner"
}

launch_autosomes() {
  local log_path="outputs/logs/hail_gwas_autosomes_broad_${RUN_ID}.log"
  if [ "$AUTOSOMES_BACKGROUND" = "1" ]; then
    local runner
    runner="$(write_autosomes_runner)"
    echo ""
    echo "===== hail_gwas_autosomes_broad ====="
    echo "Launching background autosomal GWAS."
    echo "Runner: $runner"
    echo "Log: $log_path"
    nohup "$runner" > "$log_path" 2>&1 &
    echo "$!" > "outputs/logs/hail_gwas_autosomes_broad_${RUN_ID}.pid"
    echo "PID: $(cat "outputs/logs/hail_gwas_autosomes_broad_${RUN_ID}.pid")"
    echo "Monitor with: tail -f $REPO_DIR/$log_path"
  else
    run_logged hail_gwas_autosomes_broad \
      run_hail_gwas_command "$AUTOSOME_CHROMOSOMES" "$AUTOSOME_LABEL" "$AUTOSOME_PREVIEW_N"
    review_gwas_outputs "AUTOSOMES" "$AUTOSOME_LABEL"
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
  echo "RUN_CLINICAL=$RUN_CLINICAL"
  echo "RUN_SMOKE=$RUN_SMOKE"
  echo "RUN_AUTOSOMES=$RUN_AUTOSOMES"
  echo "AUTOSOMES_BACKGROUND=$AUTOSOMES_BACKGROUND"
  echo "COHORT_CONFIG=$COHORT_CONFIG"
  echo "ANALYSIS_CONFIG=$ANALYSIS_CONFIG"
  echo "GENOTYPE_SOURCE=$GENOTYPE_SOURCE"
  echo "SMOKE_CHROMOSOMES=$SMOKE_CHROMOSOMES"
  echo "AUTOSOME_CHROMOSOMES=$AUTOSOME_CHROMOSOMES"
  echo "SMOKE_LABEL=$SMOKE_LABEL"
  echo "AUTOSOME_LABEL=$AUTOSOME_LABEL"
  echo "MIN_MAF=$MIN_MAF"
  echo "MIN_MAC=$MIN_MAC"
  echo "MIN_CALL_RATE=$MIN_CALL_RATE"
  echo "HWE_P_CONTROL=$HWE_P_CONTROL"
  echo "HWE_FILTER_MODE=$HWE_FILTER_MODE"
  echo "ANALYSIS_SPLIT=$ANALYSIS_SPLIT"
  echo "ELIGIBILITY_FLAG=$ELIGIBILITY_FLAG"

  if [ "$RUN_CLINICAL" = "1" ]; then
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

    review_clinical_outputs
  fi

  if [ "$RUN_SMOKE" = "1" ]; then
    run_logged hail_gwas_chr5_broad \
      run_hail_gwas_command "$SMOKE_CHROMOSOMES" "$SMOKE_LABEL" "$SMOKE_PREVIEW_N"
    review_gwas_outputs "CHR5" "$SMOKE_LABEL"
  fi

  if [ "$RUN_AUTOSOMES" = "1" ]; then
    launch_autosomes
  else
    echo ""
    echo "Chr5 smoke/default run complete."
    echo "If the chr5 outputs look good, launch the full autosomal run with:"
    echo "RUN_CLINICAL=0 RUN_SMOKE=0 RUN_AUTOSOMES=1 bash scripts/aou_broad_wgs_acaf_hail_gwas.sh"
  fi
}

main "$@"
