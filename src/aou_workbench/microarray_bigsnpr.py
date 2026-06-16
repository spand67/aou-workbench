"""Sparse microarray prediction using PLINK2-prepared data and bigsnpr."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess

import pandas as pd

from .config import ProjectConfig
from .io_utils import ensure_parent_dir, slugify, write_dataframe, write_json
from .microarray_plink_gwas import (
    PLINK_COVARIATES,
    _bim_row_counts,
    _chromosome_arg,
    _count_snplist,
    _expand_local_prefix,
    _plink_common_options,
    _plink_files_exist,
    _run_command,
    read_plink_fam,
    stage_microarray_plink_files,
)
from .paths import ProjectPaths, join_path
from .reporting import write_stage_report
from .stage4_hail_gwas import _hail_pilot_sample_frame, _normalize_autosomal_chromosomes


BIGSNPR_DEFAULT_ALPHAS = (1.0, 0.5, 0.1, 0.01)


def microarray_bigsnpr_default_label(chromosomes: list[str], *, min_maf: float = 0.05) -> str:
    chrom_label = "-".join(chromosomes)
    maf_text = f"{min_maf:g}"
    maf_label = maf_text[2:] if maf_text.startswith("0.") else maf_text.replace(".", "")
    return f"microarray_bigsnpr_chr{chrom_label}_maf{maf_label}_train_qc"


def microarray_bigsnpr_output_dir(paths: ProjectPaths, label: str) -> str:
    return join_path(paths.run_root, "stage4", "microarray_bigsnpr", slugify(label))


def microarray_bigsnpr_metadata_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "sample_metadata.tsv")


def microarray_bigsnpr_train_keep_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "train_keep.tsv")


def microarray_bigsnpr_test_keep_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "test_keep.tsv")


def microarray_bigsnpr_train_test_keep_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "train_test_keep.tsv")


def microarray_bigsnpr_train_controls_keep_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "train_controls_keep.tsv")


def microarray_bigsnpr_variant_qc_summary_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "variant_qc_summary.tsv")


def microarray_bigsnpr_qc_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "bigsnpr_qc.json")


def microarray_bigsnpr_script_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "run_bigsnpr_model.R")


def microarray_bigsnpr_report_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "report.md")


def microarray_bigsnpr_metrics_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "metrics.tsv")


def microarray_bigsnpr_predictions_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "predictions.tsv")


def microarray_bigsnpr_selected_variants_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "selected_variants.tsv")


def microarray_bigsnpr_score_distribution_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "score_distribution.svg")


def microarray_bigsnpr_risk_decile_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), "risk_decile.svg")


def _output_prefix(paths: ProjectPaths, label: str, name: str) -> str:
    return join_path(microarray_bigsnpr_output_dir(paths, label), name)


def _write_space_delimited(frame: pd.DataFrame, path: str, *, header: bool = False) -> None:
    ensure_parent_dir(path)
    frame.to_csv(path, sep=" ", index=False, header=header)


def _flag_value(frame: pd.DataFrame, column: str) -> pd.Series:
    values = frame[column]
    numeric = pd.to_numeric(values, errors="coerce").fillna(0)
    string = values.astype("string").str.lower().str.strip()
    return numeric.ne(0) | string.isin({"true", "t", "yes", "y"})


def _sample_frame_for_split(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    fam_ids: set[str],
    *,
    split: str,
    eligibility_flag: str,
) -> tuple[pd.DataFrame, list[str], list[str], dict[str, int]]:
    sample, kept_covariates, _, dropped_covariates, counts, _ = _hail_pilot_sample_frame(
        matched_df,
        config,
        analysis_split=split,
        eligibility_flag=eligibility_flag,
        covariates=PLINK_COVARIATES,
        restrict_to_wgs_manifest=False,
    )
    outcome = config.analysis.matched_outcome_column
    sample = sample[sample["person_id"].astype(str).isin(fam_ids)].copy()
    sample["FID"] = "0"
    sample["IID"] = sample["person_id"].astype(str)
    sample["analysis_split"] = split
    sample[outcome] = pd.to_numeric(sample[outcome], errors="coerce").astype(int)
    if "train_cv_fold" not in sample.columns:
        source = matched_df[["person_id", "train_cv_fold"]].drop_duplicates("person_id") if "train_cv_fold" in matched_df.columns else None
        if source is not None:
            sample = sample.merge(source, on="person_id", how="left")
        else:
            sample["train_cv_fold"] = ""
    counts[f"after_microarray_fam_overlap_{split}_participants"] = int(sample["IID"].nunique())
    counts[f"after_microarray_fam_overlap_{split}_cases"] = int(sample[outcome].eq(1).sum())
    counts[f"after_microarray_fam_overlap_{split}_controls"] = int(sample[outcome].eq(0).sum())
    return sample, kept_covariates, dropped_covariates, counts


def write_microarray_bigsnpr_inputs(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    fam_df: pd.DataFrame,
    paths: ProjectPaths,
    *,
    label: str,
    eligibility_flag: str = "primary_model_eligible",
) -> tuple[pd.DataFrame, list[str], list[str], dict[str, int]]:
    """Write train/test sample metadata and PLINK keep files for bigsnpr modeling."""

    fam_ids = set(fam_df["IID"].astype(str))
    train, train_covariates, train_dropped, train_counts = _sample_frame_for_split(
        config,
        matched_df,
        fam_ids,
        split="train",
        eligibility_flag=eligibility_flag,
    )
    test, test_covariates, test_dropped, test_counts = _sample_frame_for_split(
        config,
        matched_df,
        fam_ids,
        split="test",
        eligibility_flag=eligibility_flag,
    )
    covariates = [column for column in train_covariates if column in test_covariates]
    dropped_covariates = sorted(set(train_dropped).union(test_dropped).union(set(train_covariates).symmetric_difference(test_covariates)))
    outcome = config.analysis.matched_outcome_column
    metadata_columns = ["FID", "IID", "person_id", "analysis_split", "train_cv_fold", outcome, *covariates]
    metadata = pd.concat([train[metadata_columns], test[metadata_columns]], ignore_index=True)

    train_keep = train[["FID", "IID"]].drop_duplicates()
    test_keep = test[["FID", "IID"]].drop_duplicates()
    train_test_keep = metadata[["FID", "IID"]].drop_duplicates()
    train_controls_keep = train.loc[pd.to_numeric(train[outcome], errors="coerce").eq(0), ["FID", "IID"]].drop_duplicates()

    _write_space_delimited(train_keep, microarray_bigsnpr_train_keep_path(paths, label))
    _write_space_delimited(test_keep, microarray_bigsnpr_test_keep_path(paths, label))
    _write_space_delimited(train_test_keep, microarray_bigsnpr_train_test_keep_path(paths, label))
    _write_space_delimited(train_controls_keep, microarray_bigsnpr_train_controls_keep_path(paths, label))
    write_dataframe(metadata, microarray_bigsnpr_metadata_path(paths, label))

    counts: dict[str, int] = {}
    counts.update(train_counts)
    counts.update(test_counts)
    counts["train_test_microarray_fam_overlap_participants"] = int(metadata["IID"].nunique())
    counts["train_test_microarray_fam_overlap_cases"] = int(pd.to_numeric(metadata[outcome], errors="coerce").eq(1).sum())
    counts["train_test_microarray_fam_overlap_controls"] = int(pd.to_numeric(metadata[outcome], errors="coerce").eq(0).sum())
    return metadata, covariates, dropped_covariates, counts


def build_microarray_bigsnpr_plink_commands(
    *,
    plink2_bin: str,
    plink_prefix: str,
    paths: ProjectPaths,
    label: str,
    chromosomes: list[str],
    min_maf: float,
    min_mac: int,
    min_call_rate: float,
    hwe_p_control: float,
    threads: int | None = None,
    memory_mb: int | None = None,
) -> tuple[list[str], list[str], list[str]]:
    common = _plink_common_options(threads, memory_mb)
    chr_arg = _chromosome_arg(chromosomes)
    max_missing = max(0.0, min(1.0, 1.0 - float(min_call_rate)))
    analysis_qc_prefix = _output_prefix(paths, label, "train_analysis_qc")
    control_hwe_prefix = _output_prefix(paths, label, "train_control_hwe_qc")
    subset_prefix = _output_prefix(paths, label, "bigsnpr_array_subset")
    analysis_qc_cmd = [
        plink2_bin,
        "--bfile",
        plink_prefix,
        "--chr",
        chr_arg,
        "--snps-only",
        "just-acgt",
        "--keep",
        microarray_bigsnpr_train_keep_path(paths, label),
        "--maf",
        str(float(min_maf)),
        "--mac",
        str(int(min_mac)),
        "--geno",
        f"{max_missing:.12g}",
        "--write-snplist",
        "--out",
        analysis_qc_prefix,
        *common,
    ]
    control_hwe_cmd = [
        plink2_bin,
        "--bfile",
        plink_prefix,
        "--chr",
        chr_arg,
        "--snps-only",
        "just-acgt",
        "--keep",
        microarray_bigsnpr_train_controls_keep_path(paths, label),
        "--extract",
        f"{analysis_qc_prefix}.snplist",
        "--hwe",
        str(float(hwe_p_control)),
        "keep-fewhet",
        "midp",
        "--write-snplist",
        "--out",
        control_hwe_prefix,
        *common,
    ]
    make_bed_cmd = [
        plink2_bin,
        "--bfile",
        plink_prefix,
        "--chr",
        chr_arg,
        "--keep",
        microarray_bigsnpr_train_test_keep_path(paths, label),
        "--extract",
        f"{control_hwe_prefix}.snplist",
        "--make-bed",
        "--out",
        subset_prefix,
        *common,
    ]
    return analysis_qc_cmd, control_hwe_cmd, make_bed_cmd


def _subset_prefix(paths: ProjectPaths, label: str) -> str:
    return _output_prefix(paths, label, "bigsnpr_array_subset")


def _subset_files_exist(paths: ProjectPaths, label: str) -> bool:
    prefix = _subset_prefix(paths, label)
    return all(Path(f"{prefix}.{suffix}").exists() for suffix in ("bed", "bim", "fam"))


def _variant_qc_summary(
    *,
    plink_prefix: str,
    paths: ProjectPaths,
    label: str,
    chromosomes: list[str],
    min_maf: float,
    min_mac: int,
    min_call_rate: float,
    hwe_p_control: float,
) -> pd.DataFrame:
    total_bim_rows, chromosome_rows = _bim_row_counts(plink_prefix, chromosomes)
    analysis_qc_rows = _count_snplist(f"{_output_prefix(paths, label, 'train_analysis_qc')}.snplist")
    hwe_rows = _count_snplist(f"{_output_prefix(paths, label, 'train_control_hwe_qc')}.snplist")
    subset_bim = Path(f"{_subset_prefix(paths, label)}.bim")
    subset_rows = 0
    if subset_bim.exists():
        with subset_bim.open("r", encoding="utf-8") as handle:
            subset_rows = sum(1 for line in handle if line.strip())
    return pd.DataFrame(
        [
            {
                "filter": "chromosome_restriction",
                "threshold": f"chr {','.join(chromosomes)}",
                "rows_before": total_bim_rows,
                "rows_after": chromosome_rows,
                "rows_removed": max(total_bim_rows - chromosome_rows, 0),
            },
            {
                "filter": "train_sample_variant_qc",
                "threshold": f"biallelic A/C/G/T SNP, MAF >= {min_maf}, MAC >= {min_mac}, call rate >= {min_call_rate}",
                "rows_before": chromosome_rows,
                "rows_after": analysis_qc_rows,
                "rows_removed": max(chromosome_rows - analysis_qc_rows, 0),
            },
            {
                "filter": "train_control_hwe_p",
                "threshold": f">= {hwe_p_control}",
                "rows_before": analysis_qc_rows,
                "rows_after": hwe_rows,
                "rows_removed": max(analysis_qc_rows - hwe_rows, 0),
            },
            {
                "filter": "train_test_plink_subset",
                "threshold": "BED/BIM/FAM subset written for bigsnpr",
                "rows_before": hwe_rows,
                "rows_after": subset_rows,
                "rows_removed": max(hwe_rows - subset_rows, 0),
            },
        ]
    )


def _format_r_vector(values: tuple[float, ...] | list[float]) -> str:
    return "c(" + ", ".join(f"{float(value):.12g}" for value in values) + ")"


def write_bigsnpr_r_script(
    paths: ProjectPaths,
    *,
    label: str,
    outcome_column: str,
    covariates: list[str],
    alphas: tuple[float, ...],
    folds: int,
    nlambda: int,
    dfmax: int,
    ncores: int,
) -> str:
    """Write the R script that imports the PLINK subset and runs big_spLogReg."""

    output_dir = Path(microarray_bigsnpr_output_dir(paths, label))
    subset_prefix = _subset_prefix(paths, label)
    script_path = microarray_bigsnpr_script_path(paths, label)
    covar_r = "c(" + ", ".join(json.dumps(column) for column in covariates) + ")"
    script = f"""#!/usr/bin/env Rscript
# Keep BLAS/OpenMP single-threaded; `big_spLogReg(ncores=...)` provides the
# outer parallelism. This avoids bigstatsr's nested-parallelism guard.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

required <- c("bigsnpr", "bigstatsr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {{
  stop("Missing required R package(s): ", paste(missing, collapse = ", "),
       ". Install in the AoU app with: install.packages(c('bigsnpr','bigstatsr'), repos='https://cloud.r-project.org')")
}}

out_dir <- {json.dumps(str(output_dir))}
subset_prefix <- {json.dumps(subset_prefix)}
metadata_path <- {json.dumps(microarray_bigsnpr_metadata_path(paths, label))}
outcome_column <- {json.dumps(outcome_column)}
covariate_columns <- {covar_r}
backingfile <- file.path(out_dir, "bigsnpr_array_subset")
set.seed(20260616)

if (!file.exists(paste0(backingfile, ".rds"))) {{
  bigsnpr::snp_readBed(paste0(subset_prefix, ".bed"), backingfile = backingfile)
}}
obj <- bigsnpr::snp_attach(paste0(backingfile, ".rds"))
Gna <- obj$genotypes
G <- bigsnpr::snp_fastImputeSimple(Gna, method = "mean2", ncores = {int(ncores)})
fam <- obj$fam
map <- obj$map
metadata <- read.delim(metadata_path, sep = "\\t", stringsAsFactors = FALSE, check.names = FALSE)
iid_col <- if ("sample.ID" %in% names(fam)) "sample.ID" else if ("sample_id" %in% names(fam)) "sample_id" else names(fam)[2]
metadata$row_index <- match(as.character(metadata$IID), as.character(fam[[iid_col]]))
metadata <- metadata[!is.na(metadata$row_index), , drop = FALSE]
metadata$row_index <- as.integer(metadata$row_index)
metadata[[outcome_column]] <- as.integer(metadata[[outcome_column]])

train_meta <- metadata[metadata$analysis_split == "train", , drop = FALSE]
test_meta <- metadata[metadata$analysis_split == "test", , drop = FALSE]
if (nrow(train_meta) == 0 || nrow(test_meta) == 0) stop("Need non-empty train and test rows in sample metadata.")

train_ind <- train_meta$row_index
test_ind <- test_meta$row_index
y_train <- as.integer(train_meta[[outcome_column]])
y_test <- as.integer(test_meta[[outcome_column]])
covar_train <- as.matrix(train_meta[, covariate_columns, drop = FALSE])
covar_test <- as.matrix(test_meta[, covariate_columns, drop = FALSE])
storage.mode(covar_train) <- "double"
storage.mode(covar_test) <- "double"

ind_sets <- NULL
if ("train_cv_fold" %in% names(train_meta)) {{
  fold_text <- as.character(train_meta$train_cv_fold)
  fold_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", fold_text)))
  if (all(!is.na(fold_num)) && length(unique(fold_num)) >= 2) ind_sets <- fold_num
}}

fit <- bigstatsr::big_spLogReg(
  G,
  y01.train = y_train,
  ind.train = train_ind,
  ind.col = seq_len(ncol(G)),
  covar.train = covar_train,
  pf.covar = rep(0, ncol(covar_train)),
  alphas = {_format_r_vector(alphas)},
  K = {int(folds)},
  ind.sets = ind_sets,
  nlambda = {int(nlambda)},
  dfmax = {int(dfmax)},
  warn = TRUE,
  ncores = {int(ncores)}
)

train_score <- as.numeric(predict(fit, G, ind.row = train_ind, covar.row = covar_train))
test_score <- as.numeric(predict(fit, G, ind.row = test_ind, covar.row = covar_test))

auc_rank <- function(y, score) {{
  y <- as.integer(y)
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(score, ties.method = "average")
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}}
average_precision <- function(y, score) {{
  y <- as.integer(y)
  if (sum(y == 1) == 0) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  y_ord <- y[ord]
  precision <- cumsum(y_ord == 1) / seq_along(y_ord)
  sum(precision[y_ord == 1]) / sum(y_ord == 1)
}}
metric_row <- function(name, y, score) {{
  data.frame(
    evaluation_set = name,
    n = length(y),
    cases = sum(y == 1),
    controls = sum(y == 0),
    prevalence = mean(y == 1),
    roc_auc = auc_rank(y, score),
    average_precision = average_precision(y, score),
    mean_score_case = mean(score[y == 1]),
    mean_score_control = mean(score[y == 0]),
    sd_score_case = stats::sd(score[y == 1]),
    sd_score_control = stats::sd(score[y == 0])
  )
}}
metrics <- rbind(metric_row("train", y_train, train_score), metric_row("test", y_test, test_score))
write.table(metrics, file.path(out_dir, "metrics.tsv"), sep = "\\t", row.names = FALSE, quote = FALSE)

predictions <- rbind(
  data.frame(FID = train_meta$FID, IID = train_meta$IID, analysis_split = "train", outcome = y_train, bigsnpr_score = train_score),
  data.frame(FID = test_meta$FID, IID = test_meta$IID, analysis_split = "test", outcome = y_test, bigsnpr_score = test_score)
)
write.table(predictions, file.path(out_dir, "predictions.tsv"), sep = "\\t", row.names = FALSE, quote = FALSE)

summary_fit <- summary(fit, best.only = TRUE)
beta <- unlist(summary_fit$beta[[1]])
summary_out <- summary_fit
for (nm in names(summary_out)) {{
  if (is.list(summary_out[[nm]])) {{
    summary_out[[nm]] <- vapply(summary_out[[nm]], function(x) paste(as.character(x), collapse = ";"), character(1))
  }}
}}
write.table(summary_out, file.path(out_dir, "fit_summary.tsv"), sep = "\\t", row.names = FALSE, quote = FALSE)
snp_beta <- beta[seq_len(min(length(beta), ncol(G)))]
selected_idx <- which(abs(snp_beta) > 0)
if (length(selected_idx) > 0) {{
  marker_col <- if ("marker.ID" %in% names(map)) "marker.ID" else names(map)[2]
  pos_col <- if ("physical.pos" %in% names(map)) "physical.pos" else if ("physical_pos" %in% names(map)) "physical_pos" else names(map)[4]
  chr_col <- if ("chromosome" %in% names(map)) "chromosome" else names(map)[1]
  selected <- data.frame(
    variant_index = selected_idx,
    variant_id = as.character(map[[marker_col]][selected_idx]),
    chromosome = as.character(map[[chr_col]][selected_idx]),
    position = map[[pos_col]][selected_idx],
    beta = snp_beta[selected_idx]
  )
}} else {{
  selected <- data.frame(variant_index = integer(), variant_id = character(), chromosome = character(), position = numeric(), beta = numeric())
}}
write.table(selected, file.path(out_dir, "selected_variants.tsv"), sep = "\\t", row.names = FALSE, quote = FALSE)

svg(file.path(out_dir, "score_distribution.svg"), width = 7, height = 4)
hist(test_score[y_test == 0], breaks = 30, col = rgb(0.2, 0.45, 0.75, 0.45), border = NA,
     main = "bigsnpr Score Distribution: Held-out Test", xlab = "Sparse SNP model score")
hist(test_score[y_test == 1], breaks = 30, col = rgb(0.85, 0.25, 0.2, 0.45), border = NA, add = TRUE)
legend("topright", legend = c("Controls", "Cases"), fill = c(rgb(0.2, 0.45, 0.75, 0.45), rgb(0.85, 0.25, 0.2, 0.45)), bty = "n")
dev.off()

svg(file.path(out_dir, "risk_decile.svg"), width = 7, height = 4)
breaks <- unique(quantile(test_score, probs = seq(0, 1, 0.1), na.rm = TRUE))
if (length(breaks) >= 2) {{
  decile <- cut(test_score, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  decile_df <- aggregate(y_test, by = list(decile = decile), FUN = mean)
  names(decile_df)[2] <- "case_rate"
  barplot(decile_df$case_rate, names.arg = decile_df$decile, col = "#496A81",
          main = "Held-out Test Case Rate by bigsnpr Score Decile", xlab = "Score decile", ylab = "Case rate")
}} else {{
  plot.new()
  title("Held-out Test Case Rate by bigsnpr Score Decile")
  text(0.5, 0.5, "Score deciles unavailable")
}}
dev.off()

report <- c(
  "# Microarray bigsnpr Sparse Model",
  "",
  paste0("- Label: `{label}`"),
  "- Data source: AoU v8 microarray PLINK subset after train-sample variant QC.",
  "- Model: `bigstatsr::big_spLogReg` sparse logistic regression with age, sex, and PC covariates included unpenalized.",
  paste0("- Train rows: ", nrow(train_meta), "; test rows: ", nrow(test_meta)),
  "- Missing genotypes were imputed with `bigsnpr::snp_fastImputeSimple(method = 'mean2')` before sparse logistic regression.",
  paste0("- Selected SNP coefficients: ", nrow(selected)),
  "",
  "## Metrics",
  "",
  paste(capture.output(print(metrics)), collapse = "\\n")
)
writeLines(report, file.path(out_dir, "report.md"))
"""
    ensure_parent_dir(script_path)
    Path(script_path).write_text(script, encoding="utf-8")
    return script_path


def _run_rscript_command(rscript_bin: str, script: str) -> None:
    env = os.environ.copy()
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLAS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        env[name] = "1"
    print("Running: " + " ".join([rscript_bin, script]), flush=True)
    result = subprocess.run([rscript_bin, script], check=False, env=env)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {' '.join([rscript_bin, script])}")


def run_microarray_bigsnpr_model(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
    *,
    chromosomes: list[str] | None = None,
    min_maf: float = 0.05,
    min_mac: int = 20,
    min_call_rate: float = 0.98,
    hwe_p_control: float = 1e-6,
    eligibility_flag: str = "primary_model_eligible",
    label: str | None = None,
    plink_prefix: str | None = None,
    copy_plink_to: str | None = None,
    overwrite_plink_copy: bool = False,
    threads: int | None = None,
    memory_mb: int | None = None,
    plink2_bin: str = "plink2",
    rscript_bin: str = "Rscript",
    alphas: tuple[float, ...] = BIGSNPR_DEFAULT_ALPHAS,
    folds: int = 5,
    nlambda: int = 100,
    dfmax: int = 50000,
    prepare_only: bool = False,
    reuse_plink_subset: bool = True,
) -> dict[str, str | pd.DataFrame | dict[str, int]]:
    chromosome_values = _normalize_autosomal_chromosomes(chromosomes or ["22"])
    label = label or microarray_bigsnpr_default_label(chromosome_values, min_maf=min_maf)
    ensure_parent_dir(microarray_bigsnpr_qc_path(paths, label))

    if copy_plink_to:
        local_prefix = stage_microarray_plink_files(config, copy_plink_to, overwrite=overwrite_plink_copy)
    else:
        local_prefix = _expand_local_prefix(plink_prefix or str(Path.home() / "plink_microarray" / "arrays"))
    if not _plink_files_exist(local_prefix):
        raise RuntimeError(
            f"Missing local PLINK BED/BIM/FAM at prefix {local_prefix}. "
            "Pass --copy-plink-to $HOME/plink_microarray or --plink-prefix pointing at an existing local arrays prefix."
        )
    if shutil.which(plink2_bin) is None and not Path(plink2_bin).exists():
        raise RuntimeError(f"Could not find PLINK2 binary `{plink2_bin}` on PATH.")

    fam_df = read_plink_fam(local_prefix)
    metadata, covariates, dropped_covariates, sample_counts = write_microarray_bigsnpr_inputs(
        config,
        matched_df,
        fam_df,
        paths,
        label=label,
        eligibility_flag=eligibility_flag,
    )
    commands = build_microarray_bigsnpr_plink_commands(
        plink2_bin=plink2_bin,
        plink_prefix=local_prefix,
        paths=paths,
        label=label,
        chromosomes=chromosome_values,
        min_maf=min_maf,
        min_mac=min_mac,
        min_call_rate=min_call_rate,
        hwe_p_control=hwe_p_control,
        threads=threads,
        memory_mb=memory_mb,
    )
    if reuse_plink_subset and _subset_files_exist(paths, label):
        print(f"Using existing bigsnpr PLINK subset for {label}.", flush=True)
    else:
        for command in commands:
            _run_command(command)

    variant_qc = _variant_qc_summary(
        plink_prefix=local_prefix,
        paths=paths,
        label=label,
        chromosomes=chromosome_values,
        min_maf=min_maf,
        min_mac=min_mac,
        min_call_rate=min_call_rate,
        hwe_p_control=hwe_p_control,
    )
    write_dataframe(variant_qc, microarray_bigsnpr_variant_qc_summary_path(paths, label))
    ncores = int(threads or 1)
    script = write_bigsnpr_r_script(
        paths,
        label=label,
        outcome_column=config.analysis.matched_outcome_column,
        covariates=covariates,
        alphas=alphas,
        folds=folds,
        nlambda=nlambda,
        dfmax=dfmax,
        ncores=ncores,
    )
    qc = {
        "label": label,
        "genotype_source": "microarray_plink_bigsnpr",
        "microarray_plink_prefix": local_prefix,
        "chromosomes": chromosome_values,
        "eligibility_flag": eligibility_flag,
        "sample_counts": sample_counts,
        "min_maf": float(min_maf),
        "min_mac": int(min_mac),
        "min_call_rate": float(min_call_rate),
        "hwe_p_control": float(hwe_p_control),
        "covariates_used": covariates,
        "dropped_covariates": dropped_covariates,
        "alphas": list(alphas),
        "folds": int(folds),
        "nlambda": int(nlambda),
        "dfmax": int(dfmax),
        "prepare_only": bool(prepare_only),
        "r_script": script,
    }
    write_json(qc, microarray_bigsnpr_qc_path(paths, label))
    write_stage_report(
        title="Stage 4: Microarray bigsnpr Sparse Model",
        summary_lines=[
            f"- Label: {label}",
            "- Genotype source: AoU v8 microarray PLINK BED/BIM/FAM",
            "- PLINK2 prepares a train-QC'd train+test subset before R import.",
            f"- Chromosomes: {', '.join(chromosome_values)}",
            (
                "- Train FAM-overlap cases/controls: "
                f"{sample_counts.get('after_microarray_fam_overlap_train_cases', 0)} / "
                f"{sample_counts.get('after_microarray_fam_overlap_train_controls', 0)}"
            ),
            (
                "- Test FAM-overlap cases/controls: "
                f"{sample_counts.get('after_microarray_fam_overlap_test_cases', 0)} / "
                f"{sample_counts.get('after_microarray_fam_overlap_test_controls', 0)}"
            ),
            f"- Variant QC: MAF >= {min_maf}, MAC >= {min_mac}, call rate >= {min_call_rate}, train-control HWE p >= {hwe_p_control}",
            f"- R script: {script}",
        ],
        preview_df=variant_qc,
        preview_columns=["filter", "threshold", "rows_before", "rows_after", "rows_removed"],
        path=microarray_bigsnpr_report_path(paths, label),
    )

    if not prepare_only:
        if shutil.which(rscript_bin) is None and not Path(rscript_bin).exists():
            raise RuntimeError(f"Could not find Rscript binary `{rscript_bin}`. Re-run with --prepare-only or install R.")
        _run_rscript_command(rscript_bin, script)

    return {
        "metadata": metadata,
        "variant_qc": variant_qc,
        "sample_counts": sample_counts,
        "qc": qc,
        "label": label,
        "r_script": script,
    }


__all__ = [
    "BIGSNPR_DEFAULT_ALPHAS",
    "build_microarray_bigsnpr_plink_commands",
    "microarray_bigsnpr_default_label",
    "microarray_bigsnpr_metrics_path",
    "microarray_bigsnpr_output_dir",
    "microarray_bigsnpr_predictions_path",
    "microarray_bigsnpr_qc_path",
    "microarray_bigsnpr_report_path",
    "microarray_bigsnpr_risk_decile_path",
    "microarray_bigsnpr_score_distribution_path",
    "microarray_bigsnpr_script_path",
    "microarray_bigsnpr_selected_variants_path",
    "microarray_bigsnpr_variant_qc_summary_path",
    "run_microarray_bigsnpr_model",
    "write_bigsnpr_r_script",
    "write_microarray_bigsnpr_inputs",
]
