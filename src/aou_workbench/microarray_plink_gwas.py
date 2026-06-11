"""Microarray GWAS pilot using local PLINK files."""

from __future__ import annotations

import math
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pandas as pd

from .config import ProjectConfig
from .io_utils import ensure_parent_dir, slugify, write_dataframe, write_json
from .paths import ProjectPaths, join_path
from .reporting import write_stage_report
from .stage4_hail_gwas import _hail_pilot_sample_frame, _lead_hit_subset, _normalize_autosomal_chromosomes
from .statistics import bh_fdr, genomic_control_lambda
from .svg import write_manhattan_svg, write_qq_svg

PLINK_COVARIATES = ("age_at_index", "is_female", "pc1", "pc2", "pc3", "pc4", "pc5")
PLINK_RESULT_COLUMNS = [
    "variant_id",
    "chromosome",
    "position",
    "ref",
    "alt",
    "effect_allele",
    "a1_freq",
    "n_samples",
    "beta",
    "se",
    "z_stat",
    "regression_p",
    "odds_ratio",
    "ci_lower",
    "ci_upper",
]


def microarray_plink_default_label(chromosomes: list[str], *, min_maf: float = 0.05) -> str:
    chrom_label = "-".join(chromosomes)
    maf_text = f"{min_maf:g}"
    maf_label = maf_text[2:] if maf_text.startswith("0.") else maf_text.replace(".", "")
    return f"microarray_plink_chr{chrom_label}_maf{maf_label}_train_qc"


def microarray_plink_output_dir(paths: ProjectPaths, label: str) -> str:
    return join_path(paths.run_root, "stage4", "microarray_plink", slugify(label))


def microarray_plink_results_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_plink_output_dir(paths, label), "gwas_results.tsv")


def microarray_plink_lead_hits_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_plink_output_dir(paths, label), "lead_hits.tsv")


def microarray_plink_qc_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_plink_output_dir(paths, label), "gwas_qc.json")


def microarray_plink_variant_qc_summary_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_plink_output_dir(paths, label), "variant_qc_summary.tsv")


def microarray_plink_report_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_plink_output_dir(paths, label), "report.md")


def microarray_plink_manhattan_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_plink_output_dir(paths, label), "manhattan.svg")


def microarray_plink_qq_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_plink_output_dir(paths, label), "qq.svg")


def microarray_plink_keep_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_plink_output_dir(paths, label), "plink_keep.tsv")


def microarray_plink_controls_keep_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_plink_output_dir(paths, label), "plink_controls_keep.tsv")


def microarray_plink_pheno_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_plink_output_dir(paths, label), "plink_pheno.tsv")


def microarray_plink_covar_path(paths: ProjectPaths, label: str) -> str:
    return join_path(microarray_plink_output_dir(paths, label), "plink_covar.tsv")


def _output_prefix(paths: ProjectPaths, label: str, name: str) -> str:
    return join_path(microarray_plink_output_dir(paths, label), name)


def _requester_pays_project(config: ProjectConfig) -> str | None:
    return (
        config.workbench.requester_pays_project
        or os.getenv("GCS_REQUESTER_PAYS_PROJECT")
        or os.getenv("GOOGLE_PROJECT")
        or os.getenv("GOOGLE_CLOUD_PROJECT")
    )


def _run_command(cmd: list[str]) -> None:
    print("Running: " + " ".join(cmd), flush=True)
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {' '.join(cmd)}")


def _expand_local_prefix(prefix: str) -> str:
    if prefix.startswith("gs://"):
        raise RuntimeError(
            "PLINK2 needs local PLINK files. Use --copy-plink-to to copy the AoU microarray BED/BIM/FAM first, "
            "or pass --plink-prefix pointing to an existing local prefix."
        )
    return str(Path(prefix).expanduser())


def _plink_files_exist(prefix: str) -> bool:
    return all(Path(f"{prefix}.{suffix}").exists() for suffix in ("bed", "bim", "fam"))


def stage_microarray_plink_files(
    config: ProjectConfig,
    local_dir: str,
    *,
    overwrite: bool = False,
) -> str:
    """Copy AoU microarray PLINK files to a local directory and return the local prefix."""

    destination = Path(local_dir).expanduser()
    destination.mkdir(parents=True, exist_ok=True)
    local_prefix = str(destination / "arrays")
    if _plink_files_exist(local_prefix) and not overwrite:
        print(f"Using existing local PLINK files at {local_prefix}.", flush=True)
        return local_prefix

    remote_prefix = config.workbench.microarray_plink_prefix
    if not remote_prefix.startswith("gs://"):
        raise RuntimeError(f"microarray_plink_prefix must be a gs:// prefix when copying: {remote_prefix}")

    cmd = ["gsutil", "-m"]
    project = _requester_pays_project(config)
    if project:
        cmd.extend(["-u", project])
    cmd.extend(
        [
            "cp",
            f"{remote_prefix}.bed",
            f"{remote_prefix}.bim",
            f"{remote_prefix}.fam",
            str(destination),
        ]
    )
    _run_command(cmd)
    if not _plink_files_exist(local_prefix):
        raise RuntimeError(f"Expected PLINK files were not copied to {local_prefix}.")
    return local_prefix


def read_plink_fam(prefix: str) -> pd.DataFrame:
    fam_path = Path(f"{prefix}.fam")
    if not fam_path.exists():
        raise RuntimeError(f"Missing PLINK FAM file: {fam_path}")
    return pd.read_csv(
        fam_path,
        sep=r"\s+",
        header=None,
        names=["FID", "IID", "father", "mother", "sex", "phenotype"],
        dtype=str,
    )


def _write_space_delimited(frame: pd.DataFrame, path: str, *, header: bool) -> None:
    ensure_parent_dir(path)
    frame.to_csv(path, sep=" ", index=False, header=header)


def write_microarray_plink_sample_files(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    fam_df: pd.DataFrame,
    paths: ProjectPaths,
    *,
    label: str,
    analysis_split: str = "train",
    eligibility_flag: str = "primary_model_eligible",
    covariates: tuple[str, ...] = PLINK_COVARIATES,
) -> tuple[pd.DataFrame, list[str], list[str], list[str], dict[str, int]]:
    """Create PLINK keep, phenotype, and covariate files for the pilot sample."""

    sample_df, kept_covariates, raw_covariates, dropped_covariates, counts, _ = _hail_pilot_sample_frame(
        matched_df,
        config,
        analysis_split=analysis_split,
        eligibility_flag=eligibility_flag,
        covariates=covariates,
        restrict_to_wgs_manifest=False,
    )
    outcome = config.analysis.matched_outcome_column
    fam_ids = set(fam_df["IID"].astype(str))
    sample_df = sample_df[sample_df["person_id"].astype(str).isin(fam_ids)].copy()
    sample_df["FID"] = "0"
    sample_df["IID"] = sample_df["person_id"].astype(str)
    sample_df[outcome] = pd.to_numeric(sample_df[outcome], errors="coerce").astype(int)
    counts["after_microarray_fam_overlap_participants"] = int(sample_df["IID"].nunique())
    counts["after_microarray_fam_overlap_cases"] = int(sample_df[outcome].eq(1).sum())
    counts["after_microarray_fam_overlap_controls"] = int(sample_df[outcome].eq(0).sum())

    keep = sample_df[["FID", "IID"]].copy()
    controls_keep = sample_df.loc[sample_df[outcome] == 0, ["FID", "IID"]].copy()
    pheno = sample_df[["FID", "IID", outcome]].copy()
    covar = sample_df[["FID", "IID", *kept_covariates]].copy()
    _write_space_delimited(keep, microarray_plink_keep_path(paths, label), header=False)
    _write_space_delimited(controls_keep, microarray_plink_controls_keep_path(paths, label), header=False)
    _write_space_delimited(pheno, microarray_plink_pheno_path(paths, label), header=True)
    _write_space_delimited(covar, microarray_plink_covar_path(paths, label), header=True)
    return sample_df, kept_covariates, raw_covariates, dropped_covariates, counts


def _plink_common_options(threads: int | None, memory_mb: int | None) -> list[str]:
    options: list[str] = []
    if threads:
        options.extend(["--threads", str(int(threads))])
    if memory_mb:
        options.extend(["--memory", str(int(memory_mb))])
    return options


def _chromosome_arg(chromosomes: list[str]) -> str:
    return ",".join(chromosomes)


def build_microarray_plink_commands(
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
    outcome_column: str,
    covariates: list[str],
    threads: int | None = None,
    memory_mb: int | None = None,
) -> tuple[list[str], list[str], list[str]]:
    common = _plink_common_options(threads, memory_mb)
    chr_arg = _chromosome_arg(chromosomes)
    max_missing = max(0.0, min(1.0, 1.0 - float(min_call_rate)))
    analysis_qc_prefix = _output_prefix(paths, label, "chr_analysis_qc")
    control_hwe_prefix = _output_prefix(paths, label, "chr_control_hwe_qc")
    assoc_prefix = _output_prefix(paths, label, "microarray_rhabdo")

    analysis_qc_cmd = [
        plink2_bin,
        "--bfile",
        plink_prefix,
        "--chr",
        chr_arg,
        "--snps-only",
        "just-acgt",
        "--keep",
        microarray_plink_keep_path(paths, label),
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
        microarray_plink_controls_keep_path(paths, label),
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
    assoc_cmd = [
        plink2_bin,
        "--bfile",
        plink_prefix,
        "--chr",
        chr_arg,
        "--keep",
        microarray_plink_keep_path(paths, label),
        "--extract",
        f"{control_hwe_prefix}.snplist",
        "--pheno",
        microarray_plink_pheno_path(paths, label),
        "--pheno-name",
        outcome_column,
        "--1",
        "--covar",
        microarray_plink_covar_path(paths, label),
        "--covar-name",
        ",".join(covariates),
        "--covar-variance-standardize",
        "--glm",
        "hide-covar",
        "firth-fallback",
        "cols=+a1freq,+nobs,+orbeta,+ci",
        "--out",
        assoc_prefix,
        *common,
    ]
    return analysis_qc_cmd, control_hwe_cmd, assoc_cmd


def _count_snplist(path: str) -> int:
    file_path = Path(path)
    if not file_path.exists():
        return 0
    with file_path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def _bim_row_counts(prefix: str, chromosomes: list[str]) -> tuple[int, int]:
    bim_path = Path(f"{prefix}.bim")
    if not bim_path.exists():
        return 0, 0
    chrom_set = set(str(chrom).replace("chr", "") for chrom in chromosomes)
    total = 0
    selected = 0
    with bim_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            total += 1
            chrom = line.split(maxsplit=1)[0].replace("chr", "")
            if chrom in chrom_set:
                selected += 1
    return total, selected


def _find_plink_result_file(paths: ProjectPaths, label: str, outcome_column: str) -> str | None:
    assoc_prefix = _output_prefix(paths, label, "microarray_rhabdo")
    matches = sorted(Path(assoc_prefix).parent.glob(f"{Path(assoc_prefix).name}.{outcome_column}.glm.logistic*"))
    if not matches:
        matches = sorted(Path(assoc_prefix).parent.glob(f"{Path(assoc_prefix).name}*.glm.logistic*"))
    return str(matches[0]) if matches else None


def parse_plink_glm_results(path: str | None) -> pd.DataFrame:
    if not path or not Path(path).exists():
        return pd.DataFrame(columns=PLINK_RESULT_COLUMNS)
    raw = pd.read_csv(path, sep=r"\s+", dtype=str)
    if raw.empty:
        return pd.DataFrame(columns=PLINK_RESULT_COLUMNS)
    if "TEST" in raw.columns:
        raw = raw[raw["TEST"].astype(str).str.upper() == "ADD"].copy()
    if raw.empty:
        return pd.DataFrame(columns=PLINK_RESULT_COLUMNS)

    def first_column(*names: str) -> str | None:
        for name in names:
            if name in raw.columns:
                return name
        return None

    chrom_col = first_column("#CHROM", "CHROM")
    pos_col = first_column("POS", "BP")
    id_col = first_column("ID", "SNP")
    ref_col = first_column("REF", "A2")
    alt_col = first_column("ALT")
    a1_col = first_column("A1")
    freq_col = first_column("A1_FREQ", "ALT_FREQS")
    n_col = first_column("OBS_CT", "N")
    beta_col = first_column("BETA")
    or_col = first_column("OR")
    se_col = first_column("SE", "LOG(OR)_SE")
    z_col = first_column("Z_STAT", "T_STAT")
    p_col = first_column("P", "P_VALUE")
    ci_lower_col = first_column("L95", "OR_L95", "BETA_L95")
    ci_upper_col = first_column("U95", "OR_U95", "BETA_U95")

    full = pd.DataFrame()
    full["variant_id"] = raw[id_col].astype(str) if id_col else ""
    full["chromosome"] = raw[chrom_col].astype(str).str.replace("chr", "", regex=False) if chrom_col else ""
    full["position"] = pd.to_numeric(raw[pos_col], errors="coerce") if pos_col else np.nan
    full["ref"] = raw[ref_col].astype(str) if ref_col else ""
    full["alt"] = raw[alt_col].astype(str) if alt_col else ""
    full["effect_allele"] = raw[a1_col].astype(str) if a1_col else ""
    full["a1_freq"] = pd.to_numeric(raw[freq_col], errors="coerce") if freq_col else np.nan
    full["n_samples"] = pd.to_numeric(raw[n_col], errors="coerce") if n_col else np.nan

    odds_ratio = pd.to_numeric(raw[or_col], errors="coerce") if or_col else pd.Series(np.nan, index=raw.index)
    beta = pd.to_numeric(raw[beta_col], errors="coerce") if beta_col else np.log(odds_ratio.where(odds_ratio > 0))
    full["beta"] = beta
    full["se"] = pd.to_numeric(raw[se_col], errors="coerce") if se_col else np.nan
    full["z_stat"] = pd.to_numeric(raw[z_col], errors="coerce") if z_col else np.nan
    full["regression_p"] = pd.to_numeric(raw[p_col], errors="coerce") if p_col else np.nan
    full["odds_ratio"] = odds_ratio.fillna(np.exp(beta))
    full["ci_lower"] = pd.to_numeric(raw[ci_lower_col], errors="coerce") if ci_lower_col else np.nan
    full["ci_upper"] = pd.to_numeric(raw[ci_upper_col], errors="coerce") if ci_upper_col else np.nan
    full = full[pd.to_numeric(full["regression_p"], errors="coerce").notna()].copy()
    full["position"] = pd.to_numeric(full["position"], errors="coerce").astype("Int64")
    return full[PLINK_RESULT_COLUMNS].reset_index(drop=True)


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
    n_results: int,
) -> pd.DataFrame:
    total_bim_rows, chromosome_rows = _bim_row_counts(plink_prefix, chromosomes)
    analysis_qc_prefix = _output_prefix(paths, label, "chr_analysis_qc")
    control_hwe_prefix = _output_prefix(paths, label, "chr_control_hwe_qc")
    analysis_qc_rows = _count_snplist(f"{analysis_qc_prefix}.snplist")
    hwe_rows = _count_snplist(f"{control_hwe_prefix}.snplist")
    rows = [
        {
            "filter": "chromosome_restriction",
            "threshold": f"chr {','.join(chromosomes)}",
            "rows_before": total_bim_rows,
            "rows_after": chromosome_rows,
            "rows_removed": max(total_bim_rows - chromosome_rows, 0),
        },
        {
            "filter": "plink_analysis_qc",
            "threshold": f"biallelic A/C/G/T SNP, MAF >= {min_maf}, MAC >= {min_mac}, call rate >= {min_call_rate}",
            "rows_before": chromosome_rows,
            "rows_after": analysis_qc_rows,
            "rows_removed": max(chromosome_rows - analysis_qc_rows, 0),
        },
        {
            "filter": "control_hwe_p",
            "threshold": f">= {hwe_p_control}",
            "rows_before": analysis_qc_rows,
            "rows_after": hwe_rows,
            "rows_removed": max(analysis_qc_rows - hwe_rows, 0),
        },
        {
            "filter": "plink_glm_result_rows",
            "threshold": "additive logistic regression result available",
            "rows_before": hwe_rows,
            "rows_after": n_results,
            "rows_removed": max(hwe_rows - n_results, 0),
        },
    ]
    return pd.DataFrame(rows)


def _postprocess_plink_results(
    full: pd.DataFrame,
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
    plink_prefix: str,
    covariates_used: list[str],
    dropped_covariates: list[str],
    plink_result_file: str | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    stage = config.analysis.stage4
    full = full.copy()
    if not full.empty:
        full = full.sort_values(["regression_p", "variant_id"]).reset_index(drop=True)
        full["fdr_q"] = bh_fdr(full["regression_p"].values)
        full["minus_log10_p"] = -np.log10(full["regression_p"].clip(lower=1e-300))

    lead_hits = _lead_hit_subset(full, stage.lead_hit_window_bp) if stage is not None else full.head(0).copy()
    write_dataframe(full, microarray_plink_results_path(paths, label))
    write_dataframe(lead_hits, microarray_plink_lead_hits_path(paths, label))
    write_dataframe(variant_qc_summary, microarray_plink_variant_qc_summary_path(paths, label))
    write_manhattan_svg(full, microarray_plink_manhattan_path(paths, label))
    write_qq_svg(full, microarray_plink_qq_path(paths, label))
    write_json(
        {
            "pilot_label": label,
            "genotype_source": "microarray_plink",
            "microarray_plink_prefix": plink_prefix,
            "plink_result_file": plink_result_file,
            "chromosomes_tested": chromosomes,
            "analysis_split": analysis_split,
            "eligibility_flag": eligibility_flag,
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
        microarray_plink_qc_path(paths, label),
    )
    write_stage_report(
        title="Stage 4: Microarray PLINK GWAS",
        summary_lines=[
            f"- Pilot label: {label}",
            "- Genotype source: AoU v8 microarray PLINK BED/BIM/FAM",
            "- Case-control definition: broad rhabdomyolysis cases versus matched eligible controls",
            f"- Analysis split: {analysis_split}",
            f"- Eligibility flag: {eligibility_flag}",
            f"- Chromosomes: {', '.join(chromosomes)}",
            f"- Complete-case samples before FAM overlap: {sample_counts.get('after_complete_case_participants', 0)}",
            f"- Microarray FAM overlap samples: {sample_counts.get('after_microarray_fam_overlap_participants', 0)}",
            (
                "- Microarray FAM overlap cases/controls: "
                f"{sample_counts.get('after_microarray_fam_overlap_cases', 0)} / "
                f"{sample_counts.get('after_microarray_fam_overlap_controls', 0)}"
            ),
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
        path=microarray_plink_report_path(paths, label),
    )
    return full, lead_hits


def run_microarray_plink_gwas(
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
    plink_prefix: str | None = None,
    copy_plink_to: str | None = None,
    overwrite_plink_copy: bool = False,
    threads: int | None = None,
    memory_mb: int | None = None,
    plink2_bin: str = "plink2",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    chromosome_values = _normalize_autosomal_chromosomes(chromosomes or ["22"])
    label = label or microarray_plink_default_label(chromosome_values, min_maf=min_maf)
    ensure_parent_dir(microarray_plink_results_path(paths, label))

    if copy_plink_to:
        local_prefix = stage_microarray_plink_files(config, copy_plink_to, overwrite=overwrite_plink_copy)
    else:
        local_prefix = _expand_local_prefix(plink_prefix or str(Path.home() / "plink_microarray" / "arrays"))
    if not _plink_files_exist(local_prefix):
        raise RuntimeError(
            f"Missing local PLINK BED/BIM/FAM at prefix {local_prefix}. "
            "Either pass --copy-plink-to $HOME/plink_microarray to copy the AoU files, "
            "or pass --plink-prefix pointing at an existing local arrays prefix."
        )
    if shutil.which(plink2_bin) is None and not Path(plink2_bin).exists():
        raise RuntimeError(f"Could not find PLINK2 binary `{plink2_bin}` on PATH.")

    fam_df = read_plink_fam(local_prefix)
    sample_df, covariates_used, raw_covariates, dropped_covariates, sample_counts = write_microarray_plink_sample_files(
        config,
        matched_df,
        fam_df,
        paths,
        label=label,
        analysis_split=analysis_split,
        eligibility_flag=eligibility_flag,
    )
    if sample_df.empty:
        empty_qc = _variant_qc_summary(
            plink_prefix=local_prefix,
            paths=paths,
            label=label,
            chromosomes=chromosome_values,
            min_maf=min_maf,
            min_mac=min_mac,
            min_call_rate=min_call_rate,
            hwe_p_control=hwe_p_control,
            n_results=0,
        )
        return _postprocess_plink_results(
            pd.DataFrame(columns=PLINK_RESULT_COLUMNS),
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
            plink_prefix=local_prefix,
            covariates_used=raw_covariates,
            dropped_covariates=dropped_covariates,
            plink_result_file=None,
        )

    commands = build_microarray_plink_commands(
        plink2_bin=plink2_bin,
        plink_prefix=local_prefix,
        paths=paths,
        label=label,
        chromosomes=chromosome_values,
        min_maf=min_maf,
        min_mac=min_mac,
        min_call_rate=min_call_rate,
        hwe_p_control=hwe_p_control,
        outcome_column=config.analysis.matched_outcome_column,
        covariates=covariates_used,
        threads=threads,
        memory_mb=memory_mb,
    )
    for command in commands:
        _run_command(command)

    result_file = _find_plink_result_file(paths, label, config.analysis.matched_outcome_column)
    full = parse_plink_glm_results(result_file)
    variant_qc_summary = _variant_qc_summary(
        plink_prefix=local_prefix,
        paths=paths,
        label=label,
        chromosomes=chromosome_values,
        min_maf=min_maf,
        min_mac=min_mac,
        min_call_rate=min_call_rate,
        hwe_p_control=hwe_p_control,
        n_results=int(full.shape[0]),
    )
    return _postprocess_plink_results(
        full,
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
        plink_prefix=local_prefix,
        covariates_used=raw_covariates,
        dropped_covariates=dropped_covariates,
        plink_result_file=result_file,
    )


__all__ = [
    "PLINK_COVARIATES",
    "build_microarray_plink_commands",
    "microarray_plink_covar_path",
    "microarray_plink_default_label",
    "microarray_plink_keep_path",
    "microarray_plink_lead_hits_path",
    "microarray_plink_manhattan_path",
    "microarray_plink_pheno_path",
    "microarray_plink_qc_path",
    "microarray_plink_qq_path",
    "microarray_plink_report_path",
    "microarray_plink_results_path",
    "microarray_plink_variant_qc_summary_path",
    "parse_plink_glm_results",
    "read_plink_fam",
    "run_microarray_plink_gwas",
    "stage_microarray_plink_files",
    "write_microarray_plink_sample_files",
]
