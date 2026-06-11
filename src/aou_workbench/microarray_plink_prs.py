"""LD-clumped PRS scoring from microarray PLINK GWAS results."""

from __future__ import annotations

import math
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from .config import ProjectConfig
from .io_utils import ensure_parent_dir, read_table, slugify, write_dataframe, write_json, write_text
from .microarray_plink_gwas import (
    _expand_local_prefix,
    _plink_common_options,
    _plink_files_exist,
    _run_command,
    microarray_plink_output_dir,
    microarray_plink_results_path,
    parse_plink_glm_results,
    read_plink_fam,
)
from .paths import ProjectPaths, join_path
from .reporting import write_stage_report
from .stage4_hail_gwas import _pilot_case_control_definition_mask
from .statistics import run_binary_logistic_regression

DEFAULT_PRS_THRESHOLDS = (5e-8, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0)


def microarray_prs_default_label(score_split: str) -> str:
    return f"{score_split}_clumped_threshold_grid"


def microarray_prs_output_dir(paths: ProjectPaths, gwas_label: str, prs_label: str) -> str:
    return join_path(microarray_plink_output_dir(paths, gwas_label), "prs", slugify(prs_label))


def microarray_prs_weights_path(paths: ProjectPaths, gwas_label: str, prs_label: str) -> str:
    return join_path(microarray_prs_output_dir(paths, gwas_label, prs_label), "prs_weights.tsv")


def microarray_prs_range_path(paths: ProjectPaths, gwas_label: str, prs_label: str) -> str:
    return join_path(microarray_prs_output_dir(paths, gwas_label, prs_label), "prs_pvalue_ranges.tsv")


def microarray_prs_keep_path(paths: ProjectPaths, gwas_label: str, prs_label: str, score_split: str) -> str:
    return join_path(microarray_prs_output_dir(paths, gwas_label, prs_label), f"{score_split}_keep.tsv")


def microarray_prs_scores_path(paths: ProjectPaths, gwas_label: str, prs_label: str) -> str:
    return join_path(microarray_prs_output_dir(paths, gwas_label, prs_label), "prs_scores.tsv")


def microarray_prs_metrics_path(paths: ProjectPaths, gwas_label: str, prs_label: str) -> str:
    return join_path(microarray_prs_output_dir(paths, gwas_label, prs_label), "prs_metrics.tsv")


def microarray_prs_qc_path(paths: ProjectPaths, gwas_label: str, prs_label: str) -> str:
    return join_path(microarray_prs_output_dir(paths, gwas_label, prs_label), "prs_qc.json")


def microarray_prs_report_path(paths: ProjectPaths, gwas_label: str, prs_label: str) -> str:
    return join_path(microarray_prs_output_dir(paths, gwas_label, prs_label), "report.md")


def microarray_prs_case_status_svg_path(paths: ProjectPaths, gwas_label: str, prs_label: str) -> str:
    return join_path(microarray_prs_output_dir(paths, gwas_label, prs_label), "prs_by_case_status.svg")


def _prs_prefix(paths: ProjectPaths, gwas_label: str, prs_label: str, name: str) -> str:
    return join_path(microarray_prs_output_dir(paths, gwas_label, prs_label), name)


def parse_thresholds(value: str | tuple[float, ...] | list[float]) -> list[float]:
    if isinstance(value, str):
        thresholds = [float(item.strip()) for item in value.split(",") if item.strip()]
    else:
        thresholds = [float(item) for item in value]
    thresholds = sorted(set(thresholds))
    if not thresholds:
        raise ValueError("At least one PRS p-value threshold is required.")
    if thresholds[0] < 0 or thresholds[-1] > 1:
        raise ValueError("PRS p-value thresholds must be between 0 and 1.")
    return thresholds


def _threshold_label(value: float) -> str:
    text = f"{value:g}".replace(".", "_").replace("+", "")
    return f"p{text}"


def _load_gwas_results(paths: ProjectPaths, gwas_label: str) -> tuple[pd.DataFrame, str]:
    processed = Path(microarray_plink_results_path(paths, gwas_label))
    if processed.exists():
        frame = read_table(str(processed))
        if not frame.empty:
            return frame, str(processed)

    gwas_dir = Path(microarray_plink_output_dir(paths, gwas_label))
    matches = sorted(gwas_dir.glob("microarray_rhabdo.analysis_case.glm.logistic*"))
    if not matches:
        matches = sorted(gwas_dir.glob("microarray_rhabdo*.glm.logistic*"))
    raw_path = str(matches[0]) if matches else None
    frame = parse_plink_glm_results(raw_path)
    if frame.empty:
        raise RuntimeError(
            f"No GWAS results found for label `{gwas_label}`. Expected {processed} "
            "or a raw microarray_rhabdo*.glm.logistic* file in the GWAS output directory."
        )
    return frame, str(raw_path)


def _valid_gwas_weights(gwas: pd.DataFrame) -> pd.DataFrame:
    required = ["variant_id", "effect_allele", "beta", "regression_p"]
    missing = [column for column in required if column not in gwas.columns]
    if missing:
        raise RuntimeError("GWAS results are missing PRS weight columns: " + ", ".join(missing))
    weights = gwas[required].copy()
    weights["variant_id"] = weights["variant_id"].astype(str)
    weights["effect_allele"] = weights["effect_allele"].astype(str)
    weights["beta"] = pd.to_numeric(weights["beta"], errors="coerce")
    weights["regression_p"] = pd.to_numeric(weights["regression_p"], errors="coerce")
    bad_text = weights["variant_id"].isin({"", "nan", "None"}) | weights["effect_allele"].isin({"", "nan", "None"})
    weights = weights[~bad_text].dropna(subset=["beta", "regression_p"]).copy()
    weights = weights[(weights["regression_p"] >= 0) & (weights["regression_p"] <= 1)].copy()
    weights = weights.drop_duplicates(subset=["variant_id"], keep=False).copy()
    weights = weights.rename(columns={"variant_id": "ID", "effect_allele": "A1", "beta": "BETA", "regression_p": "P"})
    return weights.sort_values(["P", "ID"]).reset_index(drop=True)


def _write_prs_sample_keep(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    fam_df: pd.DataFrame,
    paths: ProjectPaths,
    *,
    gwas_label: str,
    prs_label: str,
    score_split: str,
    eligibility_flag: str,
) -> tuple[pd.DataFrame, dict[str, int]]:
    outcome = config.analysis.matched_outcome_column
    required = {"person_id", outcome, "analysis_split", eligibility_flag}
    missing = sorted(required.difference(matched_df.columns))
    if missing:
        raise RuntimeError("Matched cohort is missing PRS sample columns: " + ", ".join(missing))
    sample = matched_df.copy()
    sample["person_id"] = sample["person_id"].astype(str)
    counts = {"matched_input_participants": int(sample["person_id"].nunique())}
    sample = sample[sample["analysis_split"].astype(str) == str(score_split)].copy()
    counts["after_score_split_participants"] = int(sample["person_id"].nunique())
    eligible = pd.to_numeric(sample[eligibility_flag], errors="coerce").fillna(0)
    sample = sample[eligible == 1].copy()
    counts["after_eligibility_participants"] = int(sample["person_id"].nunique())
    sample = sample[_pilot_case_control_definition_mask(sample, outcome)].copy()
    outcome_values = pd.to_numeric(sample[outcome], errors="coerce")
    counts["after_case_control_definition_participants"] = int(sample["person_id"].nunique())
    counts["after_case_control_definition_cases"] = int(outcome_values.eq(1).sum())
    counts["after_case_control_definition_controls"] = int(outcome_values.eq(0).sum())
    fam_ids = set(fam_df["IID"].astype(str))
    sample = sample[sample["person_id"].isin(fam_ids)].copy()
    sample[outcome] = pd.to_numeric(sample[outcome], errors="coerce").astype(int)
    sample["FID"] = "0"
    sample["IID"] = sample["person_id"].astype(str)
    counts["after_microarray_fam_overlap_participants"] = int(sample["person_id"].nunique())
    counts["after_microarray_fam_overlap_cases"] = int(sample[outcome].eq(1).sum())
    counts["after_microarray_fam_overlap_controls"] = int(sample[outcome].eq(0).sum())
    ensure_parent_dir(microarray_prs_keep_path(paths, gwas_label, prs_label, score_split))
    sample[["FID", "IID"]].to_csv(
        microarray_prs_keep_path(paths, gwas_label, prs_label, score_split),
        sep=" ",
        index=False,
        header=False,
    )
    return sample, counts


def _write_clump_input(weights: pd.DataFrame, path: str) -> None:
    ensure_parent_dir(path)
    weights[["ID", "P", "A1"]].to_csv(path, sep="\t", index=False)


def _write_ranges(thresholds: list[float], path: str) -> pd.DataFrame:
    ranges = pd.DataFrame(
        [
            {
                "label": _threshold_label(threshold),
                "lower": 0.0,
                "upper": float(threshold),
            }
            for threshold in thresholds
        ]
    )
    ensure_parent_dir(path)
    ranges.to_csv(path, sep="\t", index=False, header=False)
    return ranges


def build_microarray_prs_commands(
    *,
    plink2_bin: str,
    plink_prefix: str,
    paths: ProjectPaths,
    gwas_label: str,
    prs_label: str,
    clump_split: str,
    score_split: str,
    clump_input_path: str,
    clump_r2: float,
    clump_kb: int,
    clump_p1: float,
    clump_p2: float,
    threads: int | None = None,
    memory_mb: int | None = None,
) -> tuple[list[str], list[str]]:
    common = _plink_common_options(threads, memory_mb)
    clump_prefix = _prs_prefix(paths, gwas_label, prs_label, "plink_clump")
    score_prefix = _prs_prefix(paths, gwas_label, prs_label, "plink_prs")
    clump_cmd = [
        plink2_bin,
        "--bfile",
        plink_prefix,
        "--keep",
        microarray_prs_keep_path(paths, gwas_label, prs_label, clump_split),
        "--clump",
        clump_input_path,
        "--clump-id-field",
        "ID",
        "--clump-p-field",
        "P",
        "--clump-a1-field",
        "A1",
        "--clump-r2",
        str(float(clump_r2)),
        "--clump-kb",
        str(int(clump_kb)),
        "--clump-p1",
        str(float(clump_p1)),
        "--clump-p2",
        str(float(clump_p2)),
        "--out",
        clump_prefix,
        *common,
    ]
    score_cmd = [
        plink2_bin,
        "--bfile",
        plink_prefix,
        "--keep",
        microarray_prs_keep_path(paths, gwas_label, prs_label, score_split),
        "--score",
        microarray_prs_weights_path(paths, gwas_label, prs_label),
        "1",
        "2",
        "3",
        "header-read",
        "cols=+scoresums",
        "list-variants",
        "--q-score-range",
        microarray_prs_range_path(paths, gwas_label, prs_label),
        microarray_prs_weights_path(paths, gwas_label, prs_label),
        "1",
        "4",
        "header",
        "--out",
        score_prefix,
        *common,
    ]
    return clump_cmd, score_cmd


def _clump_file(paths: ProjectPaths, gwas_label: str, prs_label: str) -> str | None:
    prefix = Path(_prs_prefix(paths, gwas_label, prs_label, "plink_clump"))
    matches = sorted(prefix.parent.glob(f"{prefix.name}.clumps*"))
    matches = [path for path in matches if "missing" not in path.name]
    return str(matches[0]) if matches else None


def parse_clumped_variant_ids(path: str | None) -> set[str]:
    if not path or not Path(path).exists():
        return set()
    frame = pd.read_csv(path, sep=r"\s+", dtype=str)
    if frame.empty:
        return set()
    id_column = "ID" if "ID" in frame.columns else "SNP" if "SNP" in frame.columns else None
    if id_column is None:
        return set()
    return set(frame[id_column].astype(str))


def _write_clumped_weights(weights: pd.DataFrame, clumped_ids: set[str], path: str) -> pd.DataFrame:
    selected = weights[weights["ID"].isin(clumped_ids)].copy() if clumped_ids else weights.head(0).copy()
    selected = selected.sort_values(["P", "ID"]).reset_index(drop=True)
    write_dataframe(selected[["ID", "A1", "BETA", "P"]], path)
    return selected


def _score_files(paths: ProjectPaths, gwas_label: str, prs_label: str) -> list[Path]:
    prefix = Path(_prs_prefix(paths, gwas_label, prs_label, "plink_prs"))
    return sorted(prefix.parent.glob(f"{prefix.name}.*.sscore"))


def _score_label_from_path(path: Path) -> str:
    prefix = path.name.split(".sscore", 1)[0]
    marker = "plink_prs."
    return prefix.split(marker, 1)[1] if marker in prefix else prefix


def _score_column(frame: pd.DataFrame) -> str:
    preferred = ["SCORE1_SUM", "SCORE1_AVG", "SCORE_SUM", "SCORE_AVG"]
    for column in preferred:
        if column in frame.columns:
            return column
    score_columns = [column for column in frame.columns if "SCORE" in column.upper()]
    if score_columns:
        return score_columns[-1]
    raise RuntimeError("Could not find a SCORE column in PLINK .sscore output.")


def parse_prs_score_files(paths: ProjectPaths, gwas_label: str, prs_label: str) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for path in _score_files(paths, gwas_label, prs_label):
        frame = pd.read_csv(path, sep=r"\s+", dtype=str)
        if frame.empty:
            continue
        iid_column = "IID" if "IID" in frame.columns else "#IID" if "#IID" in frame.columns else None
        if iid_column is None:
            raise RuntimeError(f"Could not find IID column in {path}.")
        score_column = _score_column(frame)
        threshold_label = _score_label_from_path(path)
        subset = pd.DataFrame(
            {
                "person_id": frame[iid_column].astype(str),
                "threshold_label": threshold_label,
                "prs_score": pd.to_numeric(frame[score_column], errors="coerce"),
            }
        )
        rows.append(subset)
    if not rows:
        return pd.DataFrame(columns=["person_id", "threshold_label", "prs_score"])
    return pd.concat(rows, ignore_index=True)


def _count_score_variants(paths: ProjectPaths, gwas_label: str, prs_label: str, threshold_label: str, weights: pd.DataFrame) -> int:
    prefix = Path(_prs_prefix(paths, gwas_label, prs_label, "plink_prs"))
    candidates = [
        prefix.parent / f"{prefix.name}.{threshold_label}.sscore.vars",
        prefix.parent / f"{prefix.name}.{threshold_label}.sscore.vars.zst",
    ]
    for path in candidates:
        if path.exists() and not path.name.endswith(".zst"):
            with path.open("r", encoding="utf-8") as handle:
                return sum(1 for line in handle if line.strip())
    threshold = _threshold_from_label(threshold_label)
    return int((weights["P"] <= threshold).sum()) if threshold is not None and not weights.empty else 0


def _threshold_from_label(label: str) -> float | None:
    if not label.startswith("p"):
        return None
    text = label[1:].replace("_", ".")
    try:
        return float(text)
    except ValueError:
        try:
            return float(label[1:])
        except ValueError:
            return None


def _roc_auc(y_true: np.ndarray, score: np.ndarray) -> float:
    positives = y_true == 1
    n_pos = int(positives.sum())
    n_neg = int((~positives).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = pd.Series(score).rank(method="average").to_numpy()
    rank_sum = float(ranks[positives].sum())
    return (rank_sum - (n_pos * (n_pos + 1) / 2.0)) / (n_pos * n_neg)


def prs_metrics(
    scores: pd.DataFrame,
    sample_df: pd.DataFrame,
    paths: ProjectPaths | None = None,
    gwas_label: str | None = None,
    prs_label: str | None = None,
    weights: pd.DataFrame | None = None,
) -> pd.DataFrame:
    outcome = "analysis_case"
    if outcome in scores.columns:
        merged = scores.copy()
    else:
        merged = scores.merge(sample_df[["person_id", outcome]], on="person_id", how="inner")
    rows: list[dict[str, object]] = []
    weights = weights if weights is not None else pd.DataFrame()
    for threshold_label, group in merged.groupby("threshold_label", sort=False):
        group = group.dropna(subset=["prs_score"]).copy()
        y = group[outcome].astype(int).to_numpy()
        raw_score = group["prs_score"].astype(float)
        sd = float(raw_score.std(ddof=0))
        if math.isfinite(sd) and sd > 0:
            standardized = (raw_score - float(raw_score.mean())) / sd
        else:
            standardized = raw_score * np.nan
        exposure = pd.Series(standardized.to_numpy(), index=group["person_id"].astype(str))
        regression = run_binary_logistic_regression(
            group[["person_id", outcome]],
            exposure,
            outcome_column=outcome,
            covariates=(),
        )
        variant_count = (
            _count_score_variants(paths, gwas_label, prs_label, threshold_label, weights)
            if paths is not None and gwas_label is not None and prs_label is not None
            else int((weights["P"] <= (_threshold_from_label(threshold_label) or 0)).sum()) if not weights.empty else 0
        )
        cases = group[group[outcome] == 1]
        controls = group[group[outcome] == 0]
        rows.append(
            {
                "threshold_label": threshold_label,
                "p_threshold": _threshold_from_label(threshold_label),
                "n_participants": int(group["person_id"].nunique()),
                "n_cases": int((group[outcome] == 1).sum()),
                "n_controls": int((group[outcome] == 0).sum()),
                "n_variants_scored": int(variant_count),
                "roc_auc": _roc_auc(y, raw_score.to_numpy()) if len(group) else float("nan"),
                "case_mean_prs": float(cases["prs_score"].mean()) if not cases.empty else float("nan"),
                "control_mean_prs": float(controls["prs_score"].mean()) if not controls.empty else float("nan"),
                "case_sd_prs": float(cases["prs_score"].std(ddof=0)) if not cases.empty else float("nan"),
                "control_sd_prs": float(controls["prs_score"].std(ddof=0)) if not controls.empty else float("nan"),
                "or_per_prs_sd": regression["odds_ratio"],
                "beta_per_prs_sd": regression["beta"],
                "se_per_prs_sd": regression["se"],
                "p_value": regression["regression_p"],
            }
        )
    if not rows:
        return pd.DataFrame(
            columns=[
                "threshold_label",
                "p_threshold",
                "n_participants",
                "n_cases",
                "n_controls",
                "n_variants_scored",
                "roc_auc",
                "case_mean_prs",
                "control_mean_prs",
                "case_sd_prs",
                "control_sd_prs",
                "or_per_prs_sd",
                "beta_per_prs_sd",
                "se_per_prs_sd",
                "p_value",
            ]
        )
    return pd.DataFrame(rows).sort_values("p_threshold").reset_index(drop=True)


def _write_prs_case_status_svg(metrics: pd.DataFrame, path: str) -> None:
    if metrics.empty:
        write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="240"><rect width="100%" height="100%" fill="#fffdf7"/><text x="40" y="120" font-size="20" font-family="Menlo, monospace">No PRS scores available.</text></svg>',
            path,
        )
        return
    frame = metrics.copy()
    frame["case_mean_prs"] = pd.to_numeric(frame["case_mean_prs"], errors="coerce")
    frame["control_mean_prs"] = pd.to_numeric(frame["control_mean_prs"], errors="coerce")
    values = pd.concat([frame["case_mean_prs"], frame["control_mean_prs"]]).dropna()
    if values.empty:
        y_min, y_max = -1.0, 1.0
    else:
        margin = max((float(values.max()) - float(values.min())) * 0.15, 1e-6)
        y_min, y_max = float(values.min()) - margin, float(values.max()) + margin
    width, height = 980, 430
    left, right, top, bottom = 80, 40, 40, 95
    inner_w, inner_h = width - left - right, height - top - bottom

    def x_at(index: int) -> float:
        return left + (index + 0.5) * inner_w / max(len(frame), 1)

    def y_at(value: float) -> float:
        if y_max == y_min or not math.isfinite(value):
            return top + inner_h / 2
        return top + inner_h - ((value - y_min) / (y_max - y_min)) * inner_h

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="#fffdf7"/>',
        f'<text x="{left}" y="24" font-size="22" font-family="Menlo, monospace">PRS Mean Score by Case Status</text>',
        f'<line x1="{left}" y1="{top + inner_h}" x2="{width - right}" y2="{top + inner_h}" stroke="#1f2937"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + inner_h}" stroke="#1f2937"/>',
        f'<text x="{left + 10}" y="{top + 18}" font-size="12" font-family="Menlo, monospace">mean PRS</text>',
    ]
    for index, row in enumerate(frame.itertuples(index=False)):
        x = x_at(index)
        case_y = y_at(float(row.case_mean_prs))
        control_y = y_at(float(row.control_mean_prs))
        elements.append(f'<line x1="{x:.1f}" y1="{case_y:.1f}" x2="{x:.1f}" y2="{control_y:.1f}" stroke="#9ca3af"/>')
        elements.append(f'<circle cx="{x:.1f}" cy="{case_y:.1f}" r="4.5" fill="#991b1b"/>')
        elements.append(f'<circle cx="{x:.1f}" cy="{control_y:.1f}" r="4.5" fill="#0f766e"/>')
        label = str(row.threshold_label).replace("p", "p<=")
        elements.append(
            f'<text x="{x:.1f}" y="{height - 38}" font-size="10" transform="rotate(-35 {x:.1f} {height - 38})" text-anchor="end" font-family="Menlo, monospace">{label}</text>'
        )
    elements.append(f'<circle cx="{width - 210}" cy="22" r="4.5" fill="#991b1b"/><text x="{width - 198}" y="26" font-size="12" font-family="Menlo, monospace">cases</text>')
    elements.append(f'<circle cx="{width - 135}" cy="22" r="4.5" fill="#0f766e"/><text x="{width - 123}" y="26" font-size="12" font-family="Menlo, monospace">controls</text>')
    elements.append("</svg>")
    write_text("".join(elements), path)


def run_microarray_plink_prs(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
    paths: ProjectPaths,
    *,
    gwas_label: str,
    plink_prefix: str | None = None,
    plink2_bin: str = "plink2",
    score_split: str = "test",
    eligibility_flag: str = "primary_model_eligible",
    clump_r2: float = 0.1,
    clump_kb: int = 250,
    clump_p1: float = 1.0,
    clump_p2: float = 1.0,
    thresholds: str | tuple[float, ...] | list[float] = DEFAULT_PRS_THRESHOLDS,
    label: str | None = None,
    threads: int | None = None,
    memory_mb: int | None = None,
) -> dict[str, pd.DataFrame | str]:
    prs_label = label or microarray_prs_default_label(score_split)
    output_dir = microarray_prs_output_dir(paths, gwas_label, prs_label)
    ensure_parent_dir(join_path(output_dir, "placeholder.txt"))
    local_prefix = _expand_local_prefix(plink_prefix or str(Path.home() / "plink_microarray" / "arrays"))
    if not _plink_files_exist(local_prefix):
        raise RuntimeError(f"Missing local PLINK BED/BIM/FAM at prefix {local_prefix}.")
    if shutil.which(plink2_bin) is None and not Path(plink2_bin).exists():
        raise RuntimeError(f"Could not find PLINK2 binary `{plink2_bin}` on PATH.")

    threshold_values = parse_thresholds(thresholds)
    fam_df = read_plink_fam(local_prefix)
    _, clump_sample_counts = _write_prs_sample_keep(
        config,
        matched_df,
        fam_df,
        paths,
        gwas_label=gwas_label,
        prs_label=prs_label,
        score_split="train",
        eligibility_flag=eligibility_flag,
    )
    sample_df, sample_counts = _write_prs_sample_keep(
        config,
        matched_df,
        fam_df,
        paths,
        gwas_label=gwas_label,
        prs_label=prs_label,
        score_split=score_split,
        eligibility_flag=eligibility_flag,
    )
    gwas, gwas_source = _load_gwas_results(paths, gwas_label)
    weights = _valid_gwas_weights(gwas)
    clump_input_path = _prs_prefix(paths, gwas_label, prs_label, "plink_clump_input.tsv")
    _write_clump_input(weights, clump_input_path)
    ranges = _write_ranges(threshold_values, microarray_prs_range_path(paths, gwas_label, prs_label))
    clump_cmd, score_cmd = build_microarray_prs_commands(
        plink2_bin=plink2_bin,
        plink_prefix=local_prefix,
        paths=paths,
        gwas_label=gwas_label,
        prs_label=prs_label,
        clump_split="train",
        score_split=score_split,
        clump_input_path=clump_input_path,
        clump_r2=clump_r2,
        clump_kb=clump_kb,
        clump_p1=clump_p1,
        clump_p2=clump_p2,
        threads=threads,
        memory_mb=memory_mb,
    )
    _run_command(clump_cmd)
    clumped_ids = parse_clumped_variant_ids(_clump_file(paths, gwas_label, prs_label))
    clumped_weights = _write_clumped_weights(weights, clumped_ids, microarray_prs_weights_path(paths, gwas_label, prs_label))
    if not clumped_weights.empty and not sample_df.empty:
        _run_command(score_cmd)
    scores = parse_prs_score_files(paths, gwas_label, prs_label)
    score_sample = sample_df[["person_id", config.analysis.matched_outcome_column]].rename(
        columns={config.analysis.matched_outcome_column: "analysis_case"}
    )
    if not scores.empty:
        scores = scores.merge(score_sample, on="person_id", how="left")
    metrics = prs_metrics(
        scores,
        score_sample,
        paths,
        gwas_label,
        prs_label,
        clumped_weights,
    )
    write_dataframe(scores, microarray_prs_scores_path(paths, gwas_label, prs_label))
    write_dataframe(metrics, microarray_prs_metrics_path(paths, gwas_label, prs_label))
    _write_prs_case_status_svg(metrics, microarray_prs_case_status_svg_path(paths, gwas_label, prs_label))
    write_json(
        {
            "gwas_label": gwas_label,
            "prs_label": prs_label,
            "gwas_source": gwas_source,
            "plink_prefix": local_prefix,
            "score_split": score_split,
            "clump_split": "train",
            "eligibility_flag": eligibility_flag,
            "clump_r2": float(clump_r2),
            "clump_kb": int(clump_kb),
            "clump_p1": float(clump_p1),
            "clump_p2": float(clump_p2),
            "thresholds": threshold_values,
            "sample_counts": sample_counts,
            "clump_sample_counts": clump_sample_counts,
            "n_gwas_weight_candidates": int(weights.shape[0]),
            "n_clumped_variants": int(len(clumped_ids)),
            "n_weights_written": int(clumped_weights.shape[0]),
            "n_score_rows": int(scores.shape[0]),
        },
        microarray_prs_qc_path(paths, gwas_label, prs_label),
    )
    write_stage_report(
        title="Microarray PLINK PRS",
        summary_lines=[
            f"- GWAS label: {gwas_label}",
            f"- PRS label: {prs_label}",
            f"- Scored split: {score_split}",
            "- LD clumping split: train",
            f"- Eligibility flag: {eligibility_flag}",
            "- PRS method: PLINK2 LD clumping plus p-value threshold grid",
            f"- Clumping: r2 <= {clump_r2}, window {clump_kb} kb, p1 <= {clump_p1}, p2 <= {clump_p2}",
            f"- Scored participants: {sample_counts.get('after_microarray_fam_overlap_participants', 0)}",
            f"- Scored cases/controls: {sample_counts.get('after_microarray_fam_overlap_cases', 0)} / {sample_counts.get('after_microarray_fam_overlap_controls', 0)}",
            f"- Clumping participants: {clump_sample_counts.get('after_microarray_fam_overlap_participants', 0)}",
            f"- Clumped variants written as weights: {clumped_weights.shape[0]}",
            "- Threshold comparisons are exploratory; no threshold is selected using the held-out test set.",
        ],
        preview_df=metrics,
        preview_columns=["threshold_label", "n_variants_scored", "n_participants", "roc_auc", "or_per_prs_sd", "p_value"],
        path=microarray_prs_report_path(paths, gwas_label, prs_label),
    )
    return {
        "scores": scores,
        "metrics": metrics,
        "weights": clumped_weights,
        "ranges": ranges,
        "report_path": microarray_prs_report_path(paths, gwas_label, prs_label),
        "qc_path": microarray_prs_qc_path(paths, gwas_label, prs_label),
    }


__all__ = [
    "DEFAULT_PRS_THRESHOLDS",
    "build_microarray_prs_commands",
    "microarray_prs_case_status_svg_path",
    "microarray_prs_default_label",
    "microarray_prs_keep_path",
    "microarray_prs_metrics_path",
    "microarray_prs_qc_path",
    "microarray_prs_range_path",
    "microarray_prs_report_path",
    "microarray_prs_scores_path",
    "microarray_prs_weights_path",
    "parse_clumped_variant_ids",
    "parse_prs_score_files",
    "parse_thresholds",
    "prs_metrics",
    "run_microarray_plink_prs",
]
