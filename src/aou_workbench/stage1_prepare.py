"""Stage 1 input preparation from AoU smaller genomics callsets."""

from __future__ import annotations

import io
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import pandas as pd

from .config import ProjectConfig, TargetVariant
from .io_utils import ensure_parent_dir, write_dataframe, write_json
from .preflight import apply_runtime_defaults


def _callset_root(path: str) -> str:
    clean = path.rstrip("/")
    for suffix in ("/splitMT/hail.mt", "/multiMT/hail.mt", "/hail.mt"):
        if clean.endswith(suffix):
            return clean[: -len(suffix)]
    return clean


def _split_mt_path(path: str) -> str:
    return f"{_callset_root(path)}/splitMT/hail.mt"


def _vcf_root(path: str) -> str:
    return f"{_callset_root(path)}/vcf"


def _panel_targets_frame(targets: tuple[TargetVariant, ...]) -> pd.DataFrame:
    rows = [
        {
            "label": target.label,
            "gene": target.gene,
            "rsid": target.rsid,
            "contig": target.contig,
            "position": int(target.position),
            "ref": target.ref,
            "alt": target.alt,
            "source": target.source,
            "evidence_tier": target.evidence_tier,
            "exact_test_model": target.exact_test_model,
            "variant_id": target.variant_id,
        }
        for target in targets
    ]
    return pd.DataFrame(rows)


def _matched_person_ids(matched_df: pd.DataFrame) -> list[str]:
    if "person_id" not in matched_df.columns:
        raise ValueError("Matched cohort is missing the required 'person_id' column.")
    ids = matched_df["person_id"].astype(str).dropna().drop_duplicates().sort_values()
    return ids.tolist()


def _vcf_callset_paths(config: ProjectConfig) -> dict[str, str]:
    return {
        "acaf": _vcf_root(config.workbench.acaf_mt_path),
        "exome": _vcf_root(config.workbench.exome_mt_path),
        "clinvar": _vcf_root(config.workbench.clinvar_mt_path),
    }


def _run_checked(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip() or "Command failed")
    return result


def _bcftools_binary() -> str | None:
    override = os.getenv("AOU_WORKBENCH_BCFTOOLS")
    if override:
        return override
    return shutil.which("bcftools")


def _requester_pays_project(config: ProjectConfig) -> str | None:
    return (
        config.workbench.requester_pays_project
        or os.getenv("GCS_REQUESTER_PAYS_PROJECT")
        or os.getenv("GOOGLE_PROJECT")
        or os.getenv("GOOGLE_CLOUD_PROJECT")
    )


def _gsutil_base_cmd(config: ProjectConfig) -> list[str]:
    cmd = ["gsutil"]
    requester_project = _requester_pays_project(config)
    if requester_project:
        cmd.extend(["-u", requester_project])
    return cmd


def _discover_remote_vcf_shards(config: ProjectConfig, root: str, contigs: list[str]) -> list[str]:
    listing = _run_checked([*_gsutil_base_cmd(config), "ls", f"{root}/**"]).stdout.splitlines()
    wanted: list[str] = []
    for line in listing:
        candidate = line.strip()
        lowered = candidate.lower()
        if not (lowered.endswith(".vcf.bgz") or lowered.endswith(".vcf.gz")):
            continue
        for contig in contigs:
            token = contig.lower()
            chrom = token[3:] if token.startswith("chr") else token
            if f"chr{chrom}" in lowered:
                wanted.append(candidate)
                break
            if re.search(rf"(^|[^0-9]){re.escape(chrom)}([^0-9]|$)", lowered):
                wanted.append(candidate)
                break
    return sorted(set(wanted))


def _copy_remote_objects(config: ProjectConfig, sources: list[str], dest_dir: str) -> None:
    for source in sources:
        _run_checked([*_gsutil_base_cmd(config), "cp", source, dest_dir])


def _gt_to_dosage(gt: str) -> float:
    if gt in {".", "./.", ".|."}:
        return 0.0
    tokens = re.split(r"[\/|]", gt.split(":")[0])
    if any(token == "." for token in tokens):
        return 0.0
    return float(sum(1 for token in tokens if token != "0"))


def _normalize_contig(value: str) -> str:
    text = str(value)
    return text if text.startswith("chr") else f"chr{text}"


def _extract_from_vcf_callset(
    *,
    config: ProjectConfig,
    callset_name: str,
    vcf_root: str,
    target_df: pd.DataFrame,
    sample_ids: list[str],
) -> pd.DataFrame:
    if target_df.empty or not sample_ids:
        return pd.DataFrame()
    bcftools = _bcftools_binary()
    if bcftools is None:
        raise RuntimeError(
            "bcftools is required for Stage 1 smaller-callset extraction. "
            "Set AOU_WORKBENCH_BCFTOOLS to a dedicated bcftools binary if it is installed outside PATH."
        )
    if shutil.which("gsutil") is None:
        raise RuntimeError("gsutil is required for Stage 1 smaller-callset extraction.")

    contigs = sorted(target_df["contig"].dropna().astype(str).unique())
    remote_vcfs = _discover_remote_vcf_shards(config, vcf_root, contigs)
    if not remote_vcfs:
        raise RuntimeError(f"No VCF shards found under {vcf_root} for contigs {', '.join(contigs)}")

    temp_root = tempfile.mkdtemp(prefix=f"stage1-vcf-{callset_name}-")
    try:
        sample_path = os.path.join(temp_root, "samples.txt")
        regions_path = os.path.join(temp_root, "targets.tsv")
        with open(sample_path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(sample_ids) + "\n")
        regions = target_df.loc[:, ["contig", "position"]].copy()
        regions["end"] = regions["position"]
        regions.to_csv(
            regions_path,
            sep="\t",
            index=False,
            header=False,
        )

        to_copy: list[str] = []
        for vcf in remote_vcfs:
            to_copy.append(vcf)
            to_copy.append(f"{vcf}.tbi")
        _copy_remote_objects(config, to_copy, temp_root)

        frames: list[pd.DataFrame] = []
        for remote_vcf in remote_vcfs:
            local_vcf = os.path.join(temp_root, os.path.basename(remote_vcf))
            subset_vcf = os.path.join(temp_root, f"{os.path.basename(remote_vcf)}.subset.vcf.gz")
            split_vcf = os.path.join(temp_root, f"{os.path.basename(remote_vcf)}.split.vcf.gz")
            _run_checked(
                [
                    bcftools,
                    "view",
                    "-S",
                    sample_path,
                    "-R",
                    regions_path,
                    "-Oz",
                    "-o",
                    subset_vcf,
                    local_vcf,
                ]
            )
            _run_checked(
                [
                    bcftools,
                    "norm",
                    "-m",
                    "-any",
                    "-Oz",
                    "-o",
                    split_vcf,
                    subset_vcf,
                ]
            )
            query = _run_checked(
                [
                    bcftools,
                    "query",
                    "-f",
                    "[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]",
                    split_vcf,
                ]
            ).stdout
            if not query.strip():
                continue
            frame = pd.read_csv(
                io.StringIO(query),
                sep="\t",
                names=["contig", "position", "ref", "alt", "person_id", "gt"],
                dtype={"contig": str, "position": int, "ref": str, "alt": str, "person_id": str, "gt": str},
            )
            frame["dosage"] = frame["gt"].map(_gt_to_dosage)
            frame = frame[frame["dosage"] > 0].copy()
            if frame.empty:
                continue
            frame["contig"] = frame["contig"].map(_normalize_contig)
            frame["variant_id"] = (
                frame["contig"].str.replace("chr", "", regex=False)
                + "-"
                + frame["position"].astype(str)
                + "-"
                + frame["ref"]
                + "-"
                + frame["alt"]
            )
            annotated = frame.merge(
                target_df,
                on=["contig", "position", "ref", "alt", "variant_id"],
                how="inner",
                suffixes=("", "_target"),
            )
            if annotated.empty:
                continue
            annotated["callset"] = callset_name
            frames.append(
                annotated[
                    [
                        "person_id",
                        "variant_id",
                        "gene",
                        "label",
                        "rsid",
                        "source",
                        "evidence_tier",
                        "exact_test_model",
                        "dosage",
                        "callset",
                    ]
                ].copy()
            )
        if not frames:
            return pd.DataFrame()
        return pd.concat(frames, ignore_index=True)
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _collapse_stage1_rows(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "person_id",
                "variant_id",
                "gene",
                "label",
                "rsid",
                "source",
                "evidence_tier",
                "exact_test_model",
                "dosage",
                "callset",
            ]
        )

    collapsed = (
        frame.sort_values(["person_id", "variant_id", "dosage"], ascending=[True, True, False])
        .groupby(["person_id", "variant_id"], as_index=False)
        .agg(
            gene=("gene", "first"),
            label=("label", "first"),
            rsid=("rsid", "first"),
            source=("source", "first"),
            evidence_tier=("evidence_tier", "first"),
            exact_test_model=("exact_test_model", "first"),
            dosage=("dosage", "max"),
            callset=("callset", lambda values: ",".join(sorted(set(str(value) for value in values if value)))),
        )
    )
    return collapsed


def prepare_stage1_variant_table(
    config: ProjectConfig,
    matched_df: pd.DataFrame,
) -> pd.DataFrame:
    config = apply_runtime_defaults(config)
    stage = config.analysis.stage1
    if stage is None:
        return pd.DataFrame()

    target_df = _panel_targets_frame(config.panel.a_priori_variants)
    sample_ids = _matched_person_ids(matched_df)
    output_path = stage.variant_table
    ensure_parent_dir(output_path)
    print(
        f"Preparing Stage 1 extraction for {len(sample_ids)} matched participants and {len(target_df)} target variants.",
        flush=True,
    )

    if target_df.empty:
        empty = _collapse_stage1_rows(pd.DataFrame())
        write_dataframe(empty, output_path)
        print(f"No target variants configured; wrote empty Stage 1 table to {output_path}", flush=True)
        return empty

    frames: list[pd.DataFrame] = []
    attempted: dict[str, str] = {}
    failures: dict[str, str] = {}
    vcf_paths = _vcf_callset_paths(config)
    for callset_name, vcf_root in vcf_paths.items():
        attempted[callset_name] = vcf_root
        print(f"Reading {callset_name} callset from {vcf_root}", flush=True)
        try:
            frame = _extract_from_vcf_callset(
                config=config,
                callset_name=callset_name,
                vcf_root=vcf_root,
                target_df=target_df,
                sample_ids=sample_ids,
            )
        except Exception as exc:
            failures[callset_name] = f"{type(exc).__name__}: {exc}"
            print(f"{callset_name} failed: {type(exc).__name__}: {exc}", flush=True)
            continue
        if not frame.empty:
            frames.append(frame)
            print(f"{callset_name} yielded {len(frame)} non-reference genotype rows.", flush=True)
        else:
            print(f"{callset_name} yielded 0 non-reference genotype rows.", flush=True)

        if failures and not frames:
            details = "\n".join(f"- {name}: {message}" for name, message in failures.items())
            raise RuntimeError(
                "Stage 1 extraction requires the AoU smaller-callset VCF workflow and a tool-ready environment.\n"
                "Install or use an environment with `bcftools` and `gsutil`, then retry.\n"
                f"{details}"
            )
        if failures:
            print("Proceeding with partial Stage 1 extraction after callset failures:", flush=True)
            for name, message in failures.items():
                print(f"- {name}: {message}", flush=True)

    combined = _collapse_stage1_rows(pd.concat(frames, ignore_index=True) if frames else pd.DataFrame())
    write_dataframe(combined, output_path)

    metadata_path = str(Path(output_path).with_suffix(Path(output_path).suffix + ".meta.json"))
    write_json(
        {
            "requested_variants": int(target_df["variant_id"].nunique()),
            "matched_people": int(len(sample_ids)),
            "variants_with_hits": int(combined["variant_id"].nunique()) if not combined.empty else 0,
            "rows_written": int(len(combined)),
            "callsets": attempted,
            "failures": failures,
        },
        metadata_path,
    )
    print(f"Wrote {len(combined)} collapsed Stage 1 genotype rows to {output_path}", flush=True)
    return combined


__all__ = [
    "prepare_stage1_variant_table",
    "_callset_root",
    "_bcftools_binary",
    "_collapse_stage1_rows",
    "_gt_to_dosage",
    "_normalize_contig",
    "_panel_targets_frame",
    "_vcf_root",
]
