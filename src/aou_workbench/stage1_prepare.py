"""Stage 1 input preparation from AoU smaller genomics callsets."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import ProjectConfig, TargetVariant
from .hail_utils import ensure_hail
from .io_utils import ensure_parent_dir, write_dataframe, write_json
from .preflight import apply_runtime_defaults


HAIL_REFERENCE = "GRCh38"


def _callset_root(path: str) -> str:
    clean = path.rstrip("/")
    for suffix in ("/splitMT/hail.mt", "/multiMT/hail.mt", "/hail.mt"):
        if clean.endswith(suffix):
            return clean[: -len(suffix)]
    return clean


def _split_mt_path(path: str) -> str:
    return f"{_callset_root(path)}/splitMT/hail.mt"


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


def _callset_paths(config: ProjectConfig) -> dict[str, str]:
    return {
        "acaf": _split_mt_path(config.workbench.acaf_mt_path),
        "exome": _split_mt_path(config.workbench.exome_mt_path),
        "clinvar": _split_mt_path(config.workbench.clinvar_mt_path),
    }


def _extract_from_callset(
    hl,
    *,
    callset_name: str,
    mt_path: str,
    target_df: pd.DataFrame,
    sample_ids: list[str],
) -> pd.DataFrame:
    if target_df.empty or not sample_ids:
        return pd.DataFrame()

    mt = hl.read_matrix_table(mt_path)
    sample_set = hl.literal(set(sample_ids))
    mt = mt.filter_cols(sample_set.contains(hl.str(mt.s)))

    target_ht = hl.Table.from_pandas(target_df)
    target_ht = target_ht.annotate(
        locus=hl.locus(target_ht.contig, hl.int32(target_ht.position), reference_genome=HAIL_REFERENCE),
        alleles=[target_ht.ref, target_ht.alt],
    ).key_by("locus", "alleles")

    mt = mt.filter_rows(hl.is_defined(target_ht[mt.locus, mt.alleles]))
    mt = mt.annotate_rows(target=target_ht[mt.locus, mt.alleles])
    mt = mt.filter_rows(hl.is_defined(mt.target))

    entries = mt.entries()
    entries = entries.filter(hl.is_defined(entries.GT) & (entries.GT.n_alt_alleles() > 0))
    frame = entries.select(
        person_id=hl.str(entries.s),
        variant_id=entries.target.variant_id,
        gene=entries.target.gene,
        label=entries.target.label,
        rsid=entries.target.rsid,
        source=entries.target.source,
        evidence_tier=entries.target.evidence_tier,
        exact_test_model=entries.target.exact_test_model,
        dosage=hl.float64(entries.GT.n_alt_alleles()),
        callset=hl.str(callset_name),
    ).to_pandas()
    if frame.empty:
        return frame
    frame["person_id"] = frame["person_id"].astype(str)
    frame["variant_id"] = frame["variant_id"].astype(str)
    return frame


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

    hl = ensure_hail(
        HAIL_REFERENCE,
        requester_pays_project=config.workbench.requester_pays_project,
        requester_pays_buckets=config.workbench.requester_pays_buckets,
    )

    frames: list[pd.DataFrame] = []
    attempted: dict[str, str] = {}
    failures: dict[str, str] = {}
    for callset_name, mt_path in _callset_paths(config).items():
        attempted[callset_name] = mt_path
        print(f"Reading {callset_name} callset from {mt_path}", flush=True)
        try:
            frame = _extract_from_callset(
                hl,
                callset_name=callset_name,
                mt_path=mt_path,
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

    if not frames and failures:
        details = "\n".join(f"- {name}: {message}" for name, message in failures.items())
        raise RuntimeError(f"Stage 1 extraction failed for all callsets.\n{details}")

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
    "_collapse_stage1_rows",
    "_panel_targets_frame",
    "_split_mt_path",
]
