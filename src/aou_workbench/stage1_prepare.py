"""Stage 1 input preparation from direct AoU WGS VDS extraction."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import ProjectConfig, TargetVariant
from .hail_utils import ensure_hail
from .io_utils import ensure_parent_dir, write_dataframe, write_json
from .preflight import apply_runtime_defaults


def stage1_sample_manifest_path(variant_table_path: str) -> str:
    path = Path(variant_table_path)
    if path.suffix:
        return str(path.with_suffix(f".samples{path.suffix}"))
    return f"{variant_table_path}.samples.tsv"


def _panel_targets_frame(targets: tuple[TargetVariant, ...]) -> pd.DataFrame:
    rows = [
        {
            "label": target.label,
            "gene": target.gene,
            "rsid": target.rsid or "",
            "contig": target.contig,
            "position": int(target.position),
            "ref": target.ref,
            "alt": target.alt,
            "source": target.source or "",
            "evidence_tier": target.evidence_tier or "",
            "exact_test_model": target.exact_test_model,
            "variant_id": target.variant_id,
        }
        for target in targets
    ]
    return pd.DataFrame(rows)


def _analysis_person_ids(sample_df: pd.DataFrame) -> list[str]:
    if "person_id" not in sample_df.columns:
        raise ValueError("Analysis cohort is missing the required 'person_id' column.")
    ids = sample_df["person_id"].astype(str).dropna().drop_duplicates().sort_values()
    return ids.tolist()


def _matched_person_ids(matched_df: pd.DataFrame) -> list[str]:
    return _analysis_person_ids(matched_df)


def _target_interval_strings(target_df: pd.DataFrame, pad_bp: int = 5) -> list[str]:
    intervals = {
        f"{row.contig}:{max(1, int(row.position) - pad_bp)}-{int(row.position) + pad_bp}"
        for row in target_df.itertuples(index=False)
    }
    return sorted(intervals)


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


def _target_hail_table(hl: object, target_df: pd.DataFrame):
    frame = target_df.copy()
    for column in ("label", "gene", "rsid", "contig", "ref", "alt", "source", "evidence_tier", "exact_test_model", "variant_id"):
        frame[column] = frame[column].fillna("").astype(str)
    frame["position"] = frame["position"].astype(int)
    target_ht = hl.Table.from_pandas(frame)
    target_ht = target_ht.annotate(
        locus=hl.locus(target_ht.contig, hl.int32(target_ht.position), "GRCh38"),
        alleles=[target_ht.ref, target_ht.alt],
    )
    return target_ht.key_by("locus", "alleles").select(
        "label",
        "gene",
        "rsid",
        "source",
        "evidence_tier",
        "exact_test_model",
        "variant_id",
    )


def _extract_from_vds(
    *,
    config: ProjectConfig,
    target_df: pd.DataFrame,
    sample_ids: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    if target_df.empty or not sample_ids:
        return pd.DataFrame(), pd.DataFrame(columns=["person_id"]), {"rows": 0, "cols": 0, "entries": 0}

    hl = ensure_hail(
        requester_pays_project=config.workbench.requester_pays_project,
        requester_pays_buckets=config.workbench.requester_pays_buckets,
    )
    target_ht = _target_hail_table(hl, target_df)
    intervals = [hl.parse_locus_interval(value, reference_genome="GRCh38") for value in _target_interval_strings(target_df)]

    print(f"Reading direct WGS VDS from {config.workbench.wgs_vds_path}", flush=True)
    vds = hl.vds.read_vds(config.workbench.wgs_vds_path)
    vds = hl.vds.filter_intervals(vds, intervals)
    vds = hl.vds.split_multi(vds)

    mt = vds.variant_data
    sample_set = hl.literal(set(sample_ids))
    mt = mt.filter_cols(sample_set.contains(mt.s))
    mt = mt.filter_rows(hl.is_defined(target_ht[mt.locus, mt.alleles]))
    mt = mt.annotate_rows(target=target_ht[mt.locus, mt.alleles])

    row_count = mt.count_rows()
    col_count = mt.count_cols()
    print(f"Stage 1 VDS slice contains {row_count} exact variant rows across {col_count} selected samples.", flush=True)
    col_ht = mt.cols()
    col_ht = col_ht.key_by()
    col_ht = col_ht.select(person_id=col_ht.s)
    present_samples = pd.DataFrame({"person_id": [str(row.person_id) for row in col_ht.collect()]})

    mt = mt.filter_entries(hl.is_defined(mt.GT) & mt.GT.is_non_ref())
    entries = mt.entries()
    entries = entries.key_by()
    entries = entries.select(
        person_id=entries.s,
        variant_id=entries.target.variant_id,
        gene=entries.target.gene,
        label=entries.target.label,
        rsid=entries.target.rsid,
        source=entries.target.source,
        evidence_tier=entries.target.evidence_tier,
        exact_test_model=entries.target.exact_test_model,
        dosage=hl.float64(entries.GT.n_alt_alleles()),
        callset=hl.str("wgs_vds"),
    )

    rows = entries.collect()
    if not rows:
        return pd.DataFrame(), present_samples, {"rows": row_count, "cols": col_count, "entries": 0}

    frame = pd.DataFrame(
        [
            {
                "person_id": row.person_id,
                "variant_id": row.variant_id,
                "gene": row.gene,
                "label": row.label,
                "rsid": row.rsid,
                "source": row.source,
                "evidence_tier": row.evidence_tier,
                "exact_test_model": row.exact_test_model,
                "dosage": float(row.dosage),
                "callset": row.callset,
            }
            for row in rows
        ]
    )
    return frame, present_samples, {"rows": row_count, "cols": col_count, "entries": len(frame)}


def prepare_stage1_variant_table(
    config: ProjectConfig,
    sample_df: pd.DataFrame,
) -> pd.DataFrame:
    config = apply_runtime_defaults(config)
    stage = config.analysis.stage1
    if stage is None:
        return pd.DataFrame()

    target_df = _panel_targets_frame(config.panel.a_priori_variants)
    sample_ids = _analysis_person_ids(sample_df)
    output_path = stage.variant_table
    sample_manifest_path = stage1_sample_manifest_path(output_path)
    ensure_parent_dir(output_path)
    ensure_parent_dir(sample_manifest_path)
    print(
        f"Preparing Stage 1 extraction for {len(sample_ids)} analysis participants and {len(target_df)} target variants.",
        flush=True,
    )

    if target_df.empty:
        empty = _collapse_stage1_rows(pd.DataFrame())
        write_dataframe(empty, output_path)
        write_dataframe(pd.DataFrame(columns=["person_id"]), sample_manifest_path)
        print(f"No target variants configured; wrote empty Stage 1 table to {output_path}", flush=True)
        return empty

    raw, present_samples, stats = _extract_from_vds(config=config, target_df=target_df, sample_ids=sample_ids)
    combined = _collapse_stage1_rows(raw)
    write_dataframe(combined, output_path)
    write_dataframe(present_samples, sample_manifest_path)

    metadata_path = str(Path(output_path).with_suffix(Path(output_path).suffix + ".meta.json"))
    write_json(
        {
            "requested_variants": int(target_df["variant_id"].nunique()),
            "analysis_people_requested": int(len(sample_ids)),
            "analysis_people_in_vds_slice": int(stats["cols"]),
            "analysis_people_manifest": sample_manifest_path,
            "exact_variant_rows_in_vds_slice": int(stats["rows"]),
            "non_reference_entries": int(stats["entries"]),
            "variants_with_hits": int(combined["variant_id"].nunique()) if not combined.empty else 0,
            "rows_written": int(len(combined)),
            "wgs_vds_path": config.workbench.wgs_vds_path,
        },
        metadata_path,
    )
    print(f"Wrote {len(combined)} collapsed Stage 1 genotype rows to {output_path}", flush=True)
    return combined


__all__ = [
    "prepare_stage1_variant_table",
    "_collapse_stage1_rows",
    "_analysis_person_ids",
    "_matched_person_ids",
    "_panel_targets_frame",
    "stage1_sample_manifest_path",
    "_target_interval_strings",
]
