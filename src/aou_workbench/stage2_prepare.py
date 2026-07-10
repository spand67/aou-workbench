"""Stage 2 input preparation from VAT annotations plus targeted VDS genotypes."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from .annotations import annotate_variant_masks
from .config import ProjectConfig
from .hail_utils import ensure_hail
from .io_utils import ensure_parent_dir, stable_hash, slugify, write_dataframe, write_json
from .paths import join_path
from .preflight import apply_runtime_defaults
from .stage1_prepare import _analysis_person_ids


def stage2_sample_manifest_path(variant_table_path: str) -> str:
    path = Path(variant_table_path)
    if path.suffix:
        return str(path.with_suffix(f".samples{path.suffix}"))
    return f"{variant_table_path}.samples.tsv"


def _struct_fields(expr: object) -> set[str]:
    dtype = getattr(expr, "dtype", None)
    fields = getattr(dtype, "fields", ())
    return {str(field) for field in fields}


def _pick_expr(struct_expr: object, candidates: tuple[str, ...]):
    fields = _struct_fields(struct_expr)
    for name in candidates:
        if name in fields:
            return struct_expr[name]
    return None


VAT_COLUMNS = (
    "vid",
    "contig",
    "position",
    "ref_allele",
    "alt_allele",
    "gene_symbol",
    "consequence",
    "revel",
    "gnomad_max_af",
    "gvs_max_af",
    "gvs_all_af",
    "clinvar_classification",
    "aa_change",
    "dbsnp_rsid",
)


def _join_unique(values: pd.Series) -> str:
    unique = sorted({str(value).strip() for value in values.dropna().tolist() if str(value).strip()})
    return ";".join(unique)


def _max_numeric(values: pd.Series) -> float | None:
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().any():
        return float(numeric.max())
    return None


def _variant_id_from_vat(row: pd.Series) -> str:
    vid = str(row.get("vid", "")).strip()
    if vid:
        return vid.replace("chr", "")
    contig = str(row["contig"]).replace("chr", "")
    return f"{contig}-{int(row['position'])}-{row['ref_allele']}-{row['alt_allele']}"


def _stage2_candidate_cache_path(
    config: ProjectConfig,
    *,
    variant_table_path: str,
    genes: list[str],
    max_af: float,
    revel_min: float,
    plof_terms: tuple[str, ...],
    clinvar_plp_terms: tuple[str, ...],
) -> str:
    cache_id = stable_hash(
        {
            "vat_path": config.workbench.vat_path,
            "genes": genes,
            "max_af": max_af,
            "revel_min": revel_min,
            "plof_terms": plof_terms,
            "clinvar_plp_terms": clinvar_plp_terms,
        }
    )
    bucket = config.workbench.workspace_bucket or os.getenv("WORKSPACE_BUCKET")
    if bucket:
        return join_path(
            bucket,
            "aou-workbench",
            slugify(config.analysis.analysis_name),
            "stage2",
            f"vat-candidates-{cache_id}.ht",
        )
    if variant_table_path.startswith("gs://"):
        return join_path(variant_table_path.rsplit("/", 1)[0], f"vat-candidates-{cache_id}.ht")
    return f"{Path(variant_table_path).with_suffix('')}.vat_candidates_{cache_id}.ht"


def _hail_path_exists(hl: object, path: str) -> bool:
    try:
        return bool(hl.current_backend().fs.exists(path))
    except Exception:
        return False


def _hail_table_to_pandas(table: object) -> pd.DataFrame:
    rows = table.key_by().collect()
    return pd.DataFrame([dict(row) for row in rows])


def _read_vds_intervals(hl: object, path: str, intervals: list[object]):
    try:
        return hl.vds.read_vds(path, intervals=intervals)
    except TypeError:
        vds = hl.vds.read_vds(path)
        return hl.vds.filter_intervals(vds, intervals)


def _collapse_vat_annotations(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "variant_id",
                "contig",
                "position",
                "ref_allele",
                "alt_allele",
                "gene",
                "clinvar_significance",
                "consequence",
                "revel",
                "max_af",
                "aa_change",
                "dbsnp_rsid",
            ]
        )
    out = frame.copy()
    for column in VAT_COLUMNS:
        if column not in out.columns:
            out[column] = ""
    out = out[out["gene_symbol"].fillna("").astype(str).str.strip() != ""].copy()
    if out.empty:
        return _collapse_vat_annotations(pd.DataFrame())

    out["position"] = pd.to_numeric(out["position"], errors="coerce")
    out = out.dropna(subset=["position", "contig", "ref_allele", "alt_allele"])
    out["position"] = out["position"].astype(int)
    out["contig"] = out["contig"].astype(str)
    out["ref_allele"] = out["ref_allele"].astype(str)
    out["alt_allele"] = out["alt_allele"].astype(str)
    out["variant_id"] = out.apply(_variant_id_from_vat, axis=1)
    out["gene"] = out["gene_symbol"].astype(str)
    out["clinvar_significance"] = out["clinvar_classification"].fillna("").astype(str)
    out["max_af_input"] = out[["gnomad_max_af", "gvs_max_af", "gvs_all_af"]].apply(_max_numeric, axis=1)

    grouped = (
        out.groupby(["variant_id", "gene"], as_index=False)
        .agg(
            contig=("contig", "first"),
            position=("position", "first"),
            ref_allele=("ref_allele", "first"),
            alt_allele=("alt_allele", "first"),
            clinvar_significance=("clinvar_significance", _join_unique),
            consequence=("consequence", _join_unique),
            revel=("revel", _max_numeric),
            max_af=("max_af_input", _max_numeric),
            aa_change=("aa_change", _join_unique),
            dbsnp_rsid=("dbsnp_rsid", _join_unique),
        )
        .sort_values(["gene", "variant_id"])
        .reset_index(drop=True)
    )
    return grouped


def _vat_candidate_annotations(
    *,
    config: ProjectConfig,
    genes: list[str],
    max_af: float,
    revel_min: float,
    plof_terms: tuple[str, ...],
    clinvar_plp_terms: tuple[str, ...],
    cache_path: str | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    if not genes:
        return pd.DataFrame(), {"vat_gene_rows": 0, "vat_variant_gene_rows": 0, "candidate_variant_gene_rows": 0}

    hl = ensure_hail(
        requester_pays_project=config.workbench.requester_pays_project,
        requester_pays_buckets=config.workbench.requester_pays_buckets,
    )
    if cache_path and _hail_path_exists(hl, cache_path):
        print(f"Using cached Stage 2 VAT candidates from {cache_path}", flush=True)
        candidates = _hail_table_to_pandas(hl.read_table(cache_path))
        return candidates, {
            "vat_gene_rows": None,
            "vat_variant_gene_rows": None,
            "candidate_variant_gene_rows": int(len(candidates)),
            "vat_candidate_cache_path": cache_path,
            "vat_candidate_cache_used": True,
        }

    print(f"Reading VAT annotations from {config.workbench.vat_path}", flush=True)
    vat = hl.import_table(
        config.workbench.vat_path,
        impute=False,
        force_bgz=True,
        missing="",
    )
    fields = _struct_fields(vat.row_value)
    required = {"contig", "position", "ref_allele", "alt_allele", "gene_symbol"}
    missing_required = sorted(required - fields)
    if missing_required:
        raise RuntimeError(
            "Stage 2 could not find required VAT fields "
            f"{missing_required}. Available VAT fields: {sorted(fields)}"
        )

    gene_set = hl.literal(set(genes))
    vat = vat.filter(gene_set.contains(vat["gene_symbol"]))
    select_exprs = {column: vat[column] for column in VAT_COLUMNS if column in fields}
    vat = vat.select(**select_exprs)

    rows = vat.collect()
    raw = pd.DataFrame([dict(row) for row in rows])
    collapsed = _collapse_vat_annotations(raw)
    annotated = annotate_variant_masks(
        collapsed,
        clinvar_column="clinvar_significance",
        consequence_column="consequence",
        revel_column="revel",
        af_column="max_af",
        max_af=max_af,
        revel_min=revel_min,
        plof_terms=plof_terms,
        clinvar_plp_terms=clinvar_plp_terms,
    )
    candidates = annotated[annotated["mask_primary"]].copy()
    print(
        "Stage 2 VAT candidates: "
        f"{len(raw)} transcript rows, {len(collapsed)} variant-gene rows, "
        f"{len(candidates)} primary-mask variant-gene rows.",
        flush=True,
    )
    if cache_path and not candidates.empty:
        print(f"Writing Stage 2 VAT candidate cache to {cache_path}", flush=True)
        hl.Table.from_pandas(candidates).write(cache_path, overwrite=True)
    return candidates, {
        "vat_gene_rows": int(len(raw)),
        "vat_variant_gene_rows": int(len(collapsed)),
        "candidate_variant_gene_rows": int(len(candidates)),
        "vat_candidate_cache_path": cache_path,
        "vat_candidate_cache_used": False,
    }


def _target_table_from_annotations(hl: object, annotations: pd.DataFrame):
    targets = annotations[
        ["variant_id", "contig", "position", "ref_allele", "alt_allele"]
    ].drop_duplicates("variant_id")
    target_ht = hl.Table.from_pandas(targets)
    target_ht = target_ht.annotate(
        locus=hl.locus(target_ht.contig, hl.int32(target_ht.position), "GRCh38"),
        alleles=[target_ht.ref_allele, target_ht.alt_allele],
    )
    return target_ht.key_by("locus", "alleles").select("variant_id")


def _target_intervals_from_annotations(hl: object, annotations: pd.DataFrame):
    intervals = {
        f"{row.contig}:{int(row.position)}-{int(row.position)}"
        for row in annotations[["contig", "position"]].drop_duplicates().itertuples(index=False)
    }
    return [hl.parse_locus_interval(value, reference_genome="GRCh38") for value in sorted(intervals)]


def _extract_candidate_genotypes_from_vds(
    *,
    config: ProjectConfig,
    sample_ids: list[str],
    annotations: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    if not sample_ids or annotations.empty:
        return pd.DataFrame(), pd.DataFrame(columns=["person_id"]), {"rows": 0, "cols": 0, "entries": 0}

    hl = ensure_hail(
        requester_pays_project=config.workbench.requester_pays_project,
        requester_pays_buckets=config.workbench.requester_pays_buckets,
    )
    target_ht = _target_table_from_annotations(hl, annotations)
    intervals = _target_intervals_from_annotations(hl, annotations)

    print(f"Reading direct WGS VDS from {config.workbench.wgs_vds_path}", flush=True)
    vds = _read_vds_intervals(hl, config.workbench.wgs_vds_path, intervals)
    vds = hl.vds.split_multi(vds)
    mt = vds.variant_data

    sample_set = hl.literal(set(sample_ids))
    mt = mt.filter_cols(sample_set.contains(mt.s))
    mt = mt.filter_rows(hl.is_defined(target_ht[mt.locus, mt.alleles]))
    mt = mt.annotate_rows(target=target_ht[mt.locus, mt.alleles])

    col_ht = mt.cols()
    col_ht = col_ht.key_by()
    col_ht = col_ht.select(person_id=col_ht.s)
    present_samples = pd.DataFrame({"person_id": [str(row.person_id) for row in col_ht.collect()]})
    col_count = len(present_samples)
    print(
        "Stage 2 VDS candidate slice restricted to "
        f"{col_count} selected samples; skipping separate variant row count for speed.",
        flush=True,
    )

    entry_fields = _struct_fields(mt.entry)
    if "GT" in entry_fields:
        mt = mt.annotate_entries(_stage2_dosage=hl.float64(mt.GT.n_alt_alleles()))
    elif "LGT" in entry_fields:
        mt = mt.annotate_entries(_stage2_dosage=hl.float64(mt.LGT.n_alt_alleles()))
    else:
        raise RuntimeError(f"Stage 2 VDS genotype entries have no GT or LGT field. Entry fields: {sorted(entry_fields)}")
    mt = mt.filter_entries(hl.is_defined(mt._stage2_dosage) & (mt._stage2_dosage > 0))
    entries = mt.entries().key_by()
    entries = entries.select(
        person_id=entries.s,
        variant_id=entries.target.variant_id,
        dosage=entries._stage2_dosage,
    )

    rows = entries.collect()
    if not rows:
        return pd.DataFrame(), present_samples, {"rows": None, "cols": col_count, "entries": 0}

    frame = pd.DataFrame(
        [
            {
                "person_id": row.person_id,
                "variant_id": row.variant_id,
                "dosage": float(row.dosage),
            }
            for row in rows
        ]
    )
    merged = frame.merge(annotations, on="variant_id", how="left")
    merged = merged.rename(columns={"ref_allele": "ref", "alt_allele": "alt"})
    keep_columns = [
        "person_id",
        "variant_id",
        "gene",
        "clinvar_significance",
        "consequence",
        "revel",
        "max_af",
        "dosage",
        "contig",
        "position",
        "ref",
        "alt",
        "aa_change",
        "dbsnp_rsid",
    ]
    for column in keep_columns:
        if column not in merged.columns:
            merged[column] = pd.NA
    return merged[keep_columns], present_samples, {"rows": None, "cols": col_count, "entries": len(merged)}


def prepare_stage2_variant_table(
    config: ProjectConfig,
    sample_df: pd.DataFrame,
) -> pd.DataFrame:
    config = apply_runtime_defaults(config)
    stage = config.analysis.stage2
    if stage is None:
        return pd.DataFrame()

    sample_ids = _analysis_person_ids(sample_df)
    genes = sorted(set(config.panel.genes_of_interest))
    output_path = stage.variant_table
    sample_manifest_path = stage2_sample_manifest_path(output_path)
    ensure_parent_dir(output_path)
    ensure_parent_dir(sample_manifest_path)
    print(
        f"Preparing Stage 2 VAT/VDS extraction for {len(sample_ids)} analysis participants across {len(genes)} genes.",
        flush=True,
    )

    vat_cache_path = _stage2_candidate_cache_path(
        config,
        variant_table_path=output_path,
        genes=genes,
        max_af=stage.max_af,
        revel_min=stage.revel_min,
        plof_terms=stage.plof_terms,
        clinvar_plp_terms=stage.clinvar_plp_terms,
    )
    annotations, vat_stats = _vat_candidate_annotations(
        config=config,
        genes=genes,
        max_af=stage.max_af,
        revel_min=stage.revel_min,
        plof_terms=stage.plof_terms,
        clinvar_plp_terms=stage.clinvar_plp_terms,
        cache_path=vat_cache_path,
    )
    raw, present_samples, stats = _extract_candidate_genotypes_from_vds(
        config=config,
        sample_ids=sample_ids,
        annotations=annotations,
    )
    write_dataframe(raw, output_path)
    write_dataframe(present_samples, sample_manifest_path)

    metadata_path = str(Path(output_path).with_suffix(Path(output_path).suffix + ".meta.json"))
    write_json(
        {
            "analysis_people_requested": int(len(sample_ids)),
            "analysis_people_in_vds_slice": int(stats["cols"]),
            "analysis_people_manifest": sample_manifest_path,
            "genes_requested": genes,
            **vat_stats,
            "exact_variant_rows_in_vds_slice": (
                int(stats["rows"]) if stats.get("rows") is not None else None
            ),
            "non_reference_entries": int(stats["entries"]),
            "variants_with_hits": int(raw["variant_id"].nunique()) if not raw.empty else 0,
            "rows_written": int(len(raw)),
            "vat_path": config.workbench.vat_path,
            "wgs_vds_path": config.workbench.wgs_vds_path,
        },
        metadata_path,
    )
    print(f"Wrote {len(raw)} Stage 2 VAT-classified VDS genotype rows to {output_path}", flush=True)
    return raw


def prepare_stage2_vat_candidate_cache(config: ProjectConfig) -> tuple[pd.DataFrame, dict[str, object]]:
    config = apply_runtime_defaults(config)
    stage = config.analysis.stage2
    if stage is None:
        return pd.DataFrame(), {"stage2_configured": False}

    genes = sorted(set(config.panel.genes_of_interest))
    cache_path = _stage2_candidate_cache_path(
        config,
        variant_table_path=stage.variant_table,
        genes=genes,
        max_af=stage.max_af,
        revel_min=stage.revel_min,
        plof_terms=stage.plof_terms,
        clinvar_plp_terms=stage.clinvar_plp_terms,
    )
    print(
        f"Preparing Stage 2 VAT candidate cache for {len(genes)} genes "
        f"(max_af={stage.max_af}, revel_min={stage.revel_min}).",
        flush=True,
    )
    candidates, stats = _vat_candidate_annotations(
        config=config,
        genes=genes,
        max_af=stage.max_af,
        revel_min=stage.revel_min,
        plof_terms=stage.plof_terms,
        clinvar_plp_terms=stage.clinvar_plp_terms,
        cache_path=cache_path,
    )
    stats = {
        **stats,
        "stage2_configured": True,
        "genes_requested": genes,
        "rows_available": int(len(candidates)),
    }
    return candidates, stats


__all__ = [
    "prepare_stage2_variant_table",
    "prepare_stage2_vat_candidate_cache",
    "stage2_sample_manifest_path",
    "_collapse_vat_annotations",
    "_vat_candidate_annotations",
]
