"""Stage 2 input preparation from the AoU ClinVar smaller callset."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import ProjectConfig
from .hail_utils import ensure_hail
from .io_utils import ensure_parent_dir, write_dataframe, write_json
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


def _gene_symbol_expr(hl: object, mt: object):
    direct = _pick_expr(mt.row_value, ("gene", "gene_symbol", "symbol"))
    if direct is not None:
        return hl.str(direct)
    info = getattr(mt, "info", None)
    if info is None:
        return None
    geneinfo = _pick_expr(info, ("GENEINFO", "geneinfo", "gene_info"))
    if geneinfo is not None:
        text = hl.str(geneinfo)
        return hl.if_else(
            hl.len(text) > 0,
            hl.split(hl.split(text, ",")[0], ":")[0],
            hl.missing(hl.tstr),
        )
    nested = _pick_expr(info, ("GENE", "gene", "gene_symbol", "symbol"))
    if nested is not None:
        return hl.str(nested)
    return None


def _clinvar_significance_expr(hl: object, mt: object):
    direct = _pick_expr(mt.row_value, ("clinvar_significance", "clnsig", "CLNSIG"))
    if direct is not None:
        return hl.str(direct)
    info = getattr(mt, "info", None)
    if info is None:
        return hl.missing(hl.tstr)
    nested = _pick_expr(info, ("CLNSIG", "clnsig", "clinvar_significance"))
    if nested is not None:
        return hl.str(nested)
    return hl.missing(hl.tstr)


def _consequence_expr(hl: object, mt: object):
    direct = _pick_expr(mt.row_value, ("consequence", "most_severe_consequence", "csq"))
    if direct is not None:
        return hl.str(direct)
    info = getattr(mt, "info", None)
    if info is None:
        return hl.missing(hl.tstr)
    nested = _pick_expr(info, ("CSQ", "csq", "consequence", "most_severe_consequence"))
    if nested is not None:
        return hl.str(nested)
    return hl.missing(hl.tstr)


def _max_af_expr(hl: object, mt: object):
    direct = _pick_expr(mt.row_value, ("max_af", "maf", "AF"))
    if direct is not None:
        return hl.float64(direct)
    info = getattr(mt, "info", None)
    if info is None:
        return hl.missing(hl.tfloat64)
    nested = _pick_expr(info, ("AF", "af", "max_af", "maf"))
    if nested is not None:
        return hl.float64(nested)
    return hl.missing(hl.tfloat64)


def _extract_from_clinvar_mt(
    *,
    config: ProjectConfig,
    sample_ids: list[str],
    genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    if not sample_ids or not genes:
        return pd.DataFrame(), pd.DataFrame(columns=["person_id"]), {"rows": 0, "cols": 0, "entries": 0}

    hl = ensure_hail(
        requester_pays_project=config.workbench.requester_pays_project,
        requester_pays_buckets=config.workbench.requester_pays_buckets,
    )
    print(f"Reading ClinVar MT from {config.workbench.clinvar_mt_path}", flush=True)
    mt = hl.read_matrix_table(config.workbench.clinvar_mt_path)
    mt = hl.split_multi_hts(mt)
    sample_set = hl.literal(set(sample_ids))
    mt = mt.filter_cols(sample_set.contains(mt.s))
    mt = mt.annotate_rows(
        variant_id=(
            hl.str(mt.locus.contig).replace("chr", "")
            + "-"
            + hl.str(mt.locus.position)
            + "-"
            + mt.alleles[0]
            + "-"
            + mt.alleles[1]
        )
    )
    gene_expr = _gene_symbol_expr(hl, mt)
    if gene_expr is None:
        row_fields = sorted(_struct_fields(mt.row_value))
        info_fields = sorted(_struct_fields(getattr(mt, "info", None))) if hasattr(mt, "info") else []
        raise RuntimeError(
            "Stage 2 could not find a usable gene annotation directly in the ClinVar MT. "
            f"Row fields: {row_fields}; info fields: {info_fields}"
        )
    mt = mt.annotate_rows(
        gene=gene_expr,
        clinvar_significance=_clinvar_significance_expr(hl, mt),
        consequence=_consequence_expr(hl, mt),
        max_af=_max_af_expr(hl, mt),
    )
    gene_set = hl.literal(set(genes))
    mt = mt.filter_rows(hl.is_defined(mt.gene) & gene_set.contains(mt.gene))

    row_count = mt.count_rows()
    col_count = mt.count_cols()
    print(f"Stage 2 ClinVar slice contains {row_count} annotated variant rows across {col_count} selected samples.", flush=True)

    col_ht = mt.cols()
    col_ht = col_ht.key_by()
    col_ht = col_ht.select(person_id=col_ht.s)
    present_samples = pd.DataFrame({"person_id": [str(row.person_id) for row in col_ht.collect()]})

    mt = mt.filter_entries(hl.is_defined(mt.GT) & mt.GT.is_non_ref())
    entries = mt.entries().key_by()
    entries = entries.select(
        person_id=entries.s,
        variant_id=entries.variant_id,
        gene=entries.gene,
        clinvar_significance=entries.clinvar_significance,
        consequence=entries.consequence,
        max_af=entries.max_af,
        dosage=hl.float64(entries.GT.n_alt_alleles()),
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
                "clinvar_significance": row.clinvar_significance,
                "consequence": row.consequence,
                "max_af": row.max_af,
                "dosage": float(row.dosage),
            }
            for row in rows
        ]
    )
    return frame, present_samples, {"rows": row_count, "cols": col_count, "entries": len(frame)}


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
        f"Preparing Stage 2 ClinVar extraction for {len(sample_ids)} analysis participants across {len(genes)} genes.",
        flush=True,
    )

    raw, present_samples, stats = _extract_from_clinvar_mt(config=config, sample_ids=sample_ids, genes=genes)
    write_dataframe(raw, output_path)
    write_dataframe(present_samples, sample_manifest_path)

    metadata_path = str(Path(output_path).with_suffix(Path(output_path).suffix + ".meta.json"))
    write_json(
        {
            "analysis_people_requested": int(len(sample_ids)),
            "analysis_people_in_clinvar_mt": int(stats["cols"]),
            "analysis_people_manifest": sample_manifest_path,
            "genes_requested": genes,
            "annotated_variant_rows_in_slice": int(stats["rows"]),
            "non_reference_entries": int(stats["entries"]),
            "variants_with_hits": int(raw["variant_id"].nunique()) if not raw.empty else 0,
            "rows_written": int(len(raw)),
            "clinvar_mt_path": config.workbench.clinvar_mt_path,
        },
        metadata_path,
    )
    print(f"Wrote {len(raw)} Stage 2 ClinVar genotype rows to {output_path}", flush=True)
    return raw


__all__ = ["prepare_stage2_variant_table", "stage2_sample_manifest_path"]
