"""Annotation helpers shared across stage analyses."""

from __future__ import annotations

import pandas as pd


def annotate_variant_masks(
    df: pd.DataFrame,
    *,
    clinvar_column: str,
    consequence_column: str,
    revel_column: str,
    af_column: str,
    max_af: float,
    revel_min: float,
    plof_terms: tuple[str, ...],
    clinvar_plp_terms: tuple[str, ...],
) -> pd.DataFrame:
    out = df.copy()
    clinvar = out.get(clinvar_column, pd.Series([""] * len(out), index=out.index)).astype(str).str.lower()
    consequence = out.get(consequence_column, pd.Series([""] * len(out), index=out.index)).astype(str).str.lower()
    revel = pd.to_numeric(out.get(revel_column, pd.Series([None] * len(out), index=out.index)), errors="coerce")
    af = pd.to_numeric(out.get(af_column, pd.Series([None] * len(out), index=out.index)), errors="coerce")

    out["is_clinvar_plp"] = False
    for term in clinvar_plp_terms:
        out["is_clinvar_plp"] |= clinvar.str.contains(str(term).lower(), regex=False, na=False)

    out["is_plof"] = False
    for term in plof_terms:
        out["is_plof"] |= consequence.str.contains(str(term).lower(), regex=False, na=False)

    out["is_revel_high"] = revel >= revel_min
    out["passes_af"] = af.isna() | (af <= max_af)
    out["mask_primary"] = out["passes_af"] & (
        out["is_clinvar_plp"] | out["is_plof"] | out["is_revel_high"]
    )
    out["mask_plof"] = out["passes_af"] & out["is_plof"]
    out["mask_clinvar_plp"] = out["passes_af"] & out["is_clinvar_plp"]
    out["mask_plof_or_revel"] = out["passes_af"] & (out["is_plof"] | out["is_revel_high"])
    return out


__all__ = ["annotate_variant_masks"]
