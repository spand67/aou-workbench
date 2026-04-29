"""Helpers for restricting cohorts to WGS-present and optional unrelated samples."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import ProjectConfig
from .io_utils import read_table
from .stage1_prepare import stage1_sample_manifest_path

_ID_CANDIDATES = ("person_id", "research_id", "IID", "FID", "s", "sample_id")


def _resolve_id_column(frame: pd.DataFrame) -> str:
    for column in _ID_CANDIDATES:
        if column in frame.columns:
            return column
    raise RuntimeError(f"Could not find a sample identifier column in {list(frame.columns)}.")


def _load_id_set(path: str) -> set[str]:
    frame = read_table(path)
    if frame.empty:
        raise RuntimeError(f"Sample restriction file at {path} is empty.")
    id_column = _resolve_id_column(frame)
    ids = set(frame[id_column].astype(str).dropna())
    if not ids:
        raise RuntimeError(f"Sample restriction file at {path} does not contain any usable IDs.")
    return ids


def wgs_manifest_path(config: ProjectConfig) -> str:
    stage = config.analysis.stage1
    if stage is None:
        raise RuntimeError("Stage 1 must be configured to derive the WGS-present manifest.")
    return stage1_sample_manifest_path(stage.variant_table)


def has_wgs_manifest(config: ProjectConfig) -> bool:
    return Path(wgs_manifest_path(config)).exists()


def wgs_present_ids(config: ProjectConfig, *, require: bool = True) -> set[str] | None:
    path = wgs_manifest_path(config)
    if not Path(path).exists():
        if require:
            raise RuntimeError(
                f"Missing Stage 1 WGS sample manifest at {path}. Run `aou-workbench prepare-stage1` first."
            )
        return None
    return _load_id_set(path)


def has_max_unrelated_file(config: ProjectConfig) -> bool:
    path = config.workbench.max_unrelated_path
    if not path:
        return False
    if path.startswith("gs://"):
        return True
    return Path(path).exists()


def max_unrelated_ids(config: ProjectConfig, *, require: bool = False) -> set[str] | None:
    path = config.workbench.max_unrelated_path
    if not path:
        return None
    if not path.startswith("gs://") and not Path(path).exists():
        if require:
            raise RuntimeError(f"Configured max-unrelated file does not exist: {path}")
        return None
    return _load_id_set(path)


def gwas_universe_ids(
    config: ProjectConfig,
    *,
    require_wgs: bool = True,
    require_max_unrelated: bool = False,
) -> set[str]:
    wgs_ids = wgs_present_ids(config, require=require_wgs)
    ids = set(wgs_ids or ())
    unrelated_ids = max_unrelated_ids(config, require=require_max_unrelated)
    if unrelated_ids is not None:
        ids = ids.intersection(unrelated_ids) if ids else set(unrelated_ids)
    return ids


def restrict_frame_to_ids(
    frame: pd.DataFrame,
    ids: set[str] | None,
    *,
    id_column: str = "person_id",
) -> pd.DataFrame:
    if ids is None:
        return frame.copy()
    subset = frame.copy()
    subset[id_column] = subset[id_column].astype(str)
    return subset[subset[id_column].isin(ids)].copy()


def restrict_frame_for_gwas(
    config: ProjectConfig,
    frame: pd.DataFrame,
    *,
    require_wgs: bool = True,
    require_max_unrelated: bool = False,
    id_column: str = "person_id",
) -> pd.DataFrame:
    ids = gwas_universe_ids(
        config,
        require_wgs=require_wgs,
        require_max_unrelated=require_max_unrelated,
    )
    return restrict_frame_to_ids(frame, ids, id_column=id_column)


__all__ = [
    "gwas_universe_ids",
    "has_max_unrelated_file",
    "has_wgs_manifest",
    "max_unrelated_ids",
    "restrict_frame_for_gwas",
    "restrict_frame_to_ids",
    "wgs_manifest_path",
    "wgs_present_ids",
]
