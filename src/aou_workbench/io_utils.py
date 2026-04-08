"""Shared I/O helpers for local and Workbench-friendly paths."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
import json
from pathlib import Path
import re
from typing import Any, Mapping

import pandas as pd
import yaml


def slugify(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", value.strip().lower())
    slug = re.sub(r"-{2,}", "-", slug).strip("-")
    return slug or "analysis"


def stable_hash(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, default=str).encode("utf-8").hex()[:16]


def _normalize_payload(payload: Any) -> Any:
    if is_dataclass(payload):
        return asdict(payload)
    return payload


def load_yaml(path: str) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a mapping in YAML file: {path}")
    return payload


def ensure_parent_dir(path: str) -> None:
    if path.startswith("gs://"):
        return
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_table(path: str) -> pd.DataFrame:
    suffix = path.lower()
    if suffix.endswith(".parquet"):
        return pd.read_parquet(path)
    if suffix.endswith(".tsv") or suffix.endswith(".tsv.gz") or suffix.endswith(".bgz"):
        return pd.read_csv(path, sep="\t")
    if suffix.endswith(".csv") or suffix.endswith(".csv.gz"):
        return pd.read_csv(path)
    raise ValueError(f"Unsupported table format: {path}")


def write_dataframe(df: pd.DataFrame, path: str) -> None:
    ensure_parent_dir(path)
    suffix = path.lower()
    if suffix.endswith(".parquet"):
        df.to_parquet(path, index=False)
        return
    if suffix.endswith(".tsv") or suffix.endswith(".tsv.gz") or suffix.endswith(".bgz"):
        df.to_csv(path, sep="\t", index=False)
        return
    if suffix.endswith(".csv") or suffix.endswith(".csv.gz"):
        df.to_csv(path, index=False)
        return
    raise ValueError(f"Unsupported output table format: {path}")


def write_json(payload: Mapping[str, Any], path: str) -> None:
    ensure_parent_dir(path)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(_normalize_payload(payload), handle, indent=2, default=str)


def write_text(text: str, path: str) -> None:
    ensure_parent_dir(path)
    Path(path).write_text(text, encoding="utf-8")


def parse_date(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series, errors="coerce")


__all__ = [
    "ensure_parent_dir",
    "load_yaml",
    "parse_date",
    "read_table",
    "slugify",
    "stable_hash",
    "write_dataframe",
    "write_json",
    "write_text",
]
