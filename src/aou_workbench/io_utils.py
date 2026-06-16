"""Shared I/O helpers for local and Workbench-friendly paths."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
import csv
import io
import json
import os
from pathlib import Path
import re
import subprocess
import time
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


def is_bigquery_table(path: str) -> bool:
    normalized = path[5:] if path.startswith("bq://") else path
    return bool(re.fullmatch(r"[^./`]+[.][^./`]+[.][^./`]+", normalized))


def query_bigquery_dataframe(
    sql: str,
    *,
    maximum_bytes_billed: int | None = None,
    progress_label: str | None = None,
) -> pd.DataFrame:
    try:
        from google.cloud import bigquery  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError("google-cloud-bigquery is required for BigQuery-backed inputs.") from exc
    client = bigquery.Client()
    job_config = None
    if maximum_bytes_billed is not None:
        job_config = bigquery.QueryJobConfig(maximum_bytes_billed=maximum_bytes_billed)
    job = client.query(sql, job_config=job_config)
    if progress_label:
        print(f"{progress_label}: submitted BigQuery job {job.job_id} in {job.location}", flush=True)
    started = time.monotonic()
    while not job.done(reload=True):
        if progress_label:
            elapsed = int(time.monotonic() - started)
            print(f"{progress_label}: BigQuery job {job.job_id} still running after {elapsed}s", flush=True)
        time.sleep(30)
    result = job.result()
    if progress_label:
        elapsed = int(time.monotonic() - started)
        bytes_processed = int(job.total_bytes_processed or 0)
        bytes_billed = int(job.total_bytes_billed or 0)
        print(
            f"{progress_label}: BigQuery job {job.job_id} finished after {elapsed}s "
            f"({bytes_processed} bytes processed, {bytes_billed} bytes billed)",
            flush=True,
        )
        print(f"{progress_label}: downloading result rows", flush=True)
    rows = [dict(row.items()) for row in result]
    if progress_label:
        print(f"{progress_label}: downloaded {len(rows)} rows", flush=True)
    return pd.DataFrame(rows)


def query_bigquery_to_tsv(
    sql: str,
    path: str,
    *,
    maximum_bytes_billed: int | None = None,
    progress_label: str | None = None,
    progress_every: int = 100_000,
) -> dict[str, Any]:
    try:
        from google.cloud import bigquery  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError("google-cloud-bigquery is required for BigQuery-backed inputs.") from exc
    client = bigquery.Client()
    job_config = None
    if maximum_bytes_billed is not None:
        job_config = bigquery.QueryJobConfig(maximum_bytes_billed=maximum_bytes_billed)
    job = client.query(sql, job_config=job_config)
    if progress_label:
        print(f"{progress_label}: submitted BigQuery job {job.job_id} in {job.location}", flush=True)
    started = time.monotonic()
    while not job.done(reload=True):
        if progress_label:
            elapsed = int(time.monotonic() - started)
            print(f"{progress_label}: BigQuery job {job.job_id} still running after {elapsed}s", flush=True)
        time.sleep(30)
    result = job.result()
    elapsed = int(time.monotonic() - started)
    bytes_processed = int(job.total_bytes_processed or 0)
    bytes_billed = int(job.total_bytes_billed or 0)
    if progress_label:
        print(
            f"{progress_label}: BigQuery job {job.job_id} finished after {elapsed}s "
            f"({bytes_processed} bytes processed, {bytes_billed} bytes billed)",
            flush=True,
        )
        print(f"{progress_label}: streaming result rows to {path}", flush=True)

    ensure_parent_dir(path)
    fieldnames = [field.name for field in result.schema]
    row_count = 0
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in result:
            writer.writerow({key: ("" if value is None else value) for key, value in dict(row.items()).items()})
            row_count += 1
            if progress_label and progress_every > 0 and row_count % progress_every == 0:
                print(f"{progress_label}: streamed {row_count} rows", flush=True)
    if progress_label:
        print(f"{progress_label}: streamed {row_count} rows to {path}", flush=True)
    return {
        "job_id": job.job_id,
        "location": job.location,
        "total_bytes_processed": bytes_processed,
        "total_bytes_billed": bytes_billed,
        "row_count": row_count,
        "path": path,
        "elapsed_seconds": elapsed,
    }


def dry_run_bigquery_query(sql: str, *, maximum_bytes_billed: int | None = None) -> dict[str, Any]:
    try:
        from google.cloud import bigquery  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError("google-cloud-bigquery is required for BigQuery-backed inputs.") from exc
    client = bigquery.Client()
    job_config = bigquery.QueryJobConfig(
        dry_run=True,
        use_query_cache=False,
    )
    job = client.query(sql, job_config=job_config)
    total_bytes = int(job.total_bytes_processed or 0)
    return {
        "total_bytes_processed": total_bytes,
        "total_tib_processed": total_bytes / float(1024**4),
        "maximum_bytes_billed": maximum_bytes_billed,
        "would_exceed_maximum_bytes_billed": (
            maximum_bytes_billed is not None and total_bytes > maximum_bytes_billed
        ),
    }


def _gsutil_cat(path: str) -> str:
    requester_pays_project = (
        os.getenv("GCS_REQUESTER_PAYS_PROJECT")
        or os.getenv("GOOGLE_PROJECT")
        or os.getenv("GOOGLE_CLOUD_PROJECT")
    )
    cmd = ["gsutil"]
    if requester_pays_project:
        cmd.extend(["-u", requester_pays_project])
    cmd.extend(["cat", path])
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip() or f"Failed to read {path}")
    return result.stdout


def _resolve_gcs_path(path: str) -> str:
    if not any(token in path for token in ("*", "?", "[")):
        return path
    requester_pays_project = (
        os.getenv("GCS_REQUESTER_PAYS_PROJECT")
        or os.getenv("GOOGLE_PROJECT")
        or os.getenv("GOOGLE_CLOUD_PROJECT")
    )
    cmd = ["gsutil"]
    if requester_pays_project:
        cmd.extend(["-u", requester_pays_project])
    cmd.extend(["ls", path])
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip() or f"Failed to list {path}")
    matches = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if not matches:
        raise RuntimeError(f"No URLs matched: {path}")
    return sorted(matches)[0]


def read_table(path: str) -> pd.DataFrame:
    suffix = path.lower()
    if is_bigquery_table(path):
        reference = path[5:] if path.startswith("bq://") else path
        return query_bigquery_dataframe(f"SELECT * FROM `{reference}`")
    if path.startswith("gs://"):
        resolved = _resolve_gcs_path(path)
        suffix = resolved.lower()
        if suffix.endswith(".tsv") or suffix.endswith(".tsv.gz") or suffix.endswith(".bgz") or suffix.endswith(".txt"):
            return pd.read_csv(io.StringIO(_gsutil_cat(resolved)), sep="\t", low_memory=False)
        if suffix.endswith(".csv") or suffix.endswith(".csv.gz"):
            return pd.read_csv(io.StringIO(_gsutil_cat(resolved)), low_memory=False)
        raise ValueError(f"Unsupported remote table format: {resolved}")
    if suffix.endswith(".parquet"):
        return pd.read_parquet(path)
    if suffix.endswith(".tsv") or suffix.endswith(".tsv.gz") or suffix.endswith(".bgz") or suffix.endswith(".txt"):
        return pd.read_csv(path, sep="\t", low_memory=False)
    if suffix.endswith(".csv") or suffix.endswith(".csv.gz"):
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"Unsupported table format: {path}")


def write_dataframe(df: pd.DataFrame, path: str) -> None:
    ensure_parent_dir(path)
    suffix = path.lower()
    if suffix.endswith(".parquet"):
        df.to_parquet(path, index=False)
        return
    if suffix.endswith(".tsv") or suffix.endswith(".tsv.gz") or suffix.endswith(".bgz") or suffix.endswith(".txt"):
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
    "is_bigquery_table",
    "load_yaml",
    "parse_date",
    "dry_run_bigquery_query",
    "query_bigquery_dataframe",
    "query_bigquery_to_tsv",
    "read_table",
    "slugify",
    "stable_hash",
    "write_dataframe",
    "write_json",
    "write_text",
]
