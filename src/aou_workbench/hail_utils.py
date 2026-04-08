"""Helpers that isolate optional Hail usage for genomics extraction."""

from __future__ import annotations

import os
import shlex
import subprocess
from typing import Any, Iterable


HAIL_SPARK_SERIALIZER = "org.apache.spark.serializer.KryoSerializer"
HAIL_KRYO_REGISTRATOR = "is.hail.kryo.HailKryoRegistrator"


def require_hail() -> Any:
    try:
        import hail as hl  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "Hail is required for Stage 1 extraction. Use a Hail-capable AoU environment or "
            "install Hail in a separate extraction environment."
        ) from exc
    return hl


def requester_pays_configuration(
    requester_pays_project: str | None,
    requester_pays_buckets: Iterable[str] = (),
) -> str | tuple[str, list[str]] | None:
    if not requester_pays_project:
        return None
    buckets = sorted(
        {
            bucket[5:] if bucket.startswith("gs://") else bucket
            for bucket in requester_pays_buckets
            if bucket
        }
    )
    if not buckets:
        return requester_pays_project
    return requester_pays_project, buckets


def _ensure_pyspark_submit_args() -> None:
    required = [
        "--conf",
        f"spark.serializer={HAIL_SPARK_SERIALIZER}",
        "--conf",
        f"spark.kryo.registrator={HAIL_KRYO_REGISTRATOR}",
    ]
    existing = os.environ.get("PYSPARK_SUBMIT_ARGS", "").strip()
    existing_parts = shlex.split(existing) if existing else []

    has_serializer = any(f"spark.serializer={HAIL_SPARK_SERIALIZER}" in part for part in existing_parts)
    has_registrator = any(f"spark.kryo.registrator={HAIL_KRYO_REGISTRATOR}" in part for part in existing_parts)

    parts: list[str] = []
    if not has_serializer:
        parts.extend(required[:2])
    if not has_registrator:
        parts.extend(required[2:])
    parts.extend(existing_parts or ["pyspark-shell"])
    if "pyspark-shell" not in parts:
        parts.append("pyspark-shell")
    os.environ["PYSPARK_SUBMIT_ARGS"] = " ".join(parts)


def _run_hailctl_config(key: str, value: str) -> None:
    try:
        subprocess.run(
            ["hailctl", "config", "set", key, value],
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        return


def configure_aou_hail_bootstrap(
    requester_pays_project: str | None,
    requester_pays_buckets: Iterable[str] = (),
) -> None:
    _ensure_pyspark_submit_args()
    if not requester_pays_project:
        return
    _run_hailctl_config("gcs_requester_pays/project", requester_pays_project)
    buckets = sorted(
        {
            bucket[5:] if bucket.startswith("gs://") else bucket
            for bucket in requester_pays_buckets
            if bucket
        }
    )
    if buckets:
        _run_hailctl_config("gcs_requester_pays/buckets", ",".join(buckets))


def ensure_hail(
    reference: str = "GRCh38",
    *,
    requester_pays_project: str | None = None,
    requester_pays_buckets: Iterable[str] = (),
) -> Any:
    hl = require_hail()
    configure_aou_hail_bootstrap(requester_pays_project, requester_pays_buckets)

    init_kwargs: dict[str, Any] = {}
    rp_config = requester_pays_configuration(requester_pays_project, requester_pays_buckets)
    if rp_config is not None:
        init_kwargs["gcs_requester_pays_configuration"] = rp_config

    try:
        hl.current_backend()
    except Exception:
        hl.init(**init_kwargs)
    hl.default_reference(reference)
    return hl


__all__ = [
    "configure_aou_hail_bootstrap",
    "ensure_hail",
    "requester_pays_configuration",
    "require_hail",
]
