"""Helpers that isolate optional Hail usage for genomics extraction."""

from __future__ import annotations

import os
import shlex
import subprocess
from typing import Any, Iterable


HAIL_SPARK_SERIALIZER = "org.apache.spark.serializer.KryoSerializer"
HAIL_KRYO_REGISTRATOR = "is.hail.kryo.HailKryoRegistrator"
REQUESTER_PAYS_PROJECT_ENVS = (
    "GCS_REQUESTER_PAYS_PROJECT",
    "GOOGLE_CLOUD_PROJECT",
    "GOOGLE_PROJECT",
    "GCLOUD_PROJECT",
)


def _normalize_requester_pays_buckets(requester_pays_buckets: Iterable[str] = ()) -> list[str]:
    return sorted(
        {
            bucket[5:] if bucket.startswith("gs://") else bucket
            for bucket in requester_pays_buckets
            if bucket
        }
    )


def resolve_requester_pays_project(requester_pays_project: str | None = None) -> str | None:
    if requester_pays_project:
        return requester_pays_project
    for key in REQUESTER_PAYS_PROJECT_ENVS:
        value = os.getenv(key)
        if value:
            return value
    try:
        result = subprocess.run(
            ["gcloud", "config", "get-value", "project"],
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        return None
    project = result.stdout.strip()
    if result.returncode == 0 and project and project != "(unset)":
        return project
    return None


def require_hail() -> Any:
    # Hail 0.2.137 still references deprecated NumPy aliases like np.bool.
    # Workbench images commonly carry NumPy >=1.24 where those aliases were removed.
    try:
        import numpy as np  # type: ignore

        if not hasattr(np, "bool"):
            np.bool = np.bool_  # type: ignore[attr-defined]
    except Exception:
        pass
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
    buckets = _normalize_requester_pays_buckets(requester_pays_buckets)
    if not buckets:
        return requester_pays_project
    return requester_pays_project, buckets


def _spark_requester_pays_conf(
    requester_pays_project: str | None,
    requester_pays_buckets: Iterable[str] = (),
) -> dict[str, str]:
    if not requester_pays_project:
        return {}
    buckets = _normalize_requester_pays_buckets(requester_pays_buckets)
    conf = {
        "spark.hadoop.fs.gs.requester.pays.project.id": requester_pays_project,
        "spark.hadoop.fs.gs.requester.pays.mode": "CUSTOM" if buckets else "ENABLED",
    }
    if buckets:
        conf["spark.hadoop.fs.gs.requester.pays.buckets"] = ",".join(buckets)
    return conf


def _append_missing_spark_conf(parts: list[str], key: str, value: str) -> None:
    assignment = f"{key}={value}"
    if any(part == assignment or part.startswith(f"{key}=") for part in parts):
        return
    parts.extend(["--conf", assignment])


def _ensure_pyspark_submit_args(
    requester_pays_project: str | None = None,
    requester_pays_buckets: Iterable[str] = (),
) -> None:
    existing = os.environ.get("PYSPARK_SUBMIT_ARGS", "").strip()
    existing_parts = shlex.split(existing) if existing else []

    desired = [
        ("spark.serializer", HAIL_SPARK_SERIALIZER),
        ("spark.kryo.registrator", HAIL_KRYO_REGISTRATOR),
    ]
    for key, value in _spark_requester_pays_conf(requester_pays_project, requester_pays_buckets).items():
        desired.append((key, value))

    parts: list[str] = []
    for key, value in desired:
        if any(part.startswith(f"{key}=") for part in existing_parts):
            continue
        _append_missing_spark_conf(parts, key, value)
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
    requester_pays_project = resolve_requester_pays_project(requester_pays_project)
    if requester_pays_project:
        os.environ.setdefault("GCS_REQUESTER_PAYS_PROJECT", requester_pays_project)
    _ensure_pyspark_submit_args(requester_pays_project, requester_pays_buckets)
    if not requester_pays_project:
        return
    _run_hailctl_config("gcs_requester_pays/project", requester_pays_project)
    buckets = _normalize_requester_pays_buckets(requester_pays_buckets)
    if buckets:
        _run_hailctl_config("gcs_requester_pays/buckets", ",".join(buckets))


def _set_active_spark_requester_pays_conf(
    hl: Any,
    requester_pays_project: str | None,
    requester_pays_buckets: Iterable[str] = (),
) -> None:
    if not requester_pays_project:
        return
    spark_conf = _spark_requester_pays_conf(requester_pays_project, requester_pays_buckets)
    hadoop_conf = {
        key.removeprefix("spark.hadoop."): value
        for key, value in spark_conf.items()
        if key.startswith("spark.hadoop.")
    }
    if not hadoop_conf:
        return
    try:
        backend = hl.current_backend()
    except Exception:
        return
    candidates = [
        getattr(backend, "sc", None),
        getattr(backend, "_sc", None),
    ]
    try:
        from pyspark import SparkContext  # type: ignore

        candidates.append(SparkContext._active_spark_context)
    except Exception:
        pass
    for sc in candidates:
        if sc is None:
            continue
        try:
            jconf = sc._jsc.hadoopConfiguration()
        except Exception:
            continue
        for key, value in hadoop_conf.items():
            jconf.set(key, value)
        return


def ensure_hail(
    reference: str = "GRCh38",
    *,
    requester_pays_project: str | None = None,
    requester_pays_buckets: Iterable[str] = (),
) -> Any:
    hl = require_hail()
    requester_pays_project = resolve_requester_pays_project(requester_pays_project)
    configure_aou_hail_bootstrap(requester_pays_project, requester_pays_buckets)

    init_kwargs: dict[str, Any] = {}
    rp_config = requester_pays_configuration(requester_pays_project, requester_pays_buckets)
    if rp_config is not None:
        init_kwargs["gcs_requester_pays_configuration"] = rp_config

    try:
        hl.current_backend()
    except Exception:
        hl.init(**init_kwargs)
    _set_active_spark_requester_pays_conf(hl, requester_pays_project, requester_pays_buckets)
    hl.default_reference(reference)
    return hl


__all__ = [
    "configure_aou_hail_bootstrap",
    "ensure_hail",
    "requester_pays_configuration",
    "resolve_requester_pays_project",
    "require_hail",
]
