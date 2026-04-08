"""Runtime discovery and preflight checks for AoU-style environments."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import json
import os
import subprocess
from typing import Any, Iterable

from .config import DEFAULT_WORKSPACE_CDR, ProjectConfig, WorkbenchConfig


@dataclass(frozen=True)
class RuntimeDefaults:
    workspace_bucket: str | None = None
    workspace_cdr: str | None = None
    genomics_bucket: str | None = None
    requester_pays_project: str | None = None
    requester_pays_buckets: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class PreflightCheck:
    name: str
    status: str
    message: str
    detail: str | None = None

    @property
    def ok(self) -> bool:
        return self.status in {"PASS", "WARN"}


class PreflightError(RuntimeError):
    """Raised when preflight detects an unsafe or incomplete setup."""


def _run_command(cmd: list[str]) -> tuple[int, str, str]:
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except FileNotFoundError as exc:
        return 127, "", str(exc)
    return result.returncode, result.stdout.strip(), result.stderr.strip()


def _load_workbench_resources() -> list[dict[str, Any]]:
    code, stdout, _ = _run_command(["wb", "resource", "list", "--format=JSON"])
    if code != 0 or not stdout:
        return []
    try:
        payload = json.loads(stdout)
    except json.JSONDecodeError:
        return []
    return payload if isinstance(payload, list) else []


def _choose_bucket(resources: Iterable[dict[str, Any]], *, genomics: bool) -> str | None:
    buckets = [item for item in resources if item.get("resourceType") == "GCS_BUCKET"]
    if genomics:
        for item in buckets:
            bucket = (item.get("bucketName") or "").lower()
            rid = (item.get("id") or "").lower()
            if "aou-datasets-controlled" in bucket or "aou-datasets-controlled" in rid:
                return f"gs://{item['bucketName']}"
        return None
    for item in buckets:
        bucket = (item.get("bucketName") or "").lower()
        if "workspace" in bucket or "dataproc-staging" in bucket:
            return f"gs://{item['bucketName']}"
    return f"gs://{buckets[0]['bucketName']}" if buckets else None


def _choose_cdr(resources: Iterable[dict[str, Any]]) -> str | None:
    datasets = [item for item in resources if item.get("resourceType") == "BQ_DATASET"]
    candidates: list[tuple[str, str]] = []
    for item in datasets:
        project_id = item.get("projectId")
        dataset_id = item.get("datasetId")
        if project_id and dataset_id:
            candidates.append((dataset_id, f"{project_id}.{dataset_id}"))
    if not candidates:
        return None
    candidates.sort(key=lambda row: row[0], reverse=True)
    return candidates[0][1]


def discover_runtime_defaults() -> RuntimeDefaults:
    resources = _load_workbench_resources()
    warnings: list[str] = []
    workspace_bucket = _choose_bucket(resources, genomics=False)
    workspace_cdr = _choose_cdr(resources)
    genomics_bucket = _choose_bucket(resources, genomics=True)
    requester_pays_project = (
        os.getenv("GCS_REQUESTER_PAYS_PROJECT")
        or os.getenv("GOOGLE_PROJECT")
        or os.getenv("GOOGLE_CLOUD_PROJECT")
    )
    requester_pays_buckets = ()
    if genomics_bucket and genomics_bucket.startswith("gs://"):
        requester_pays_buckets = (genomics_bucket[5:],)
    if not workspace_cdr:
        warnings.append(
            "Could not auto-discover an attached CDR. Falling back to explicit config or defaults."
        )
    if not workspace_bucket:
        warnings.append(
            "Could not auto-discover a writable workspace bucket. Set workspace_bucket in configs/workbench.yaml."
        )
    return RuntimeDefaults(
        workspace_bucket=workspace_bucket,
        workspace_cdr=workspace_cdr or DEFAULT_WORKSPACE_CDR,
        genomics_bucket=genomics_bucket,
        requester_pays_project=requester_pays_project,
        requester_pays_buckets=requester_pays_buckets,
        warnings=tuple(warnings),
    )


def apply_runtime_defaults(config: ProjectConfig, runtime: RuntimeDefaults | None = None) -> ProjectConfig:
    runtime = runtime or discover_runtime_defaults()
    workbench = config.workbench
    updated = WorkbenchConfig(
        workspace_bucket=workbench.workspace_bucket or runtime.workspace_bucket,
        workspace_cdr=workbench.workspace_cdr or runtime.workspace_cdr,
        requester_pays_project=workbench.requester_pays_project or runtime.requester_pays_project,
        requester_pays_buckets=workbench.requester_pays_buckets or runtime.requester_pays_buckets,
        genomics_bucket=workbench.genomics_bucket or runtime.genomics_bucket or workbench.genomics_bucket,
        storage_root=workbench.storage_root,
        wgs_vds_path=workbench.wgs_vds_path,
        vat_path=workbench.vat_path,
        exome_mt_path=workbench.exome_mt_path,
        clinvar_mt_path=workbench.clinvar_mt_path,
        acaf_mt_path=workbench.acaf_mt_path,
        flagged_samples_path=workbench.flagged_samples_path,
        max_unrelated_path=workbench.max_unrelated_path,
    )
    return replace(config, workbench=updated)


def _check_local_or_gcs_path(path: str, name: str) -> PreflightCheck:
    if path.startswith("gs://"):
        code, stdout, stderr = _run_command(["gsutil", "ls", path])
        if code == 0:
            detail = stdout or path
            return PreflightCheck(name=name, status="PASS", message=f"Can see `{path}`.", detail=detail)
        if code == 127:
            return PreflightCheck(
                name=name,
                status="WARN",
                message=f"Skipped remote path check for `{path}` because gsutil is unavailable.",
                detail=stderr,
            )
        return PreflightCheck(
            name=name,
            status="FAIL",
            message=f"Cannot see `{path}`.",
            detail=stderr or stdout,
        )
    if os.path.exists(path):
        return PreflightCheck(name=name, status="PASS", message=f"Found local path `{path}`.")
    return PreflightCheck(name=name, status="FAIL", message=f"Missing local path `{path}`.")


def _bigquery_check(cdr: str | None) -> PreflightCheck:
    if not cdr:
        return PreflightCheck(
            name="bigquery:person",
            status="WARN",
            message="Skipped BigQuery check because no workspace CDR is configured.",
        )
    try:
        from google.cloud import bigquery  # type: ignore
    except ImportError as exc:
        return PreflightCheck(
            name="bigquery:person",
            status="WARN",
            message="Skipped BigQuery check because google-cloud-bigquery is unavailable.",
            detail=str(exc),
        )
    client = bigquery.Client()
    sql = f"SELECT COUNT(*) AS n FROM `{cdr}.person` LIMIT 1"
    try:
        result = client.query(sql).result().to_dataframe()
        return PreflightCheck(
            name="bigquery:person",
            status="PASS",
            message=f"Can query `{cdr}.person`.",
            detail=f"COUNT(*) probe returned {int(result.iloc[0, 0])}.",
        )
    except Exception as exc:  # pragma: no cover - environment dependent
        return PreflightCheck(
            name="bigquery:person",
            status="FAIL",
            message=f"Cannot query `{cdr}.person`.",
            detail=f"{type(exc).__name__}: {exc}",
        )


def _hail_check(config: ProjectConfig) -> PreflightCheck:
    try:
        import hail as hl  # type: ignore
    except ImportError as exc:
        return PreflightCheck(
            name="hail:init",
            status="WARN",
            message="Skipped Hail initialization because hail is unavailable.",
            detail=str(exc),
        )
    try:  # pragma: no cover - environment dependent
        kwargs: dict[str, Any] = {}
        if config.workbench.requester_pays_project:
            kwargs["gcs_requester_pays_configuration"] = (
                config.workbench.requester_pays_project,
                list(config.workbench.requester_pays_buckets),
            )
        hl.init(**kwargs)
        hl.default_reference(new_default_reference="GRCh38")
        return PreflightCheck(
            name="hail:init",
            status="PASS",
            message="Hail initialized successfully.",
        )
    except Exception as exc:
        return PreflightCheck(
            name="hail:init",
            status="FAIL",
            message="Hail could not initialize cleanly.",
            detail=f"{type(exc).__name__}: {exc}",
        )


def run_preflight_checks(config: ProjectConfig) -> list[PreflightCheck]:
    runtime = discover_runtime_defaults()
    effective = apply_runtime_defaults(config, runtime)
    checks = [
        PreflightCheck(
            name="runtime:workspace_bucket",
            status="PASS" if effective.workbench.workspace_bucket else "WARN",
            message=(
                f"Using workspace bucket `{effective.workbench.workspace_bucket}`."
                if effective.workbench.workspace_bucket
                else "No workspace bucket configured."
            ),
        ),
        PreflightCheck(
            name="runtime:workspace_cdr",
            status="PASS" if effective.workbench.workspace_cdr else "WARN",
            message=(
                f"Using workspace CDR `{effective.workbench.workspace_cdr}`."
                if effective.workbench.workspace_cdr
                else "No workspace CDR configured."
            ),
        ),
        _check_local_or_gcs_path(effective.phenotype.tables.cohort_table, "input:cohort"),
        _check_local_or_gcs_path(effective.phenotype.tables.condition_table, "input:condition"),
        _check_local_or_gcs_path(effective.phenotype.tables.measurement_table, "input:measurement"),
    ]
    if effective.phenotype.tables.ancestry_table:
        checks.append(_check_local_or_gcs_path(effective.phenotype.tables.ancestry_table, "input:ancestry"))
    if effective.analysis.stage1:
        checks.append(_check_local_or_gcs_path(effective.analysis.stage1.variant_table, "input:stage1"))
    if effective.analysis.stage2:
        checks.append(_check_local_or_gcs_path(effective.analysis.stage2.variant_table, "input:stage2"))
    if effective.analysis.stage3:
        checks.append(_check_local_or_gcs_path(effective.analysis.stage3.variant_table, "input:stage3"))
    if effective.analysis.stage4:
        checks.append(_check_local_or_gcs_path(effective.analysis.stage4.genotype_table, "input:stage4"))
    checks.append(_bigquery_check(effective.workbench.workspace_cdr))
    checks.append(_hail_check(effective))
    for warning in runtime.warnings:
        checks.append(
            PreflightCheck(
                name="runtime:warning",
                status="WARN",
                message=warning,
            )
        )
    return checks


def format_preflight_report(checks: Iterable[PreflightCheck]) -> str:
    lines = ["AoU Workbench preflight"]
    for check in checks:
        lines.append(f"[{check.status}] {check.name}: {check.message}")
        if check.detail:
            lines.append(f"  {check.detail}")
    return "\n".join(lines)


def assert_preflight_ok(checks: Iterable[PreflightCheck]) -> None:
    failures = [check for check in checks if not check.ok]
    if failures:
        raise PreflightError(format_preflight_report(failures))


__all__ = [
    "PreflightCheck",
    "PreflightError",
    "RuntimeDefaults",
    "apply_runtime_defaults",
    "assert_preflight_ok",
    "discover_runtime_defaults",
    "format_preflight_report",
    "run_preflight_checks",
]
