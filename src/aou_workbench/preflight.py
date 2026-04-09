"""Runtime discovery and preflight checks for AoU-style environments."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import json
import os
import shutil
import subprocess
from typing import Any, Iterable

from .config import DEFAULT_WORKSPACE_CDR, ProjectConfig, WorkbenchConfig
from .io_utils import is_bigquery_table


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
    candidates: list[tuple[tuple[int, str], str]] = []
    for item in datasets:
        project_id = item.get("projectId")
        dataset_id = item.get("datasetId")
        if project_id and dataset_id:
            lowered = dataset_id.lower()
            rank = 2
            if lowered.startswith("prep_"):
                rank = 1
            elif lowered.startswith("c20"):
                rank = 3
            candidates.append(((rank, dataset_id), f"{project_id}.{dataset_id}"))
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
    if not requester_pays_project:
        for item in resources:
            project_id = item.get("projectId")
            if item.get("resourceType") in {"DATAPROC_CLUSTER", "GCS_BUCKET"} and project_id:
                requester_pays_project = str(project_id)
                break
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
    if updated.requester_pays_project:
        os.environ.setdefault("GCS_REQUESTER_PAYS_PROJECT", updated.requester_pays_project)
    return replace(config, workbench=updated)


def _check_local_or_gcs_path(path: str, name: str, requester_pays_project: str | None = None) -> PreflightCheck:
    if path.startswith("gs://"):
        cmd = ["gsutil"]
        project = requester_pays_project or os.getenv("GCS_REQUESTER_PAYS_PROJECT")
        if project:
            cmd.extend(["-u", project])
        cmd.extend(["ls", path])
        code, stdout, stderr = _run_command(cmd)
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


def _qualify_cdr_table(cdr: str | None, reference: str) -> str | None:
    if reference.startswith("gs://"):
        return None
    normalized = reference[5:] if reference.startswith("bq://") else reference
    if is_bigquery_table(normalized):
        return normalized
    if "." in normalized or not cdr:
        return None
    return f"{cdr}.{normalized}"


def _bigquery_table_check(cdr: str | None, table_reference: str, name: str) -> PreflightCheck:
    qualified = _qualify_cdr_table(cdr, table_reference)
    if not qualified:
        return _check_local_or_gcs_path(table_reference, name)
    try:
        from google.cloud import bigquery  # type: ignore
    except ImportError as exc:
        return PreflightCheck(
            name=name,
            status="WARN",
            message=f"Skipped BigQuery table check for `{qualified}` because google-cloud-bigquery is unavailable.",
            detail=str(exc),
        )
    client = bigquery.Client()
    sql = f"SELECT 1 AS ok FROM `{qualified}` LIMIT 1"
    try:
        list(client.query(sql).result())
        return PreflightCheck(name=name, status="PASS", message=f"Can query `{qualified}`.")
    except Exception as exc:  # pragma: no cover - environment dependent
        return PreflightCheck(
            name=name,
            status="FAIL",
            message=f"Cannot query `{qualified}`.",
            detail=f"{type(exc).__name__}: {exc}",
        )


def _check_input_reference(cdr: str | None, reference: str, name: str) -> PreflightCheck:
    if _qualify_cdr_table(cdr, reference) or is_bigquery_table(reference):
        return _bigquery_table_check(cdr, reference, name)
    return _check_local_or_gcs_path(reference, name, os.getenv("GCS_REQUESTER_PAYS_PROJECT"))


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


def _hail_check() -> PreflightCheck:
    try:
        import hail  # type: ignore
    except ImportError as exc:
        return PreflightCheck(
            name="hail:available",
            status="WARN",
            message="Hail is unavailable in this environment.",
            detail=str(exc),
        )
    version = getattr(hail, "__version__", "unknown")
    return PreflightCheck(
        name="hail:available",
        status="PASS",
        message="Hail is importable for direct WGS VDS extraction and broader MT/VDS analyses.",
        detail=f"Version: {version}",
    )


def _planned_output_check(path: str, name: str) -> PreflightCheck:
    if path.startswith("gs://"):
        return PreflightCheck(
            name=name,
            status="PASS",
            message=f"Will write output to `{path}`.",
        )
    parent = os.path.dirname(path) or "."
    if os.path.exists(parent):
        return PreflightCheck(
            name=name,
            status="PASS",
            message=f"Will write local output to `{path}`.",
            detail=f"Parent directory exists: {parent}",
        )
    return PreflightCheck(
        name=name,
        status="PASS",
        message=f"Will write local output to `{path}`.",
        detail=f"Parent directory will be created on write: {parent}",
    )


def _tool_check(command: str, *, name: str, required_for: str, override_env: str | None = None) -> PreflightCheck:
    path = os.getenv(override_env) if override_env else None
    if not path:
        path = shutil.which(command)
    if path:
        return PreflightCheck(
            name=name,
            status="PASS",
            message=f"`{command}` is available.",
            detail=path,
        )
    return PreflightCheck(
        name=name,
        status="WARN",
        message=f"`{command}` is unavailable.",
        detail=f"Required for {required_for}. Use an AoU app image with this tool or install it before running that stage.",
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
    ]
    if effective.phenotype.tables.cohort_table:
        checks.append(_check_local_or_gcs_path(effective.phenotype.tables.cohort_table, "input:cohort"))
    else:
        checks.extend(
            [
                _bigquery_table_check(effective.workbench.workspace_cdr, effective.phenotype.tables.person_table, "input:person"),
                _bigquery_table_check(
                    effective.workbench.workspace_cdr,
                    effective.phenotype.tables.observation_table,
                    "input:observation_period",
                ),
                _bigquery_table_check(
                    effective.workbench.workspace_cdr,
                    effective.phenotype.tables.condition_table,
                    "input:condition_occurrence",
                ),
                _bigquery_table_check(
                    effective.workbench.workspace_cdr,
                    effective.phenotype.tables.measurement_table,
                    "input:measurement",
                ),
            ]
        )
    if effective.phenotype.tables.ancestry_table:
        checks.append(
            _check_input_reference(
                effective.workbench.workspace_cdr,
                effective.phenotype.tables.ancestry_table,
                "input:ancestry",
            )
        )
    if effective.phenotype.clinical_cofactors:
        checks.append(
            _bigquery_table_check(
                effective.workbench.workspace_cdr,
                effective.phenotype.tables.concept_table,
                "input:concept",
            )
        )
    if effective.analysis.stage1:
        checks.append(
            _check_local_or_gcs_path(
                effective.workbench.wgs_vds_path,
                "input:wgs_vds",
                effective.workbench.requester_pays_project,
            )
        )
        checks.append(_planned_output_check(effective.analysis.stage1.variant_table, "output:stage1"))
    if effective.analysis.stage2:
        checks.append(
            _check_local_or_gcs_path(
                effective.workbench.clinvar_mt_path,
                "input:clinvar_mt",
                effective.workbench.requester_pays_project,
            )
        )
        checks.append(_planned_output_check(effective.analysis.stage2.variant_table, "output:stage2"))
    if effective.analysis.run_stage3 and effective.analysis.stage3:
        checks.append(_check_input_reference(effective.workbench.workspace_cdr, effective.analysis.stage3.variant_table, "input:stage3"))
    if effective.analysis.stage4:
        checks.append(
            _check_local_or_gcs_path(
                effective.workbench.acaf_mt_path,
                "input:acaf_mt",
                effective.workbench.requester_pays_project,
            )
        )
        checks.append(_planned_output_check(effective.analysis.stage4.genotype_table, "output:stage4"))
    checks.append(_bigquery_check(effective.workbench.workspace_cdr))
    checks.append(_tool_check("gsutil", name="tool:gsutil", required_for="Workbench bucket and genomics bucket access"))
    checks.append(_hail_check())
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
