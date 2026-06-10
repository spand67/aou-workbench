"""Pre-index clinical data availability profiling for rhabdomyolysis cases."""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
import pandas as pd

from .config import ProjectConfig
from .io_utils import parse_date, query_bigquery_dataframe, read_table, write_dataframe, write_text
from .paths import ProjectPaths, join_path
from .phenotype_sql import _qualify_table
from .reporting import dataframe_markdown


DEFAULT_WINDOWS: tuple[int | None, ...] = (365, 1095, None)
DEFAULT_BIOMARKER_TERMS: dict[str, tuple[str, ...]] = {
    "creatine_kinase": ("creatine kinase",),
    "creatinine": ("creatinine",),
    "egfr": ("estimated glomerular filtration", "egfr"),
    "ast": ("aspartate aminotransferase",),
    "alt": ("alanine aminotransferase",),
    "potassium": ("potassium",),
    "sodium": ("sodium",),
    "calcium": ("calcium",),
    "phosphate": ("phosphate", "phosphorus"),
    "bicarbonate": ("bicarbonate", "carbon dioxide"),
    "magnesium": ("magnesium",),
    "tsh": ("thyrotropin", "thyroid stimulating hormone"),
    "glucose": ("glucose",),
    "hemoglobin_a1c": ("hemoglobin a1c", "glycated hemoglobin"),
    "lactate": ("lactate",),
}


def preindex_case_profile_root(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "clinical", "preindex_case_profile")


def preindex_summary_path(paths: ProjectPaths) -> str:
    return join_path(preindex_case_profile_root(paths), "availability_summary.tsv")


def preindex_condition_top_path(paths: ProjectPaths) -> str:
    return join_path(preindex_case_profile_root(paths), "top_conditions.tsv")


def preindex_measurement_top_path(paths: ProjectPaths) -> str:
    return join_path(preindex_case_profile_root(paths), "top_measurements.tsv")


def preindex_biomarker_path(paths: ProjectPaths) -> str:
    return join_path(preindex_case_profile_root(paths), "biomarker_availability.tsv")


def preindex_report_path(paths: ProjectPaths) -> str:
    return join_path(preindex_case_profile_root(paths), "report.md")


def _window_label(days: int | None) -> str:
    return "any_prior" if days is None else f"{days}d"


def _parse_windows(values: Iterable[str | int | None] | None) -> tuple[int | None, ...]:
    if values is None:
        return DEFAULT_WINDOWS
    windows: list[int | None] = []
    for value in values:
        if value is None:
            windows.append(None)
            continue
        text = str(value).strip().lower()
        if text in {"all", "any", "any_prior", ""}:
            windows.append(None)
        else:
            windows.append(int(text))
    return tuple(dict.fromkeys(windows))


def _case_frame(config: ProjectConfig, cohort_df: pd.DataFrame, case_tier: str | None) -> pd.DataFrame:
    tier = case_tier or config.cohort.primary_case_tier
    cases = cohort_df[cohort_df["case_tier"].astype(str) == tier].copy()
    if cases.empty:
        raise ValueError(f"No cases found with case_tier={tier!r}.")
    cases["person_id"] = cases["person_id"].astype(str)
    cases["index_date"] = parse_date(cases["index_date"])
    cases["obs_start_date"] = parse_date(cases["obs_start_date"])
    cases = cases[cases["index_date"].notna()].copy()
    if cases.empty:
        raise ValueError(f"No {tier!r} cases have an index_date.")
    cases["preindex_observation_days"] = (cases["index_date"] - cases["obs_start_date"]).dt.days
    return cases[["person_id", "case_tier", "index_date", "obs_start_date", "preindex_observation_days"]]


def _summarize_event_domain(
    *,
    events: pd.DataFrame,
    cases: pd.DataFrame,
    domain: str,
    windows: tuple[int | None, ...],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    n_cases = cases["person_id"].nunique()
    for days in windows:
        view = events.copy()
        if days is not None:
            view = view[view["days_before_index"].between(1, days, inclusive="both")]
        people = view.groupby("person_id").size() if not view.empty else pd.Series(dtype=int)
        rows.append(
            {
                "domain": domain,
                "window": _window_label(days),
                "n_cases": n_cases,
                "n_cases_with_any": int(people.shape[0]),
                "pct_cases_with_any": round(100 * people.shape[0] / n_cases, 2) if n_cases else 0.0,
                "total_events": int(view.shape[0]),
                "median_events_per_case_with_any": float(people.median()) if not people.empty else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _top_concepts(
    *,
    events: pd.DataFrame,
    domain: str,
    windows: tuple[int | None, ...],
    top_n: int,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    group_cols = ["concept_id", "concept_name"]
    for days in windows:
        view = events.copy()
        if days is not None:
            view = view[view["days_before_index"].between(1, days, inclusive="both")]
        if view.empty:
            continue
        grouped = (
            view.groupby(group_cols, dropna=False)
            .agg(
                n_cases=("person_id", "nunique"),
                total_events=("person_id", "size"),
                median_days_before_index=("days_before_index", "median"),
            )
            .reset_index()
            .sort_values(["n_cases", "total_events", "concept_name"], ascending=[False, False, True])
            .head(top_n)
        )
        grouped.insert(0, "window", _window_label(days))
        grouped.insert(0, "domain", domain)
        rows.append(grouped)
    if not rows:
        return pd.DataFrame(
            columns=[
                "domain",
                "window",
                "concept_id",
                "concept_name",
                "n_cases",
                "total_events",
                "median_days_before_index",
            ]
        )
    return pd.concat(rows, ignore_index=True)


def _condition_events_local(config: ProjectConfig, cases: pd.DataFrame) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.condition_table).copy()
    raw["person_id"] = raw[config.phenotype.person_id_column].astype(str)
    raw["event_date"] = parse_date(raw[config.phenotype.condition_date_column])
    raw["concept_id"] = raw[config.phenotype.condition_concept_column].astype(str)
    raw["concept_name"] = raw.get("condition_concept_name", raw["concept_id"]).astype(str)
    merged = raw.merge(cases[["person_id", "index_date"]], on="person_id", how="inner")
    merged["days_before_index"] = (merged["index_date"] - merged["event_date"]).dt.days
    return merged[merged["days_before_index"] > 0][
        ["person_id", "event_date", "days_before_index", "concept_id", "concept_name"]
    ].copy()


def _measurement_events_local(config: ProjectConfig, cases: pd.DataFrame) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.measurement_table).copy()
    raw["person_id"] = raw[config.phenotype.person_id_column].astype(str)
    raw["event_date"] = parse_date(raw[config.phenotype.measurement_date_column])
    raw["concept_id"] = raw[config.phenotype.measurement_concept_column].astype(str)
    raw["concept_name"] = raw.get("measurement_concept_name", raw["concept_id"]).astype(str)
    raw["value_as_number"] = pd.to_numeric(raw.get(config.phenotype.measurement_value_column), errors="coerce")
    merged = raw.merge(cases[["person_id", "index_date"]], on="person_id", how="inner")
    merged["days_before_index"] = (merged["index_date"] - merged["event_date"]).dt.days
    return merged[merged["days_before_index"] > 0][
        ["person_id", "event_date", "days_before_index", "concept_id", "concept_name", "value_as_number"]
    ].copy()


def _biomarker_availability_local(
    measurements: pd.DataFrame,
    cases: pd.DataFrame,
    windows: tuple[int | None, ...],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    n_cases = cases["person_id"].nunique()
    lowered = measurements["concept_name"].astype(str).str.lower() if not measurements.empty else pd.Series(dtype=str)
    for biomarker, terms in DEFAULT_BIOMARKER_TERMS.items():
        if measurements.empty:
            subset = measurements
        else:
            mask = pd.Series(False, index=measurements.index)
            for term in terms:
                mask |= lowered.str.contains(term.lower(), regex=False, na=False)
            subset = measurements[mask].copy()
        for days in windows:
            view = subset
            if days is not None:
                view = view[view["days_before_index"].between(1, days, inclusive="both")]
            people = view.groupby("person_id").size() if not view.empty else pd.Series(dtype=int)
            numeric = view["value_as_number"].dropna() if "value_as_number" in view else pd.Series(dtype=float)
            rows.append(
                {
                    "biomarker": biomarker,
                    "window": _window_label(days),
                    "n_cases": n_cases,
                    "n_cases_with_measurement": int(people.shape[0]),
                    "pct_cases_with_measurement": round(100 * people.shape[0] / n_cases, 2) if n_cases else 0.0,
                    "total_measurements": int(view.shape[0]),
                    "median_latest_value": float(numeric.median()) if not numeric.empty else np.nan,
                }
            )
    return pd.DataFrame(rows)


def _literal_case_structs(cases: pd.DataFrame) -> str:
    structs: list[str] = []
    for row in cases.itertuples(index=False):
        index_date = pd.Timestamp(row.index_date).date().isoformat()
        obs_start = pd.Timestamp(row.obs_start_date).date().isoformat() if pd.notna(row.obs_start_date) else index_date
        structs.append(
            "STRUCT("
            f"'{str(row.person_id).replace(chr(39), '')}' AS person_id, "
            f"DATE '{index_date}' AS index_date, "
            f"DATE '{obs_start}' AS obs_start_date)"
        )
    return ",\n    ".join(structs)


def _case_cte_sql(cases: pd.DataFrame) -> str:
    return f"""
cases AS (
  SELECT *
  FROM UNNEST([
    {_literal_case_structs(cases)}
  ])
)
""".strip()


def _window_defs_sql(windows: tuple[int | None, ...]) -> str:
    values = []
    for days in windows:
        if days is None:
            values.append("STRUCT('any_prior' AS window, CAST(NULL AS INT64) AS days)")
        else:
            values.append(f"STRUCT('{days}d' AS window, {days} AS days)")
    return "window_defs AS (\n  SELECT * FROM UNNEST([\n    " + ",\n    ".join(values) + "\n  ])\n)"


def _render_top_concepts_sql(
    config: ProjectConfig,
    cases: pd.DataFrame,
    *,
    domain: str,
    windows: tuple[int | None, ...],
    top_n: int,
) -> str:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    concept_table = _qualify_table(cdr, config.phenotype.tables.concept_table)
    if domain == "condition":
        event_table = _qualify_table(cdr, config.phenotype.tables.condition_table)
        concept_column = config.phenotype.condition_concept_column
        date_column = config.phenotype.condition_date_column
    elif domain == "measurement":
        event_table = _qualify_table(cdr, config.phenotype.tables.measurement_table)
        concept_column = config.phenotype.measurement_concept_column
        date_column = config.phenotype.measurement_date_column
    else:
        raise ValueError(f"Unsupported domain: {domain}")
    value_column = (
        f"SAFE_CAST(e.{config.phenotype.measurement_value_column} AS FLOAT64) AS value_as_number,"
        if domain == "measurement"
        else "CAST(NULL AS FLOAT64) AS value_as_number,"
    )
    return f"""
WITH
{_case_cte_sql(cases)},
{_window_defs_sql(windows)},
events AS (
  SELECT
    c.person_id,
    CAST(e.{concept_column} AS STRING) AS concept_id,
    concept.concept_name AS concept_name,
    DATE_DIFF(c.index_date, DATE(e.{date_column}), DAY) AS days_before_index,
    {value_column}
  FROM cases c
  JOIN `{event_table}` e
    ON CAST(e.{config.phenotype.person_id_column} AS STRING) = c.person_id
  LEFT JOIN `{concept_table}` concept
    ON e.{concept_column} = concept.concept_id
  WHERE DATE(e.{date_column}) < c.index_date
)
SELECT
  '{domain}' AS domain,
  window_defs.window,
  events.concept_id,
  COALESCE(events.concept_name, events.concept_id) AS concept_name,
  COUNT(DISTINCT events.person_id) AS n_cases,
  COUNT(*) AS total_events,
  APPROX_QUANTILES(events.days_before_index, 2)[OFFSET(1)] AS median_days_before_index
FROM events
JOIN window_defs
  ON window_defs.days IS NULL OR events.days_before_index BETWEEN 1 AND window_defs.days
GROUP BY domain, window_defs.window, events.concept_id, concept_name
QUALIFY ROW_NUMBER() OVER (
  PARTITION BY window_defs.window
  ORDER BY n_cases DESC, total_events DESC, concept_name
) <= {top_n}
""".strip()


def _render_event_summary_sql(
    config: ProjectConfig,
    cases: pd.DataFrame,
    *,
    domain: str,
    windows: tuple[int | None, ...],
) -> str:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    if domain == "condition":
        event_table = _qualify_table(cdr, config.phenotype.tables.condition_table)
        date_column = config.phenotype.condition_date_column
    elif domain == "measurement":
        event_table = _qualify_table(cdr, config.phenotype.tables.measurement_table)
        date_column = config.phenotype.measurement_date_column
    else:
        raise ValueError(f"Unsupported domain: {domain}")
    n_cases = cases["person_id"].nunique()
    return f"""
WITH
{_case_cte_sql(cases)},
{_window_defs_sql(windows)},
events AS (
  SELECT
    c.person_id,
    DATE_DIFF(c.index_date, DATE(e.{date_column}), DAY) AS days_before_index
  FROM cases c
  JOIN `{event_table}` e
    ON CAST(e.{config.phenotype.person_id_column} AS STRING) = c.person_id
  WHERE DATE(e.{date_column}) < c.index_date
),
per_person AS (
  SELECT
    window_defs.window,
    events.person_id,
    COUNT(*) AS event_count
  FROM events
  JOIN window_defs
    ON window_defs.days IS NULL OR events.days_before_index BETWEEN 1 AND window_defs.days
  GROUP BY window_defs.window, events.person_id
)
SELECT
  '{domain}' AS domain,
  window_defs.window,
  {n_cases} AS n_cases,
  COUNT(per_person.person_id) AS n_cases_with_any,
  ROUND(100 * COUNT(per_person.person_id) / {n_cases}, 2) AS pct_cases_with_any,
  COALESCE(SUM(per_person.event_count), 0) AS total_events,
  APPROX_QUANTILES(per_person.event_count, 2 IGNORE NULLS)[OFFSET(1)] AS median_events_per_case_with_any
FROM window_defs
LEFT JOIN per_person
  ON window_defs.window = per_person.window
GROUP BY window_defs.window
ORDER BY window_defs.window
""".strip()


def _render_biomarker_sql(config: ProjectConfig, cases: pd.DataFrame, windows: tuple[int | None, ...]) -> str:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    measurement_table = _qualify_table(cdr, config.phenotype.tables.measurement_table)
    concept_table = _qualify_table(cdr, config.phenotype.tables.concept_table)
    term_structs = []
    for biomarker, terms in DEFAULT_BIOMARKER_TERMS.items():
        for term in terms:
            term_structs.append(f"STRUCT('{biomarker}' AS biomarker, '{term.lower()}' AS term)")
    return f"""
WITH
{_case_cte_sql(cases)},
{_window_defs_sql(windows)},
biomarker_terms AS (
  SELECT * FROM UNNEST([
    {", ".join(term_structs)}
  ])
),
events AS (
  SELECT
    c.person_id,
    concept.concept_name AS concept_name,
    DATE_DIFF(c.index_date, DATE(m.{config.phenotype.measurement_date_column}), DAY) AS days_before_index,
    SAFE_CAST(m.{config.phenotype.measurement_value_column} AS FLOAT64) AS value_as_number
  FROM cases c
  JOIN `{measurement_table}` m
    ON CAST(m.{config.phenotype.person_id_column} AS STRING) = c.person_id
  LEFT JOIN `{concept_table}` concept
    ON m.{config.phenotype.measurement_concept_column} = concept.concept_id
  WHERE DATE(m.{config.phenotype.measurement_date_column}) < c.index_date
)
SELECT
  biomarker_terms.biomarker,
  window_defs.window,
  {cases['person_id'].nunique()} AS n_cases,
  COUNT(DISTINCT events.person_id) AS n_cases_with_measurement,
  ROUND(100 * COUNT(DISTINCT events.person_id) / {cases['person_id'].nunique()}, 2) AS pct_cases_with_measurement,
  COUNT(*) AS total_measurements,
  APPROX_QUANTILES(events.value_as_number, 2 IGNORE NULLS)[OFFSET(1)] AS median_latest_value
FROM events
JOIN biomarker_terms
  ON LOWER(COALESCE(events.concept_name, '')) LIKE CONCAT('%', biomarker_terms.term, '%')
JOIN window_defs
  ON window_defs.days IS NULL OR events.days_before_index BETWEEN 1 AND window_defs.days
GROUP BY biomarker_terms.biomarker, window_defs.window
ORDER BY biomarker_terms.biomarker, window_defs.window
""".strip()


def _profile_bigquery(
    config: ProjectConfig,
    paths: ProjectPaths,
    cases: pd.DataFrame,
    windows: tuple[int | None, ...],
    top_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    condition_summary_sql = _render_event_summary_sql(config, cases, domain="condition", windows=windows)
    measurement_summary_sql = _render_event_summary_sql(config, cases, domain="measurement", windows=windows)
    condition_sql = _render_top_concepts_sql(config, cases, domain="condition", windows=windows, top_n=top_n)
    measurement_sql = _render_top_concepts_sql(config, cases, domain="measurement", windows=windows, top_n=top_n)
    biomarker_sql = _render_biomarker_sql(config, cases, windows)
    write_text(condition_summary_sql, join_path(preindex_case_profile_root(paths), "sql", "condition_summary.sql"))
    write_text(measurement_summary_sql, join_path(preindex_case_profile_root(paths), "sql", "measurement_summary.sql"))
    write_text(condition_sql, join_path(preindex_case_profile_root(paths), "sql", "top_conditions.sql"))
    write_text(measurement_sql, join_path(preindex_case_profile_root(paths), "sql", "top_measurements.sql"))
    write_text(biomarker_sql, join_path(preindex_case_profile_root(paths), "sql", "biomarker_availability.sql"))
    event_summary = pd.concat(
        [
            query_bigquery_dataframe(condition_summary_sql),
            query_bigquery_dataframe(measurement_summary_sql),
        ],
        ignore_index=True,
    )
    return (
        event_summary,
        query_bigquery_dataframe(condition_sql),
        query_bigquery_dataframe(measurement_sql),
        query_bigquery_dataframe(biomarker_sql),
    )


def _profile_local(
    config: ProjectConfig,
    cases: pd.DataFrame,
    windows: tuple[int | None, ...],
    top_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    conditions = _condition_events_local(config, cases)
    measurements = _measurement_events_local(config, cases)
    summary = pd.concat(
        [
            _summarize_event_domain(events=conditions, cases=cases, domain="condition", windows=windows),
            _summarize_event_domain(events=measurements, cases=cases, domain="measurement", windows=windows),
        ],
        ignore_index=True,
    )
    return (
        summary,
        _top_concepts(events=conditions, domain="condition", windows=windows, top_n=top_n),
        _top_concepts(events=measurements, domain="measurement", windows=windows, top_n=top_n),
        _biomarker_availability_local(measurements, cases, windows),
        measurements,
    )


def _observation_summary(cases: pd.DataFrame, windows: tuple[int | None, ...]) -> pd.DataFrame:
    n_cases = cases["person_id"].nunique()
    rows: list[dict[str, object]] = [
        {
            "domain": "observation_history",
            "window": "preindex",
            "n_cases": n_cases,
            "n_cases_with_any": int(cases["preindex_observation_days"].gt(0).sum()),
            "pct_cases_with_any": round(100 * cases["preindex_observation_days"].gt(0).sum() / n_cases, 2),
            "total_events": np.nan,
            "median_events_per_case_with_any": np.nan,
        }
    ]
    for days in windows:
        if days is None:
            continue
        count = int(cases["preindex_observation_days"].ge(days).sum())
        rows.append(
            {
                "domain": "observation_history",
                "window": f">={days}d",
                "n_cases": n_cases,
                "n_cases_with_any": count,
                "pct_cases_with_any": round(100 * count / n_cases, 2) if n_cases else 0.0,
                "total_events": np.nan,
                "median_events_per_case_with_any": np.nan,
            }
        )
    return pd.DataFrame(rows)


def _write_report(
    *,
    paths: ProjectPaths,
    cases: pd.DataFrame,
    summary: pd.DataFrame,
    condition_top: pd.DataFrame,
    measurement_top: pd.DataFrame,
    biomarker: pd.DataFrame,
    case_tier: str,
) -> None:
    obs = cases["preindex_observation_days"].dropna()
    median_obs = float(obs.median()) if not obs.empty else math.nan
    body = [
        "# Pre-Index Case Data Availability",
        "",
        f"- Case tier profiled: `{case_tier}`",
        f"- Cases with index date: {cases['person_id'].nunique()}",
        f"- Median pre-index observation days: {median_obs:.1f}" if not math.isnan(median_obs) else "- Median pre-index observation days: NA",
        "",
        "## Availability Summary",
        "",
        dataframe_markdown(summary, limit=20),
        "",
        "## Top Pre-Index Conditions",
        "",
        dataframe_markdown(condition_top, limit=15),
        "",
        "## Top Pre-Index Measurements",
        "",
        dataframe_markdown(measurement_top, limit=15),
        "",
        "## Candidate Biomarker Coverage",
        "",
        dataframe_markdown(
            biomarker.sort_values(["window", "pct_cases_with_measurement", "biomarker"], ascending=[True, False, True]),
            limit=25,
        ),
    ]
    write_text("\n".join(body), preindex_report_path(paths))


def profile_preindex_case_data(
    config: ProjectConfig,
    cohort_df: pd.DataFrame,
    paths: ProjectPaths,
    *,
    case_tier: str | None = None,
    windows: Iterable[str | int | None] | None = None,
    top_n: int = 25,
) -> dict[str, pd.DataFrame]:
    parsed_windows = _parse_windows(windows)
    tier = case_tier or config.cohort.primary_case_tier
    cases = _case_frame(config, cohort_df, tier)
    if config.phenotype.tables.cohort_table:
        event_summary, condition_top, measurement_top, biomarker, _measurements = _profile_local(
            config,
            cases,
            parsed_windows,
            top_n,
        )
        summary = pd.concat([_observation_summary(cases, parsed_windows), event_summary], ignore_index=True)
    else:
        event_summary, condition_top, measurement_top, biomarker = _profile_bigquery(
            config,
            paths,
            cases,
            parsed_windows,
            top_n,
        )
        summary = pd.concat([_observation_summary(cases, parsed_windows), event_summary], ignore_index=True)
    write_dataframe(summary, preindex_summary_path(paths))
    write_dataframe(condition_top, preindex_condition_top_path(paths))
    write_dataframe(measurement_top, preindex_measurement_top_path(paths))
    write_dataframe(biomarker, preindex_biomarker_path(paths))
    _write_report(
        paths=paths,
        cases=cases,
        summary=summary,
        condition_top=condition_top,
        measurement_top=measurement_top,
        biomarker=biomarker,
        case_tier=tier,
    )
    return {
        "summary": summary,
        "condition_top": condition_top,
        "measurement_top": measurement_top,
        "biomarker": biomarker,
    }


__all__ = [
    "DEFAULT_BIOMARKER_TERMS",
    "DEFAULT_WINDOWS",
    "preindex_biomarker_path",
    "preindex_condition_top_path",
    "preindex_measurement_top_path",
    "preindex_report_path",
    "preindex_summary_path",
    "profile_preindex_case_data",
]
