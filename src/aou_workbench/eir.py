"""Incident EIR-enriched phenotype and clinical-only prediction workflow."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np
import pandas as pd

from .clinical_model import (
    ClinicalModel,
    FeatureSpec,
    _calibration_table,
    _choose_threshold,
    _coefficient_table,
    _curve_points,
    _feature_matrix,
    _fit_logistic,
    _metrics,
    _predict,
    _roc_curve,
    _sigmoid,
    _write_calibration_svg,
    _write_line_svg,
    clinical_model_calibration_path,
    clinical_model_calibration_svg_path,
    clinical_model_coefficients_path,
    clinical_model_cv_metrics_path,
    clinical_model_metrics_path,
    clinical_model_pr_svg_path,
    clinical_model_predictions_path,
    clinical_model_report_path,
    clinical_model_roc_svg_path,
)
from .cohort import _normalize_sex, _normalize_sex_category, _prepare_ancestry_table
from .cohort_summary import (
    _binary_p_value,
    _binary_smd,
    _categorical_p_value,
    _categorical_smd,
    _continuous_p_value,
    _continuous_smd,
    _format_binary,
    _format_category_from_series,
    _format_ci,
    _format_count,
    _format_mean_sd,
    _format_p_value,
    _format_smd,
    _mean_difference_ci,
    _missing_mask,
    _odds_ratio_ci,
    _ordered_categories,
    _sex_category_series,
    _stable_hash_score,
)
from .config import CaseTierRule, ClinicalCofactorRule, ProjectConfig
from .io_utils import parse_date, query_bigquery_dataframe, read_table, write_dataframe, write_json, write_text
from .paths import ProjectPaths, build_output_paths, join_path
from .phenotype_sql import _escape_term, _match_predicate, _qualify_table
from .preflight import apply_runtime_defaults


LOOKBACK_DAYS = 365
TIME_AT_RISK_DAYS = 730
TRAUMA_EXCLUSION_START_DAYS = -7
TRAUMA_EXCLUSION_END_DAYS = 7
SEPSIS_EXCLUSION_START_DAYS = -30
SEPSIS_EXCLUSION_END_DAYS = 0
SUPPORTIVE_EIR_START_DAYS = -30
SUPPORTIVE_EIR_END_DAYS = 7
TRAIN_FRACTION = 0.8
CV_FOLDS = 5

OUTCOME_COLUMN = "eir_primary_case"
TRAUMA_COFACTORS = ("crush_injury", "major_trauma")
SEPSIS_COFACTOR = "sepsis"
SUPPORTIVE_EIR_COFATORS = ("exertion_exercise", "heat_illness", "dehydration")

LAB_TERMS: dict[str, tuple[str, ...]] = {
    "ck": ("creatine kinase",),
    "creatinine": ("creatinine",),
    "ast": ("aspartate aminotransferase", "ast"),
    "alt": ("alanine aminotransferase", "alt"),
    "potassium": ("potassium",),
    "sodium": ("sodium",),
    "calcium": ("calcium",),
    "tsh": ("thyroid stimulating hormone", "thyrotropin", "tsh"),
}


def eir_consort_counts_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "consort_counts.tsv")


def eir_consort_counts_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "consort_counts.md")


def eir_model_input_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "clinical_model_input.tsv")


def eir_split_summary_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "model_split_summary.tsv")


def eir_table1_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "table1_clinical_by_split.tsv")


def eir_missingness_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "missingness_summary.tsv")


def eir_characterization_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "cohort", "eir_characterization_report.md")


def eir_sparse_status_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "clinical", "sparse_omop_model", "status.json")


def eir_sparse_report_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "clinical", "sparse_omop_model", "report.md")


def eir_risk_decile_path(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "clinical", "model", "risk_deciles.svg")


def _copy_with_eir_outcome(config: ProjectConfig) -> ProjectConfig:
    analysis = replace(config.analysis, matched_outcome_column=OUTCOME_COLUMN)
    return replace(config, analysis=analysis)


def _empty(columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=columns)


def _md_table(frame: pd.DataFrame) -> str:
    try:
        return frame.to_markdown(index=False)
    except ImportError:
        return frame.to_string(index=False)


def _suppress_small_counts(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    count_names = {"n", "count", "n_rows", "cases", "controls", "n_cases", "n_controls", "missing_n", "observed_n"}
    for column in output.columns:
        if column.lower() in count_names:
            values = pd.to_numeric(output[column], errors="coerce")
            output[column] = output[column].astype(object)
            output.loc[values.between(1, 20, inclusive="both"), column] = "<=20"
    return output


def _write_md_table(frame: pd.DataFrame, path: str, title: str) -> None:
    write_text(f"# {title}\n\n{_md_table(_suppress_small_counts(frame))}\n", path)


def _condition_match_mask(raw: pd.DataFrame, config: ProjectConfig, rule: ClinicalCofactorRule | CaseTierRule) -> pd.Series:
    mask = pd.Series(False, index=raw.index)
    ids = tuple(getattr(rule, "condition_concept_ids", ()))
    terms = tuple(getattr(rule, "condition_terms", ()))
    if ids and config.phenotype.condition_concept_column in raw.columns:
        mask |= raw[config.phenotype.condition_concept_column].isin(ids)
    if terms and "condition_concept_name" in raw.columns:
        lowered = raw["condition_concept_name"].astype(str).str.lower()
        term_mask = pd.Series(False, index=raw.index)
        for term in terms:
            term_mask |= lowered.str.contains(term.lower(), regex=False, na=False)
        mask |= term_mask
    return mask


def _measurement_match_mask(raw: pd.DataFrame, config: ProjectConfig, terms: tuple[str, ...], ids: tuple[int, ...] = ()) -> pd.Series:
    mask = pd.Series(False, index=raw.index)
    if ids and config.phenotype.measurement_concept_column in raw.columns:
        mask |= raw[config.phenotype.measurement_concept_column].isin(ids)
    if terms and "measurement_concept_name" in raw.columns:
        lowered = raw["measurement_concept_name"].astype(str).str.lower()
        term_mask = pd.Series(False, index=raw.index)
        for term in terms:
            term_mask |= lowered.str.contains(term.lower(), regex=False, na=False)
        mask |= term_mask
    return mask


def _all_condition_dates_local(config: ProjectConfig) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.condition_table).copy()
    if raw.empty:
        return _empty(["person_id", "condition_date"])
    source_column = next(
        (
            column
            for column in ("condition_source_vocabulary_id", "source_vocabulary_id", "vocabulary_id")
            if column in raw.columns
        ),
        None,
    )
    if source_column:
        source_values = raw[source_column].astype(str).str.upper()
        if source_values.str.startswith("ICD").any():
            raw = raw[source_values.str.startswith("ICD")].copy()
    frame = pd.DataFrame(
        {
            "person_id": raw[config.phenotype.person_id_column].astype(str),
            "condition_date": parse_date(raw[config.phenotype.condition_date_column]),
        }
    )
    return frame.dropna(subset=["condition_date"]).drop_duplicates()


def _all_measurement_dates_local(config: ProjectConfig) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.measurement_table).copy()
    if raw.empty:
        return _empty(["person_id", "measurement_date"])
    frame = pd.DataFrame(
        {
            "person_id": raw[config.phenotype.person_id_column].astype(str),
            "measurement_date": parse_date(raw[config.phenotype.measurement_date_column]),
        }
    )
    return frame.dropna(subset=["measurement_date"]).drop_duplicates()


def _attach_baseline_dates(
    base: pd.DataFrame,
    condition_dates: pd.DataFrame,
    measurement_dates: pd.DataFrame,
) -> pd.DataFrame:
    output = base.copy()
    output["washout_end_date"] = output["obs_start_date"] + pd.to_timedelta(LOOKBACK_DAYS, unit="D")
    if condition_dates.empty:
        output["omop_condition_record_dates"] = 0
        output["second_condition_date"] = pd.NaT
        output["last_condition_date"] = pd.NaT
        output["baseline_date"] = pd.NaT
    else:
        grouped = condition_dates.sort_values(["person_id", "condition_date"]).groupby("person_id")["condition_date"]
        summary = grouped.agg(
            omop_condition_record_dates=lambda values: pd.Series(values).dropna().nunique(),
            second_condition_date=lambda values: sorted(pd.Series(values).dropna().unique())[1]
            if len(pd.Series(values).dropna().unique()) >= 2
            else pd.NaT,
            last_condition_date="max",
        ).reset_index()
        output = output.merge(summary, on="person_id", how="left")
        activity_frames = [
            condition_dates.rename(columns={"condition_date": "activity_date"})[["person_id", "activity_date"]]
        ]
        if not measurement_dates.empty:
            activity_frames.append(
                measurement_dates.rename(columns={"measurement_date": "activity_date"})[["person_id", "activity_date"]]
            )
        activity_dates = pd.concat(activity_frames, ignore_index=True).dropna(subset=["activity_date"]).drop_duplicates()
        candidates = activity_dates.merge(
            output[["person_id", "washout_end_date", "second_condition_date"]],
            on="person_id",
            how="inner",
        )
        candidates = candidates[
            (candidates["activity_date"] >= candidates["washout_end_date"])
            & (candidates["activity_date"] >= candidates["second_condition_date"])
        ].copy()
        baseline = candidates.groupby("person_id", as_index=False)["activity_date"].min().rename(columns={"activity_date": "baseline_date"})
        output = output.merge(baseline, on="person_id", how="left")
    output["omop_condition_record_dates"] = pd.to_numeric(output.get("omop_condition_record_dates"), errors="coerce").fillna(0).astype(int)
    if "last_condition_date" not in output.columns:
        output["last_condition_date"] = pd.NaT
    if measurement_dates.empty:
        output["last_measurement_date"] = pd.NaT
    else:
        last_measurement = measurement_dates.groupby("person_id", as_index=False)["measurement_date"].max().rename(
            columns={"measurement_date": "last_measurement_date"}
        )
        output = output.merge(last_measurement, on="person_id", how="left")
    output["horizon_end_date"] = output["baseline_date"] + pd.to_timedelta(TIME_AT_RISK_DAYS, unit="D")
    output["eligible_ehr_denominator"] = (
        output["baseline_date"].notna()
        & output["obs_start_date"].notna()
        & output["obs_end_date"].notna()
        & output["omop_condition_record_dates"].ge(2)
        & output["washout_end_date"].le(output["baseline_date"])
    )
    return output


def _baseline_local(config: ProjectConfig) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.cohort_table).copy()
    base = pd.DataFrame({"person_id": raw[config.phenotype.person_id_column].astype(str)})
    base["obs_start_date"] = parse_date(raw[config.phenotype.observation_start_column])
    base["obs_end_date"] = parse_date(raw[config.phenotype.observation_end_column])
    base["year_of_birth"] = pd.to_numeric(raw.get("year_of_birth"), errors="coerce")
    if config.phenotype.birth_date_column and config.phenotype.birth_date_column in raw.columns:
        base["birth_date"] = parse_date(raw[config.phenotype.birth_date_column])
    else:
        base["birth_date"] = pd.NaT
    base["age_raw"] = (
        pd.to_numeric(raw[config.phenotype.age_column], errors="coerce")
        if config.phenotype.age_column in raw.columns
        else np.nan
    )
    if config.phenotype.sex_column in raw.columns:
        base["gender_concept_name"] = raw[config.phenotype.sex_column].astype(str)
        base["sex_category"] = _normalize_sex_category(raw[config.phenotype.sex_column])
        base["is_female"] = _normalize_sex(raw[config.phenotype.sex_column])
    else:
        base["gender_concept_name"] = np.nan
        base["sex_category"] = "missing"
        base["is_female"] = np.nan

    base = _attach_baseline_dates(base, _all_condition_dates_local(config), _all_measurement_dates_local(config))
    ancestry = _prepare_ancestry_table(config)
    if not ancestry.empty:
        base = base.merge(ancestry, on="person_id", how="left")
    return base


def _baseline_bigquery(config: ProjectConfig) -> pd.DataFrame:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    person_table = _qualify_table(cdr, config.phenotype.tables.person_table)
    observation_table = _qualify_table(cdr, config.phenotype.tables.observation_table)
    condition_table = _qualify_table(cdr, config.phenotype.tables.condition_table)
    measurement_table = _qualify_table(cdr, config.phenotype.tables.measurement_table)
    concept_table = _qualify_table(cdr, config.phenotype.tables.concept_table)
    sql = f"""
WITH observation_windows AS (
  SELECT
    CAST(person_id AS STRING) AS person_id,
    MIN(DATE(observation_period_start_date)) AS obs_start_date,
    MAX(DATE(observation_period_end_date)) AS obs_end_date
  FROM `{observation_table}`
  GROUP BY person_id
),
condition_records AS (
  SELECT
    CAST(co.{config.phenotype.person_id_column} AS STRING) AS person_id,
    DATE(co.{config.phenotype.condition_date_column}) AS condition_date,
    UPPER(COALESCE(source_concept.vocabulary_id, '')) AS source_vocabulary_id
  FROM `{condition_table}` co
  LEFT JOIN `{concept_table}` source_concept
    ON co.condition_source_concept_id = source_concept.concept_id
  WHERE co.{config.phenotype.condition_date_column} IS NOT NULL
),
denominator_mode AS (
  SELECT
    IF(COUNTIF(STARTS_WITH(source_vocabulary_id, 'ICD')) > 0, 'icd_source_condition_records', 'all_condition_records')
      AS omop_condition_record_source
  FROM condition_records
),
filtered_condition_records AS (
  SELECT person_id, condition_date
  FROM condition_records
  CROSS JOIN denominator_mode
  WHERE denominator_mode.omop_condition_record_source = 'all_condition_records'
    OR STARTS_WITH(condition_records.source_vocabulary_id, 'ICD')
  GROUP BY person_id, condition_date
),
condition_rollup AS (
  SELECT
    person_id,
    COUNT(DISTINCT condition_date) AS omop_condition_record_dates,
    ARRAY_AGG(DISTINCT condition_date IGNORE NULLS ORDER BY condition_date) AS condition_dates,
    MAX(condition_date) AS last_condition_date
  FROM filtered_condition_records
  GROUP BY person_id
),
measurement_rollup AS (
  SELECT
    CAST({config.phenotype.person_id_column} AS STRING) AS person_id,
    MAX(DATE({config.phenotype.measurement_date_column})) AS last_measurement_date
  FROM `{measurement_table}`
  WHERE {config.phenotype.measurement_date_column} IS NOT NULL
  GROUP BY person_id
),
activity_records AS (
  SELECT person_id, condition_date AS activity_date
  FROM filtered_condition_records
  UNION DISTINCT
  SELECT
    CAST({config.phenotype.person_id_column} AS STRING) AS person_id,
    DATE({config.phenotype.measurement_date_column}) AS activity_date
  FROM `{measurement_table}`
  WHERE {config.phenotype.measurement_date_column} IS NOT NULL
)
SELECT
  CAST(p.person_id AS STRING) AS person_id,
  observation_windows.obs_start_date,
  observation_windows.obs_end_date,
  DATE_ADD(observation_windows.obs_start_date, INTERVAL {LOOKBACK_DAYS} DAY) AS washout_end_date,
  SAFE_CAST(p.year_of_birth AS FLOAT64) AS year_of_birth,
  SAFE_CAST(EXTRACT(YEAR FROM CURRENT_DATE()) - p.year_of_birth AS FLOAT64) AS age_raw,
  CAST(NULL AS DATE) AS birth_date,
  g.concept_name AS gender_concept_name,
  CASE
    WHEN LOWER(g.concept_name) = 'female' THEN 'female'
    WHEN LOWER(g.concept_name) = 'male' THEN 'male'
    WHEN g.concept_name IS NULL
      OR LOWER(g.concept_name) IN ('', 'no matching concept', 'none', 'unknown', 'skip', 'pmi: skip', 'prefer not to answer', 'pmi: prefer not to answer')
      THEN 'missing'
    ELSE 'other_or_unknown'
  END AS sex_category,
  CASE WHEN LOWER(g.concept_name) = 'female' THEN 1.0 WHEN LOWER(g.concept_name) = 'male' THEN 0.0 ELSE NULL END AS is_female,
  CAST(NULL AS STRING) AS ancestry_pred,
  CAST(NULL AS FLOAT64) AS pc1,
  CAST(NULL AS FLOAT64) AS pc2,
  CAST(NULL AS FLOAT64) AS pc3,
  CAST(NULL AS FLOAT64) AS pc4,
  CAST(NULL AS FLOAT64) AS pc5,
  CAST(NULL AS FLOAT64) AS pc6,
  CAST(NULL AS FLOAT64) AS pc7,
  CAST(NULL AS FLOAT64) AS pc8,
  CAST(NULL AS FLOAT64) AS pc9,
  CAST(NULL AS FLOAT64) AS pc10,
  COALESCE(condition_rollup.omop_condition_record_dates, 0) AS omop_condition_record_dates,
  denominator_mode.omop_condition_record_source,
  condition_rollup.condition_dates[SAFE_OFFSET(1)] AS second_condition_date,
  condition_rollup.last_condition_date,
  measurement_rollup.last_measurement_date,
  (
    SELECT MIN(activity_date)
    FROM activity_records
    WHERE activity_records.person_id = CAST(p.person_id AS STRING)
      AND activity_date >= DATE_ADD(observation_windows.obs_start_date, INTERVAL {LOOKBACK_DAYS} DAY)
      AND activity_date >= condition_rollup.condition_dates[SAFE_OFFSET(1)]
  ) AS baseline_date
FROM `{person_table}` p
LEFT JOIN observation_windows
  ON CAST(p.person_id AS STRING) = observation_windows.person_id
LEFT JOIN `{concept_table}` g
  ON p.gender_concept_id = g.concept_id
LEFT JOIN condition_rollup
  ON CAST(p.person_id AS STRING) = condition_rollup.person_id
LEFT JOIN measurement_rollup
  ON CAST(p.person_id AS STRING) = measurement_rollup.person_id
CROSS JOIN denominator_mode
""".strip()
    frame = query_bigquery_dataframe(sql)
    for column in (
        "obs_start_date",
        "obs_end_date",
        "washout_end_date",
        "birth_date",
        "second_condition_date",
        "last_condition_date",
        "last_measurement_date",
        "baseline_date",
    ):
        frame[column] = parse_date(frame[column])
    frame["horizon_end_date"] = frame["baseline_date"] + pd.to_timedelta(TIME_AT_RISK_DAYS, unit="D")
    frame["omop_condition_record_dates"] = pd.to_numeric(frame["omop_condition_record_dates"], errors="coerce").fillna(0).astype(int)
    frame["eligible_ehr_denominator"] = (
        frame["baseline_date"].notna()
        & frame["obs_start_date"].notna()
        & frame["obs_end_date"].notna()
        & frame["omop_condition_record_dates"].ge(2)
        & frame["washout_end_date"].le(frame["baseline_date"])
    )
    ancestry = _prepare_ancestry_table(config)
    if not ancestry.empty:
        frame = frame.drop(columns=[column for column in ancestry.columns if column != "person_id" and column in frame.columns])
        frame = frame.merge(ancestry, on="person_id", how="left")
    return frame


def _condition_events_local(
    config: ProjectConfig,
    rule: ClinicalCofactorRule | CaseTierRule,
    label: str,
) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.condition_table).copy()
    if raw.empty:
        return _empty(["person_id", "event", "event_date", "concept_name"])
    filtered = raw[_condition_match_mask(raw, config, rule)].copy()
    if filtered.empty:
        return _empty(["person_id", "event", "event_date", "concept_name"])
    return pd.DataFrame(
        {
            "person_id": filtered[config.phenotype.person_id_column].astype(str),
            "event": label,
            "event_date": parse_date(filtered[config.phenotype.condition_date_column]),
            "concept_name": filtered.get("condition_concept_name", pd.Series("", index=filtered.index)).astype(str),
        }
    ).dropna(subset=["event_date"]).drop_duplicates()


def _condition_events_bigquery(
    config: ProjectConfig,
    rule: ClinicalCofactorRule | CaseTierRule,
    label: str,
) -> pd.DataFrame:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    condition_table = _qualify_table(cdr, config.phenotype.tables.condition_table)
    concept_table = _qualify_table(cdr, config.phenotype.tables.concept_table)
    predicate = _match_predicate(
        concept_id_column=f"co.{config.phenotype.condition_concept_column}",
        concept_name_column="condition_concept.concept_name",
        concept_ids=tuple(getattr(rule, "condition_concept_ids", ())),
        concept_terms=tuple(getattr(rule, "condition_terms", ())),
    )
    if predicate == "FALSE":
        return _empty(["person_id", "event", "event_date", "concept_name"])
    sql = f"""
SELECT
  CAST(co.{config.phenotype.person_id_column} AS STRING) AS person_id,
  '{_escape_term(label)}' AS event,
  DATE(co.{config.phenotype.condition_date_column}) AS event_date,
  condition_concept.concept_name AS concept_name
FROM `{condition_table}` co
LEFT JOIN `{concept_table}` condition_concept
  ON co.{config.phenotype.condition_concept_column} = condition_concept.concept_id
WHERE ({predicate})
  AND co.{config.phenotype.condition_date_column} IS NOT NULL
GROUP BY person_id, event, event_date, concept_name
""".strip()
    frame = query_bigquery_dataframe(sql)
    if frame.empty:
        return _empty(["person_id", "event", "event_date", "concept_name"])
    frame["event_date"] = parse_date(frame["event_date"])
    return frame


def _ck_measurement_events_local(config: ProjectConfig, *, high_only: bool = False) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.measurement_table).copy()
    if raw.empty:
        return _empty(["person_id", "measurement_date", "measurement_value"])
    rule = config.phenotype.definite
    mask = _measurement_match_mask(raw, config, rule.measurement_terms or LAB_TERMS["ck"], rule.measurement_concept_ids)
    filtered = raw[mask].copy()
    filtered["measurement_value"] = pd.to_numeric(filtered[config.phenotype.measurement_value_column], errors="coerce")
    if high_only and rule.measurement_min is not None:
        filtered = filtered[filtered["measurement_value"] >= rule.measurement_min].copy()
    if filtered.empty:
        return _empty(["person_id", "measurement_date", "measurement_value"])
    return pd.DataFrame(
        {
            "person_id": filtered[config.phenotype.person_id_column].astype(str),
            "measurement_date": parse_date(filtered[config.phenotype.measurement_date_column]),
            "measurement_value": filtered["measurement_value"],
        }
    ).dropna(subset=["measurement_date"]).drop_duplicates()


def _ck_measurement_events_bigquery(config: ProjectConfig, *, high_only: bool = False) -> pd.DataFrame:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    measurement_table = _qualify_table(cdr, config.phenotype.tables.measurement_table)
    concept_table = _qualify_table(cdr, config.phenotype.tables.concept_table)
    rule = config.phenotype.definite
    predicate = _match_predicate(
        concept_id_column=f"m.{config.phenotype.measurement_concept_column}",
        concept_name_column="measurement_concept.concept_name",
        concept_ids=rule.measurement_concept_ids,
        concept_terms=rule.measurement_terms or LAB_TERMS["ck"],
    )
    if predicate == "FALSE":
        return _empty(["person_id", "measurement_date", "measurement_value"])
    threshold = (
        f"AND SAFE_CAST(m.{config.phenotype.measurement_value_column} AS FLOAT64) >= {rule.measurement_min}"
        if high_only and rule.measurement_min is not None
        else ""
    )
    sql = f"""
SELECT
  CAST(m.{config.phenotype.person_id_column} AS STRING) AS person_id,
  DATE(m.{config.phenotype.measurement_date_column}) AS measurement_date,
  SAFE_CAST(m.{config.phenotype.measurement_value_column} AS FLOAT64) AS measurement_value
FROM `{measurement_table}` m
LEFT JOIN `{concept_table}` measurement_concept
  ON m.{config.phenotype.measurement_concept_column} = measurement_concept.concept_id
WHERE ({predicate})
  AND m.{config.phenotype.measurement_date_column} IS NOT NULL
  {threshold}
""".strip()
    frame = query_bigquery_dataframe(sql)
    if frame.empty:
        return _empty(["person_id", "measurement_date", "measurement_value"])
    frame["measurement_date"] = parse_date(frame["measurement_date"])
    return frame


def _selected_lab_events_local(config: ProjectConfig) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.measurement_table).copy()
    if raw.empty:
        return _empty(["person_id", "lab", "measurement_date", "measurement_value"])
    rows: list[pd.DataFrame] = []
    for lab, terms in LAB_TERMS.items():
        ids = config.phenotype.definite.measurement_concept_ids if lab == "ck" else ()
        filtered = raw[_measurement_match_mask(raw, config, terms, ids)].copy()
        if filtered.empty:
            continue
        rows.append(
            pd.DataFrame(
                {
                    "person_id": filtered[config.phenotype.person_id_column].astype(str),
                    "lab": lab,
                    "measurement_date": parse_date(filtered[config.phenotype.measurement_date_column]),
                    "measurement_value": pd.to_numeric(filtered[config.phenotype.measurement_value_column], errors="coerce"),
                }
            )
        )
    if not rows:
        return _empty(["person_id", "lab", "measurement_date", "measurement_value"])
    return pd.concat(rows, ignore_index=True).dropna(subset=["measurement_date", "measurement_value"])


def _selected_lab_events_bigquery(config: ProjectConfig) -> pd.DataFrame:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    measurement_table = _qualify_table(cdr, config.phenotype.tables.measurement_table)
    concept_table = _qualify_table(cdr, config.phenotype.tables.concept_table)
    selects: list[str] = []
    for lab, terms in LAB_TERMS.items():
        predicate = _match_predicate(
            concept_id_column=f"m.{config.phenotype.measurement_concept_column}",
            concept_name_column="measurement_concept.concept_name",
            concept_ids=config.phenotype.definite.measurement_concept_ids if lab == "ck" else (),
            concept_terms=terms,
        )
        if predicate == "FALSE":
            continue
        selects.append(
            f"""
SELECT
  CAST(m.{config.phenotype.person_id_column} AS STRING) AS person_id,
  '{_escape_term(lab)}' AS lab,
  DATE(m.{config.phenotype.measurement_date_column}) AS measurement_date,
  SAFE_CAST(m.{config.phenotype.measurement_value_column} AS FLOAT64) AS measurement_value
FROM `{measurement_table}` m
LEFT JOIN `{concept_table}` measurement_concept
  ON m.{config.phenotype.measurement_concept_column} = measurement_concept.concept_id
WHERE ({predicate})
  AND m.{config.phenotype.measurement_date_column} IS NOT NULL
  AND SAFE_CAST(m.{config.phenotype.measurement_value_column} AS FLOAT64) IS NOT NULL
""".strip()
        )
    if not selects:
        return _empty(["person_id", "lab", "measurement_date", "measurement_value"])
    frame = query_bigquery_dataframe("\nUNION ALL\n".join(selects))
    if frame.empty:
        return _empty(["person_id", "lab", "measurement_date", "measurement_value"])
    frame["measurement_date"] = parse_date(frame["measurement_date"])
    return frame


def _cofactor_events(config: ProjectConfig) -> pd.DataFrame:
    rows = []
    for rule in config.phenotype.clinical_cofactors:
        event_frame = (
            _condition_events_local(config, rule, rule.name)
            if config.phenotype.tables.cohort_table
            else _condition_events_bigquery(config, rule, rule.name)
        )
        if not event_frame.empty:
            rows.append(event_frame)
    if not rows:
        return _empty(["person_id", "event", "event_date", "concept_name"])
    return pd.concat(rows, ignore_index=True).drop_duplicates()


def _first_in_horizon(events: pd.DataFrame, cohort: pd.DataFrame, event_name: str) -> pd.Series:
    subset = events[events["event"] == event_name].merge(
        cohort[["person_id", "baseline_date", "horizon_end_date"]],
        on="person_id",
        how="inner",
    )
    subset = subset[(subset["event_date"] > subset["baseline_date"]) & (subset["event_date"] <= subset["horizon_end_date"])].copy()
    if subset.empty:
        return pd.Series(pd.NaT, index=cohort["person_id"].astype(str), dtype="datetime64[ns]")
    first = subset.groupby("person_id")["event_date"].min()
    return cohort["person_id"].astype(str).map(first)


def _has_event_between(
    events: pd.DataFrame,
    cohort: pd.DataFrame,
    event_names: tuple[str, ...],
    anchor_column: str,
    start_days: int,
    end_days: int,
) -> pd.Series:
    if events.empty or not event_names:
        return pd.Series(False, index=cohort.index)
    indexed = cohort[["person_id", anchor_column]].copy()
    indexed["_row_id"] = cohort.index
    subset = events[events["event"].isin(event_names)].merge(indexed, on="person_id", how="inner")
    subset = subset.dropna(subset=[anchor_column, "event_date"]).copy()
    if subset.empty:
        return pd.Series(False, index=cohort.index)
    delta = (subset["event_date"] - subset[anchor_column]).dt.days
    hits = subset.loc[delta.between(start_days, end_days, inclusive="both"), "_row_id"].drop_duplicates()
    out = pd.Series(False, index=cohort.index)
    out.loc[hits] = True
    return out


def _has_measurement_between(
    measurements: pd.DataFrame,
    cohort: pd.DataFrame,
    date_column: str,
    start_column: str,
    end_column: str,
) -> pd.Series:
    if measurements.empty:
        return pd.Series(False, index=cohort.index)
    indexed = cohort[["person_id", start_column, end_column]].copy()
    indexed["_row_id"] = cohort.index
    subset = measurements.merge(indexed, on="person_id", how="inner")
    subset = subset.dropna(subset=[date_column, start_column, end_column]).copy()
    hits = subset.loc[
        (subset[date_column] >= subset[start_column])
        & (subset[date_column] <= subset[end_column]),
        "_row_id",
    ].drop_duplicates()
    out = pd.Series(False, index=cohort.index)
    out.loc[hits] = True
    return out


def _first_ck_near_index(ck: pd.DataFrame, cohort: pd.DataFrame, config: ProjectConfig) -> pd.DataFrame:
    output = cohort.copy()
    output["ck_confirming_date"] = pd.NaT
    output["ck_confirming_value"] = np.nan
    if ck.empty:
        return output
    indexed = cohort[["person_id", "first_rhabdo_date"]].copy()
    indexed["_row_id"] = cohort.index
    subset = ck.merge(indexed, on="person_id", how="inner")
    subset = subset.dropna(subset=["measurement_date", "first_rhabdo_date"]).copy()
    delta = (subset["measurement_date"] - subset["first_rhabdo_date"]).dt.days
    subset = subset[
        delta.between(
            config.phenotype.definite.measurement_window_start_days,
            config.phenotype.definite.measurement_window_end_days,
            inclusive="both",
        )
    ].copy()
    if subset.empty:
        return output
    subset["_abs_delta"] = (subset["measurement_date"] - subset["first_rhabdo_date"]).dt.days.abs()
    first = subset.sort_values(["_row_id", "_abs_delta", "measurement_date"]).groupby("_row_id", as_index=False).first()
    output.loc[first["_row_id"], "ck_confirming_date"] = first["measurement_date"].to_numpy()
    output.loc[first["_row_id"], "ck_confirming_value"] = first["measurement_value"].to_numpy()
    return output


def _age_at_baseline(frame: pd.DataFrame) -> pd.Series:
    age = pd.to_numeric(frame.get("age_raw"), errors="coerce")
    birth_date = parse_date(frame.get("birth_date", pd.Series(pd.NaT, index=frame.index)))
    baseline = parse_date(frame["baseline_date"])
    has_birth = birth_date.notna() & baseline.notna()
    age = age.copy()
    age.loc[has_birth] = (baseline.loc[has_birth] - birth_date.loc[has_birth]).dt.days / 365.25
    yob = pd.to_numeric(frame.get("year_of_birth"), errors="coerce")
    has_yob = yob.notna() & baseline.notna()
    age.loc[has_yob] = baseline.loc[has_yob].dt.year - yob.loc[has_yob]
    return age


def _first_prior_event_date(events: pd.DataFrame, cohort: pd.DataFrame, event_name: str) -> pd.Series:
    if events.empty:
        return pd.Series(pd.NaT, index=cohort.index, dtype="datetime64[ns]")
    indexed = cohort[["person_id", "baseline_date"]].copy()
    indexed["_row_id"] = cohort.index
    subset = events[events["event"] == event_name].merge(indexed, on="person_id", how="inner")
    subset = subset.dropna(subset=["event_date", "baseline_date"]).copy()
    subset = subset[subset["event_date"] < subset["baseline_date"]].copy()
    out = pd.Series(pd.NaT, index=cohort.index, dtype="datetime64[ns]")
    if subset.empty:
        return out
    first = subset.groupby("_row_id")["event_date"].min()
    out.loc[first.index] = first.to_numpy()
    return out


def _has_measurement_relative_to_baseline(
    measurements: pd.DataFrame,
    cohort: pd.DataFrame,
    start_days: int,
    end_days: int,
) -> pd.Series:
    if measurements.empty:
        return pd.Series(False, index=cohort.index)
    indexed = cohort[["person_id", "baseline_date"]].copy()
    indexed["_row_id"] = cohort.index
    subset = measurements.merge(indexed, on="person_id", how="inner")
    subset = subset.dropna(subset=["measurement_date", "baseline_date"]).copy()
    if subset.empty:
        return pd.Series(False, index=cohort.index)
    delta = (subset["measurement_date"] - subset["baseline_date"]).dt.days
    hits = subset.loc[delta.between(start_days, end_days, inclusive="both"), "_row_id"].drop_duplicates()
    out = pd.Series(False, index=cohort.index)
    out.loc[hits] = True
    return out


def _add_prebaseline_features(
    cohort: pd.DataFrame,
    cofactor_events: pd.DataFrame,
    lab_events: pd.DataFrame,
) -> pd.DataFrame:
    output = cohort.copy()
    feature_events = {
        "prior_kidney_disease_or_aki": ("kidney_disease", "renal_injury"),
        "remote_sepsis_history": ("sepsis",),
        "prior_heat_illness_or_dehydration": ("heat_illness", "dehydration"),
        "prior_exertion_or_exercise_code": ("exertion_exercise",),
        "prior_myopathy_muscle_disease": ("myopathy_muscle_disease",),
        "prior_sickle_trait_or_disease": ("sickle_cell_trait_or_disease",),
        "prior_diabetes": ("diabetes",),
        "prior_thyroid_disease": ("thyroid_disease",),
        "prior_liver_disease": ("liver_disease",),
        "prior_alcohol_substance_condition": ("alcohol_substance",),
        "prior_injury_trauma_history": ("prior_injury_trauma", "major_trauma", "crush_injury"),
    }
    for column, event_names in feature_events.items():
        output[column] = _has_event_between(
            cofactor_events,
            output,
            event_names,
            "baseline_date",
            -100_000,
            -1,
        ).astype(int)

    output["ck_tested_prebaseline"] = 0
    if not lab_events.empty:
        indexed = output[["person_id", "baseline_date"]].copy()
        indexed["_row_id"] = output.index
        subset = lab_events.merge(indexed, on="person_id", how="inner")
        subset = subset.dropna(subset=["measurement_date", "baseline_date", "measurement_value"]).copy()
        subset = subset[subset["measurement_date"] < subset["baseline_date"]].copy()
        if not subset.empty:
            output.loc[subset.loc[subset["lab"] == "ck", "_row_id"].drop_duplicates(), "ck_tested_prebaseline"] = 1
            grouped = (
                subset.sort_values(["_row_id", "lab", "measurement_date"])
                .groupby(["_row_id", "lab"], as_index=False)
                .agg(
                    last_value=("measurement_value", "last"),
                    max_value=("measurement_value", "max"),
                    n=("measurement_value", "count"),
                )
            )
            for _, row in grouped.iterrows():
                lab = str(row["lab"])
                row_id = int(row["_row_id"])
                output.loc[row_id, f"prebaseline_{lab}_last"] = row["last_value"]
                output.loc[row_id, f"prebaseline_{lab}_max"] = row["max_value"]
                output.loc[row_id, f"prebaseline_{lab}_n"] = row["n"]
    for lab in LAB_TERMS:
        for suffix in ("last", "max", "n"):
            column = f"prebaseline_{lab}_{suffix}"
            if column not in output.columns:
                output[column] = np.nan
    return output


def build_eir_cohort(config: ProjectConfig) -> pd.DataFrame:
    """Build the incident EIR-enriched cohort with auditable phenotype flags."""

    effective = _copy_with_eir_outcome(apply_runtime_defaults(config))
    local_mode = bool(effective.phenotype.tables.cohort_table)
    baseline = _baseline_local(effective) if local_mode else _baseline_bigquery(effective)
    baseline["person_id"] = baseline["person_id"].astype(str)
    baseline["age_at_baseline"] = _age_at_baseline(baseline)

    rhabdo_events = (
        _condition_events_local(effective, effective.phenotype.broad, "rhabdomyolysis")
        if local_mode
        else _condition_events_bigquery(effective, effective.phenotype.broad, "rhabdomyolysis")
    )
    cofactor_events = _cofactor_events(effective)
    ck_all = _ck_measurement_events_local(effective, high_only=False) if local_mode else _ck_measurement_events_bigquery(effective, high_only=False)
    high_ck = _ck_measurement_events_local(effective, high_only=True) if local_mode else _ck_measurement_events_bigquery(effective, high_only=True)
    lab_events = _selected_lab_events_local(effective) if local_mode else _selected_lab_events_bigquery(effective)

    cohort = baseline.copy()
    cohort["horizon_start_date"] = cohort["baseline_date"] + pd.to_timedelta(1, unit="D")
    cohort["prior_rhabdo_date"] = _first_prior_event_date(rhabdo_events, cohort, "rhabdomyolysis")
    cohort["prior_high_ck"] = _has_measurement_relative_to_baseline(high_ck, cohort, -100_000, -1).astype(int)
    cohort["incident_denominator"] = (
        cohort["eligible_ehr_denominator"]
        & cohort["prior_rhabdo_date"].isna()
        & pd.to_numeric(cohort["prior_high_ck"], errors="coerce").fillna(0).eq(0)
    )

    cohort["first_rhabdo_date"] = _first_in_horizon(rhabdo_events, cohort, "rhabdomyolysis")
    cohort = _first_ck_near_index(ck_all, cohort, effective)
    ck_threshold = effective.phenotype.definite.measurement_min or 5000
    cohort["ck_confirmed_index"] = pd.to_numeric(cohort["ck_confirming_value"], errors="coerce").ge(ck_threshold).astype(int)
    cohort["rhabdo_during_horizon"] = cohort["first_rhabdo_date"].notna().astype(int)
    cohort["high_ck_during_horizon"] = _has_measurement_between(
        high_ck,
        cohort,
        "measurement_date",
        "horizon_start_date",
        "horizon_end_date",
    ).astype(int)
    cohort["ck_tested_during_horizon"] = _has_measurement_between(
        ck_all,
        cohort,
        "measurement_date",
        "horizon_start_date",
        "horizon_end_date",
    ).astype(int)

    cohort["excluded_periindex_trauma"] = _has_event_between(
        cofactor_events,
        cohort,
        TRAUMA_COFACTORS,
        "first_rhabdo_date",
        TRAUMA_EXCLUSION_START_DAYS,
        TRAUMA_EXCLUSION_END_DAYS,
    ).astype(int)
    cohort["excluded_pre_or_same_day_sepsis"] = _has_event_between(
        cofactor_events,
        cohort,
        (SEPSIS_COFACTOR,),
        "first_rhabdo_date",
        SEPSIS_EXCLUSION_START_DAYS,
        SEPSIS_EXCLUSION_END_DAYS,
    ).astype(int)
    cohort["supportive_exertion_heat_dehydration_code"] = _has_event_between(
        cofactor_events,
        cohort,
        SUPPORTIVE_EIR_COFATORS,
        "first_rhabdo_date",
        SUPPORTIVE_EIR_START_DAYS,
        SUPPORTIVE_EIR_END_DAYS,
    ).astype(int)

    cohort["eir_diagnosis_only_case"] = (
        cohort["incident_denominator"]
        & cohort["first_rhabdo_date"].notna()
        & pd.to_numeric(cohort["excluded_periindex_trauma"], errors="coerce").fillna(0).eq(0)
        & pd.to_numeric(cohort["excluded_pre_or_same_day_sepsis"], errors="coerce").fillna(0).eq(0)
    ).astype(int)
    cohort["eir_primary_case"] = (
        pd.to_numeric(cohort["eir_diagnosis_only_case"], errors="coerce").fillna(0).eq(1)
        & pd.to_numeric(cohort["ck_confirmed_index"], errors="coerce").fillna(0).eq(1)
    ).astype(int)
    cohort["eir_exertion_heat_dehydration_coded_case"] = (
        pd.to_numeric(cohort["eir_primary_case"], errors="coerce").fillna(0).eq(1)
        & pd.to_numeric(cohort["supportive_exertion_heat_dehydration_code"], errors="coerce").fillna(0).eq(1)
    ).astype(int)
    cohort["ck_only_no_rhabdo"] = (
        cohort["incident_denominator"]
        & pd.to_numeric(cohort["rhabdo_during_horizon"], errors="coerce").fillna(0).eq(0)
        & pd.to_numeric(cohort["high_ck_during_horizon"], errors="coerce").fillna(0).eq(1)
    ).astype(int)

    clinical_activity = pd.concat(
        [
            parse_date(cohort.get("last_condition_date", pd.Series(pd.NaT, index=cohort.index))),
            parse_date(cohort.get("last_measurement_date", pd.Series(pd.NaT, index=cohort.index))),
        ],
        axis=1,
    )
    cohort["last_clinical_activity_date"] = clinical_activity.max(axis=1)
    cohort["observed_followup_through_horizon"] = (
        cohort["horizon_end_date"].notna()
        & cohort["obs_end_date"].notna()
        & (cohort["obs_end_date"] >= cohort["horizon_end_date"])
        & (
            cohort["last_clinical_activity_date"].isna()
            | (cohort["last_clinical_activity_date"] >= cohort["baseline_date"])
        )
    ).astype(int)
    cohort["eligible_control"] = (
        cohort["incident_denominator"]
        & pd.to_numeric(cohort["observed_followup_through_horizon"], errors="coerce").fillna(0).eq(1)
        & pd.to_numeric(cohort["rhabdo_during_horizon"], errors="coerce").fillna(0).eq(0)
        & pd.to_numeric(cohort["high_ck_during_horizon"], errors="coerce").fillna(0).eq(0)
    ).astype(int)
    cohort["eir_ck_tested_control"] = (
        pd.to_numeric(cohort["eligible_control"], errors="coerce").fillna(0).eq(1)
        & pd.to_numeric(cohort["ck_tested_during_horizon"], errors="coerce").fillna(0).eq(1)
    ).astype(int)
    cohort["case_control_eligible"] = (
        pd.to_numeric(cohort["eir_primary_case"], errors="coerce").fillna(0).eq(1)
        | pd.to_numeric(cohort["eligible_control"], errors="coerce").fillna(0).eq(1)
    ).astype(int)
    cohort["analysis_case"] = pd.to_numeric(cohort["eir_primary_case"], errors="coerce").fillna(0).astype(int)
    cohort["case_tier"] = np.select(
        [
            pd.to_numeric(cohort["eir_primary_case"], errors="coerce").fillna(0).eq(1),
            pd.to_numeric(cohort["eir_diagnosis_only_case"], errors="coerce").fillna(0).eq(1),
            pd.to_numeric(cohort["ck_only_no_rhabdo"], errors="coerce").fillna(0).eq(1),
            pd.to_numeric(cohort["eligible_control"], errors="coerce").fillna(0).eq(1),
        ],
        ["eir_primary", "diagnosis_only", "ck_only_no_rhabdo", "control"],
        default="excluded_or_unclassified",
    )

    cohort = _add_prebaseline_features(cohort, cofactor_events, lab_events)
    return cohort.sort_values("person_id").reset_index(drop=True)


def _eir_consort_counts(cohort: pd.DataFrame) -> pd.DataFrame:
    rows = [
        ("Source cohort rows", len(cohort)),
        ("Observation-period rows", int(cohort["obs_start_date"].notna().sum())),
        (f">={LOOKBACK_DAYS}d lookback and >=2 condition dates by baseline", int(cohort["eligible_ehr_denominator"].sum())),
        ("Excluded prior rhabdo before baseline", int(cohort["prior_rhabdo_date"].notna().sum())),
        ("Excluded prior CK >=5000 before baseline", int(pd.to_numeric(cohort["prior_high_ck"], errors="coerce").fillna(0).sum())),
        ("Incident denominator", int(cohort["incident_denominator"].sum())),
        ("Rhabdo diagnosis in 2-year horizon", int(pd.to_numeric(cohort["rhabdo_during_horizon"], errors="coerce").fillna(0).sum())),
        ("Diagnosis plus CK >=5000 in [-7,+45]d", int(pd.to_numeric(cohort["ck_confirmed_index"], errors="coerce").fillna(0).sum())),
        ("Excluded peri-index trauma [-7,+7]d", int(pd.to_numeric(cohort["excluded_periindex_trauma"], errors="coerce").fillna(0).sum())),
        ("Excluded pre/same-day sepsis [-30,0]d", int(pd.to_numeric(cohort["excluded_pre_or_same_day_sepsis"], errors="coerce").fillna(0).sum())),
        ("Primary EIR-enriched cases", int(pd.to_numeric(cohort["eir_primary_case"], errors="coerce").fillna(0).sum())),
        ("Diagnosis-only sensitivity cases", int(pd.to_numeric(cohort["eir_diagnosis_only_case"], errors="coerce").fillna(0).sum())),
        ("Primary cases with exertion/heat/dehydration code", int(pd.to_numeric(cohort["eir_exertion_heat_dehydration_coded_case"], errors="coerce").fillna(0).sum())),
        ("CK >=5000 without rhabdo during horizon", int(pd.to_numeric(cohort["ck_only_no_rhabdo"], errors="coerce").fillna(0).sum())),
        ("Eligible controls", int(pd.to_numeric(cohort["eligible_control"], errors="coerce").fillna(0).sum())),
        ("Eligible controls with CK measured during horizon", int(pd.to_numeric(cohort["eir_ck_tested_control"], errors="coerce").fillna(0).sum())),
        ("Primary model-eligible rows", int(pd.to_numeric(cohort["case_control_eligible"], errors="coerce").fillna(0).sum())),
    ]
    return pd.DataFrame(rows, columns=["step", "n"])


def _assign_splits(analysis: pd.DataFrame) -> pd.DataFrame:
    output = analysis.copy()
    output["analysis_split"] = ""
    output["cv_fold"] = np.nan
    for value in (0, 1):
        indices = output.index[pd.to_numeric(output["analysis_case"], errors="coerce").fillna(-1).eq(value)].tolist()
        if not indices:
            continue
        ranked = sorted(indices, key=lambda idx: _stable_hash_score(output.loc[idx, "person_id"]))
        n_test = int(round(len(ranked) * (1.0 - TRAIN_FRACTION)))
        if len(ranked) >= 2:
            n_test = max(1, min(n_test, len(ranked) - 1))
        else:
            n_test = 0
        test_indices = set(ranked[:n_test])
        output.loc[list(test_indices), "analysis_split"] = "test"
        train_indices = [idx for idx in ranked if idx not in test_indices]
        output.loc[train_indices, "analysis_split"] = "train"
        for fold_index, idx in enumerate(train_indices):
            output.loc[idx, "cv_fold"] = (fold_index % CV_FOLDS) + 1
    return output


def _split_summary(model_input: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, frame in [("all", model_input), *[(split, model_input[model_input["analysis_split"] == split]) for split in ("train", "test")]]:
        rows.append(
            {
                "group": label,
                "n_rows": len(frame),
                "cases": int(pd.to_numeric(frame["analysis_case"], errors="coerce").fillna(0).sum()),
                "controls": int((pd.to_numeric(frame["analysis_case"], errors="coerce").fillna(0) == 0).sum()),
            }
        )
    for fold in sorted(pd.to_numeric(model_input["cv_fold"], errors="coerce").dropna().astype(int).unique().tolist()):
        frame = model_input[pd.to_numeric(model_input["cv_fold"], errors="coerce") == fold]
        rows.append(
            {
                "group": f"train_cv_fold_{fold}",
                "n_rows": len(frame),
                "cases": int(pd.to_numeric(frame["analysis_case"], errors="coerce").fillna(0).sum()),
                "controls": int((pd.to_numeric(frame["analysis_case"], errors="coerce").fillna(0) == 0).sum()),
            }
        )
    return pd.DataFrame(rows)


def _analysis_model_input(cohort: pd.DataFrame) -> pd.DataFrame:
    analysis = cohort[pd.to_numeric(cohort["case_control_eligible"], errors="coerce").fillna(0).eq(1)].copy()
    analysis = _assign_splits(analysis)
    analysis["model_sex_category"] = _sex_category_series(analysis)
    if "ancestry_pred" in analysis.columns:
        values = analysis["ancestry_pred"].astype("string").str.lower().str.strip()
        analysis["model_ancestry_pred"] = values.where(values.notna() & (values != ""), "missing").astype(str)
    else:
        analysis["model_ancestry_pred"] = "missing"
    return analysis.reset_index(drop=True)


def _table1_rows_for_split(model_input: pd.DataFrame, split: str) -> list[dict[str, object]]:
    frame = model_input[model_input["analysis_split"] == split].copy()
    cases = frame[pd.to_numeric(frame["analysis_case"], errors="coerce").fillna(0).eq(1)]
    controls = frame[pd.to_numeric(frame["analysis_case"], errors="coerce").fillna(0).eq(0)]
    rows: list[dict[str, object]] = []

    continuous = [
        ("Age at baseline", "age_at_baseline"),
        ("Prebaseline CK max", "prebaseline_ck_max"),
        ("Prebaseline creatinine last", "prebaseline_creatinine_last"),
        ("Prebaseline AST last", "prebaseline_ast_last"),
        ("Prebaseline ALT last", "prebaseline_alt_last"),
        ("Prebaseline sodium last", "prebaseline_sodium_last"),
        ("Prebaseline potassium last", "prebaseline_potassium_last"),
        ("Prebaseline calcium last", "prebaseline_calcium_last"),
        ("Prebaseline TSH last", "prebaseline_tsh_last"),
    ]
    for label, column in continuous:
        if column not in frame.columns:
            continue
        left = pd.to_numeric(cases[column], errors="coerce").dropna()
        right = pd.to_numeric(controls[column], errors="coerce").dropna()
        if left.empty and right.empty:
            continue
        diff, low, high = _mean_difference_ci(left, right)
        rows.append(
            {
                "analysis_split": split,
                "variable": label,
                "level": "",
                "cases": _format_mean_sd(cases, column),
                "controls": _format_mean_sd(controls, column),
                "effect": "mean_difference",
                "effect_size": _format_smd(diff),
                "effect_ci": _format_ci(low, high),
                "p_value": _format_p_value(_continuous_p_value(left, right)),
                "smd": _format_smd(_continuous_smd(left, right)),
            }
        )

    binary = [
        ("Prior kidney disease or AKI", "prior_kidney_disease_or_aki"),
        ("Remote sepsis history", "remote_sepsis_history"),
        ("Prior heat illness/dehydration", "prior_heat_illness_or_dehydration"),
        ("Prior exertion/exercise code", "prior_exertion_or_exercise_code"),
        ("Prior myopathy/muscle disease", "prior_myopathy_muscle_disease"),
        ("Sickle trait/disease coded", "prior_sickle_trait_or_disease"),
        ("Diabetes coded", "prior_diabetes"),
        ("Thyroid disease coded", "prior_thyroid_disease"),
        ("Liver disease coded", "prior_liver_disease"),
        ("Alcohol/substance condition coded", "prior_alcohol_substance_condition"),
        ("Prior injury/trauma history", "prior_injury_trauma_history"),
        ("CK tested prebaseline", "ck_tested_prebaseline"),
        ("CK tested during horizon", "ck_tested_during_horizon"),
    ]
    for label, column in binary:
        if column not in frame.columns:
            continue
        left = (pd.to_numeric(cases[column], errors="coerce").fillna(0) >= 1).astype(float)
        right = (pd.to_numeric(controls[column], errors="coerce").fillna(0) >= 1).astype(float)
        odds_ratio, low, high = _odds_ratio_ci(left, right)
        rows.append(
            {
                "analysis_split": split,
                "variable": label,
                "level": "yes",
                "cases": _format_binary(cases, column),
                "controls": _format_binary(controls, column),
                "effect": "odds_ratio",
                "effect_size": _format_smd(odds_ratio),
                "effect_ci": _format_ci(low, high),
                "p_value": _format_p_value(_binary_p_value(left, right)),
                "smd": _format_smd(_binary_smd(left, right)),
            }
        )

    categoricals = [
        ("Sex category", "model_sex_category", ["female", "male", "other_or_unknown", "missing"]),
        ("Genetic ancestry label", "model_ancestry_pred", ["eur", "afr", "amr", "eas", "sas", "mid", "other", "missing"]),
    ]
    for label, column, preferred in categoricals:
        if column not in frame.columns:
            continue
        left = cases[column].astype(str)
        right = controls[column].astype(str)
        categories = _ordered_categories([left, right], preferred)
        if not categories:
            continue
        p_value = _categorical_p_value(left, right, categories)
        smd = _categorical_smd(left, right, categories)
        for category in categories:
            rows.append(
                {
                    "analysis_split": split,
                    "variable": label,
                    "level": category,
                    "cases": _format_category_from_series(left, category),
                    "controls": _format_category_from_series(right, category),
                    "effect": "categorical",
                    "effect_size": "",
                    "effect_ci": "",
                    "p_value": _format_p_value(p_value),
                    "smd": _format_smd(smd),
                }
            )
    return rows


def _table1(model_input: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for split in ("train", "test"):
        rows.extend(_table1_rows_for_split(model_input, split))
    return pd.DataFrame(rows)


def _missingness(model_input: pd.DataFrame) -> pd.DataFrame:
    rows = []
    columns = [
        "age_at_baseline",
        "model_sex_category",
        "model_ancestry_pred",
        "prebaseline_ck_max",
        "prebaseline_creatinine_last",
        "prebaseline_ast_last",
        "prebaseline_alt_last",
        "prebaseline_sodium_last",
        "prebaseline_potassium_last",
        "prebaseline_calcium_last",
        "prebaseline_tsh_last",
        "ck_tested_prebaseline",
        "ck_tested_during_horizon",
    ]
    for column in columns:
        if column not in model_input.columns:
            continue
        mask = _missing_mask(model_input, column)
        rows.append(
            {
                "field": column,
                "n": len(model_input),
                "missing_n": int(mask.sum()),
                "missing_pct": float(100.0 * mask.mean()) if len(mask) else np.nan,
                "observed_n": int((~mask).sum()),
            }
        )
    return pd.DataFrame(rows)


def _write_characterization_report(
    *,
    consort: pd.DataFrame,
    split_summary: pd.DataFrame,
    table1: pd.DataFrame,
    missingness: pd.DataFrame,
    path: str,
) -> None:
    body = [
        "# EIR-Enriched Clinical Cohort Characterization",
        "",
        "Phenotype label: CK-confirmed, non-traumatic, non-septic rhabdomyolysis enriched for exertional rhabdomyolysis.",
        "",
        "## CONSORT Counts",
        "",
        _md_table(_suppress_small_counts(consort)),
        "",
        "## Train/Test Split",
        "",
        _md_table(_suppress_small_counts(split_summary)),
        "",
        "## Table 1",
        "",
        _md_table(table1),
        "",
        "## Missingness",
        "",
        _md_table(_suppress_small_counts(missingness)),
        "",
    ]
    write_text("\n".join(body), path)


def characterize_eir_cohort(config: ProjectConfig, cohort: pd.DataFrame, paths: ProjectPaths) -> dict[str, pd.DataFrame]:
    model_input = _analysis_model_input(cohort)
    consort = _eir_consort_counts(cohort)
    split_summary = _split_summary(model_input)
    table1 = _table1(model_input)
    missingness = _missingness(model_input)

    write_dataframe(consort, eir_consort_counts_path(paths))
    _write_md_table(consort, eir_consort_counts_report_path(paths), "EIR-Enriched CONSORT Counts")
    write_dataframe(model_input, eir_model_input_path(paths))
    write_dataframe(split_summary, eir_split_summary_path(paths))
    write_dataframe(table1, eir_table1_path(paths))
    write_dataframe(missingness, eir_missingness_path(paths))
    _write_characterization_report(
        consort=consort,
        split_summary=split_summary,
        table1=table1,
        missingness=missingness,
        path=eir_characterization_report_path(paths),
    )
    return {
        "consort": consort,
        "model_input": model_input,
        "split_summary": split_summary,
        "table1": table1,
        "missingness": missingness,
    }


def build_eir_cohort_artifacts(config: ProjectConfig) -> tuple[ProjectConfig, ProjectPaths, pd.DataFrame]:
    effective = _copy_with_eir_outcome(apply_runtime_defaults(config))
    paths = build_output_paths(effective)
    cohort = build_eir_cohort(effective)
    write_dataframe(cohort, paths.built_cohort_tsv)
    write_json(
        {
            "workflow": "eir_clinical_v1",
            "lookback_days": LOOKBACK_DAYS,
            "time_at_risk_days": TIME_AT_RISK_DAYS,
            "primary_label": "CK-confirmed non-traumatic/non-septic rhabdomyolysis enriched for exertional rhabdomyolysis",
            "rows": len(cohort),
            "primary_cases": int(pd.to_numeric(cohort["eir_primary_case"], errors="coerce").fillna(0).sum()),
            "eligible_controls": int(pd.to_numeric(cohort["eligible_control"], errors="coerce").fillna(0).sum()),
        },
        paths.manifest_json,
    )
    return effective, paths, cohort


def characterize_eir_artifacts(config: ProjectConfig) -> tuple[ProjectConfig, ProjectPaths, dict[str, pd.DataFrame]]:
    effective = _copy_with_eir_outcome(apply_runtime_defaults(config))
    paths = build_output_paths(effective)
    try:
        cohort = read_table(paths.built_cohort_tsv)
    except FileNotFoundError:
        _, paths, cohort = build_eir_cohort_artifacts(effective)
    for column in (
        "obs_start_date",
        "obs_end_date",
        "baseline_date",
        "horizon_start_date",
        "horizon_end_date",
        "first_rhabdo_date",
        "ck_confirming_date",
    ):
        if column in cohort.columns:
            cohort[column] = parse_date(cohort[column])
    outputs = characterize_eir_cohort(effective, cohort, paths)
    return effective, paths, outputs


def _eir_feature_specs(train_df: pd.DataFrame) -> tuple[FeatureSpec, ...]:
    specs: list[FeatureSpec] = []
    continuous = [
        ("age_at_baseline", "age_at_baseline_per_sd"),
        ("prebaseline_ck_max", "prebaseline_ck_max_per_sd"),
        ("prebaseline_creatinine_last", "prebaseline_creatinine_last_per_sd"),
        ("prebaseline_ast_last", "prebaseline_ast_last_per_sd"),
        ("prebaseline_alt_last", "prebaseline_alt_last_per_sd"),
        ("prebaseline_sodium_last", "prebaseline_sodium_last_per_sd"),
        ("prebaseline_potassium_last", "prebaseline_potassium_last_per_sd"),
        ("prebaseline_calcium_last", "prebaseline_calcium_last_per_sd"),
        ("prebaseline_tsh_last", "prebaseline_tsh_last_per_sd"),
    ]
    for source_column, name in continuous:
        if source_column not in train_df.columns:
            continue
        values = pd.to_numeric(train_df[source_column], errors="coerce")
        if values.notna().sum() < 5:
            continue
        median = float(values.median())
        filled = values.fillna(median)
        sd = float(filled.std(ddof=0))
        if not np.isfinite(sd) or sd == 0:
            continue
        specs.append(
            FeatureSpec(
                name=name,
                source_column=source_column,
                kind="continuous",
                median=median,
                mean=float(filled.mean()),
                sd=sd,
            )
        )

    binary = [
        "prior_kidney_disease_or_aki",
        "remote_sepsis_history",
        "prior_heat_illness_or_dehydration",
        "prior_exertion_or_exercise_code",
        "prior_myopathy_muscle_disease",
        "prior_sickle_trait_or_disease",
        "prior_diabetes",
        "prior_thyroid_disease",
        "prior_liver_disease",
        "prior_alcohol_substance_condition",
        "prior_injury_trauma_history",
        "ck_tested_prebaseline",
    ]
    for column in binary:
        if column not in train_df.columns:
            continue
        values = (pd.to_numeric(train_df[column], errors="coerce").fillna(0) >= 1).astype(float)
        if values.nunique(dropna=False) <= 1:
            continue
        specs.append(FeatureSpec(name=column, source_column=column, kind="binary"))

    categorical = {
        "model_sex_category": ("female", "male", "other_or_unknown", "missing"),
        "model_ancestry_pred": ("eur", "afr", "amr", "eas", "sas", "mid", "other", "missing"),
    }
    for source_column, preferred in categorical.items():
        if source_column not in train_df.columns:
            continue
        values = train_df[source_column].astype("string").str.lower().str.strip().fillna("missing").astype(str)
        observed = [category for category in preferred if category in set(values)]
        observed.extend(sorted(set(values).difference(observed)))
        if len(observed) <= 1:
            continue
        reference = observed[0]
        for category in observed[1:]:
            specs.append(
                FeatureSpec(
                    name=f"{source_column}={category}",
                    source_column=source_column,
                    kind="categorical",
                    reference=reference,
                )
            )
    return tuple(specs)


def _fit_eir_model(train_df: pd.DataFrame, *, l2_penalty: float) -> tuple[ClinicalModel, np.ndarray]:
    specs = _eir_feature_specs(train_df)
    if not specs:
        raise ValueError("No usable EIR clinical model features were found.")
    x = _feature_matrix(train_df, specs)
    y = pd.to_numeric(train_df[OUTCOME_COLUMN], errors="coerce").astype(int).to_numpy()
    beta = _fit_logistic(x, y, l2_penalty=l2_penalty)
    train_prob = _sigmoid(x @ beta)
    threshold = _choose_threshold(y, train_prob)
    return ClinicalModel(beta=beta, features=specs, threshold=threshold, l2_penalty=l2_penalty), train_prob


def _eir_cv_metrics(train_df: pd.DataFrame, *, l2_penalty: float) -> pd.DataFrame:
    rows = []
    for fold in sorted(pd.to_numeric(train_df["cv_fold"], errors="coerce").dropna().astype(int).unique().tolist()):
        fit_df = train_df[pd.to_numeric(train_df["cv_fold"], errors="coerce") != fold].copy()
        valid_df = train_df[pd.to_numeric(train_df["cv_fold"], errors="coerce") == fold].copy()
        if fit_df.empty or valid_df.empty or fit_df[OUTCOME_COLUMN].nunique() < 2 or valid_df[OUTCOME_COLUMN].nunique() < 2:
            continue
        model, _ = _fit_eir_model(fit_df, l2_penalty=l2_penalty)
        valid_prob = _predict(model, valid_df)
        y_valid = pd.to_numeric(valid_df[OUTCOME_COLUMN], errors="coerce").astype(int).to_numpy()
        row = _metrics(y_valid, valid_prob, threshold=model.threshold, label=f"cv_fold_{fold}")
        row["fold"] = fold
        rows.append(row)
    return pd.DataFrame(rows)


def _risk_deciles(predictions: pd.DataFrame) -> pd.DataFrame:
    test = predictions[predictions["analysis_split"] == "test"].copy()
    if test.empty:
        return pd.DataFrame(columns=["decile", "n", "mean_predicted_probability", "observed_case_fraction"])
    test = test.sort_values("predicted_probability").copy()
    groups = np.array_split(test.index.to_numpy(), min(10, len(test)))
    rows = []
    for index, group in enumerate(groups, start=1):
        frame = test.loc[group]
        rows.append(
            {
                "decile": index,
                "n": len(frame),
                "mean_predicted_probability": float(pd.to_numeric(frame["predicted_probability"], errors="coerce").mean()),
                "observed_case_fraction": float(pd.to_numeric(frame[OUTCOME_COLUMN], errors="coerce").mean()),
            }
        )
    return pd.DataFrame(rows)


def _write_risk_decile_svg(deciles: pd.DataFrame, path: str) -> None:
    width, height = 640, 420
    left, right, top, bottom = 70, 30, 48, 60
    inner_w = width - left - right
    inner_h = height - top - bottom
    bars = []
    if not deciles.empty:
        max_y = max(float(deciles["observed_case_fraction"].max()), float(deciles["mean_predicted_probability"].max()), 0.01)
        bar_w = inner_w / max(len(deciles), 1) * 0.35
        for idx, row in deciles.reset_index(drop=True).iterrows():
            x0 = left + idx * inner_w / len(deciles) + inner_w / len(deciles) * 0.2
            observed = float(row["observed_case_fraction"]) / max_y
            predicted = float(row["mean_predicted_probability"]) / max_y
            bars.append(
                f'<rect x="{x0:.1f}" y="{top + inner_h - observed * inner_h:.1f}" width="{bar_w:.1f}" height="{observed * inner_h:.1f}" fill="#0f766e" />'
            )
            bars.append(
                f'<rect x="{x0 + bar_w + 2:.1f}" y="{top + inner_h - predicted * inner_h:.1f}" width="{bar_w:.1f}" height="{predicted * inner_h:.1f}" fill="#475569" />'
            )
            bars.append(
                f'<text x="{x0 + bar_w:.1f}" y="{height - 38}" text-anchor="middle" font-size="10" font-family="Menlo, monospace">{int(row["decile"])}</text>'
            )
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
        '<rect width="100%" height="100%" fill="#fffdf7" />'
        f'<text x="{left}" y="28" font-size="20" font-family="Menlo, monospace">Test Risk Deciles</text>'
        f'<line x1="{left}" y1="{top + inner_h}" x2="{left + inner_w}" y2="{top + inner_h}" stroke="#1f2937" />'
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + inner_h}" stroke="#1f2937" />'
        + "".join(bars)
        + f'<text x="{left + inner_w - 160}" y="28" font-size="11" font-family="Menlo, monospace" fill="#0f766e">observed</text>'
        + f'<text x="{left + inner_w - 82}" y="28" font-size="11" font-family="Menlo, monospace" fill="#475569">predicted</text>'
        + f'<text x="{width / 2:.1f}" y="{height - 15}" text-anchor="middle" font-size="12" font-family="Menlo, monospace">Predicted-risk decile</text>'
        + "</svg>"
    )
    write_text(svg, path)


def _write_eir_model_report(
    *,
    metrics: pd.DataFrame,
    cv_metrics: pd.DataFrame,
    coefficients: pd.DataFrame,
    sparse_status: dict[str, Any],
    path: str,
) -> None:
    body = [
        "# EIR-Enriched Clinical-Only Model",
        "",
        "Model target: primary EIR-enriched incident rhabdomyolysis within 2 years after baseline.",
        "",
        "Predictors are restricted to pre-baseline demographic, comorbidity, injury, and lab features. Observation depth, total condition-record count, post-baseline data, peri-index sepsis, and peri-index renal injury are not predictors.",
        "",
        "## Metrics",
        "",
        _md_table(metrics),
        "",
        "## Cross-Validation Metrics",
        "",
        _md_table(cv_metrics) if not cv_metrics.empty else "_No cross-validation folds with both classes were available._",
        "",
        "## Coefficients",
        "",
        _md_table(coefficients.head(40)),
        "",
        "## Sparse OMOP Model",
        "",
        sparse_status.get("message", "Sparse OMOP model was not run."),
        "",
    ]
    write_text("\n".join(body), path)


def _write_sparse_status(paths: ProjectPaths, *, run_sparse: bool) -> dict[str, Any]:
    status = {
        "status": "skipped",
        "run_sparse_requested": bool(run_sparse),
        "message": (
            "Sparse OMOP model was requested but is deferred in this v1 implementation; "
            "the curated leakage-controlled clinical model was run."
            if run_sparse
            else "Sparse OMOP model not requested; curated leakage-controlled clinical model was run."
        ),
    }
    write_json(status, eir_sparse_status_path(paths))
    write_text("# Sparse OMOP Model\n\n" + status["message"] + "\n", eir_sparse_report_path(paths))
    return status


def run_eir_clinical_model(
    config: ProjectConfig,
    paths: ProjectPaths,
    *,
    l2_penalty: float = 1.0,
    run_sparse: bool = False,
) -> dict[str, pd.DataFrame | str]:
    model_input = read_table(eir_model_input_path(paths))
    for column in ("baseline_date", "horizon_start_date", "horizon_end_date", "first_rhabdo_date"):
        if column in model_input.columns:
            model_input[column] = parse_date(model_input[column])
    if OUTCOME_COLUMN not in model_input.columns:
        model_input[OUTCOME_COLUMN] = pd.to_numeric(model_input.get("eir_primary_case"), errors="coerce").fillna(0).astype(int)
    analysis = model_input[pd.to_numeric(model_input["case_control_eligible"], errors="coerce").fillna(0).eq(1)].copy()
    analysis[OUTCOME_COLUMN] = pd.to_numeric(analysis[OUTCOME_COLUMN], errors="coerce")
    analysis = analysis[analysis[OUTCOME_COLUMN].isin([0, 1])].copy()
    train = analysis[analysis["analysis_split"] == "train"].copy()
    test = analysis[analysis["analysis_split"] == "test"].copy()
    if train.empty or test.empty:
        raise ValueError("EIR clinical model requires non-empty train and test splits.")
    if train[OUTCOME_COLUMN].nunique() < 2 or test[OUTCOME_COLUMN].nunique() < 2:
        raise ValueError("EIR clinical model requires both cases and controls in train and test splits.")

    model, train_prob = _fit_eir_model(train, l2_penalty=l2_penalty)
    test_prob = _predict(model, test)
    y_train = train[OUTCOME_COLUMN].astype(int).to_numpy()
    y_test = test[OUTCOME_COLUMN].astype(int).to_numpy()
    metrics = pd.DataFrame(
        [
            _metrics(y_train, train_prob, threshold=model.threshold, label="train"),
            _metrics(y_test, test_prob, threshold=model.threshold, label="test"),
        ]
    )
    cv_metrics = _eir_cv_metrics(train, l2_penalty=l2_penalty)
    coefficients = _coefficient_table(model)
    calibration = _calibration_table(y_test, test_prob)
    predictions = pd.concat(
        [
            train.assign(predicted_probability=train_prob, predicted_label=(train_prob >= model.threshold).astype(int)),
            test.assign(predicted_probability=test_prob, predicted_label=(test_prob >= model.threshold).astype(int)),
        ],
        ignore_index=True,
    )
    keep_columns = [
        column
        for column in (
            "person_id",
            OUTCOME_COLUMN,
            "analysis_case",
            "analysis_split",
            "cv_fold",
            "baseline_date",
            "predicted_probability",
            "predicted_label",
        )
        if column in predictions.columns
    ]
    predictions = predictions[keep_columns].copy()
    deciles = _risk_deciles(predictions)
    sparse_status = _write_sparse_status(paths, run_sparse=run_sparse)

    write_dataframe(metrics, clinical_model_metrics_path(paths))
    write_dataframe(cv_metrics, clinical_model_cv_metrics_path(paths))
    write_dataframe(coefficients, clinical_model_coefficients_path(paths))
    write_dataframe(predictions, clinical_model_predictions_path(paths))
    write_dataframe(calibration, clinical_model_calibration_path(paths))
    write_dataframe(deciles, join_path(paths.run_root, "clinical", "model", "risk_deciles.tsv"))
    write_json(
        {
            "workflow": "eir_clinical_v1",
            "outcome_column": OUTCOME_COLUMN,
            "time_at_risk_days": TIME_AT_RISK_DAYS,
            "l2_penalty": l2_penalty,
            "threshold": model.threshold,
            "features": [spec.__dict__ for spec in model.features],
        },
        join_path(paths.run_root, "clinical", "model", "model_spec.json"),
    )
    _write_line_svg(_curve_points(y_test, test_prob, "roc"), clinical_model_roc_svg_path(paths), title="EIR Clinical Model ROC", x_label="False positive rate", y_label="True positive rate", diagonal=True)
    _write_line_svg(_curve_points(y_test, test_prob, "pr"), clinical_model_pr_svg_path(paths), title="EIR Clinical Model Precision-Recall", x_label="Recall", y_label="Precision")
    _write_calibration_svg(calibration, clinical_model_calibration_svg_path(paths))
    _write_risk_decile_svg(deciles, eir_risk_decile_path(paths))
    _write_eir_model_report(
        metrics=metrics,
        cv_metrics=cv_metrics,
        coefficients=coefficients,
        sparse_status=sparse_status,
        path=clinical_model_report_path(paths),
    )
    return {
        "metrics": metrics,
        "cv_metrics": cv_metrics,
        "coefficients": coefficients,
        "predictions": predictions,
        "calibration": calibration,
        "risk_deciles": deciles,
        "metrics_path": clinical_model_metrics_path(paths),
        "cv_metrics_path": clinical_model_cv_metrics_path(paths),
        "coefficients_path": clinical_model_coefficients_path(paths),
        "predictions_path": clinical_model_predictions_path(paths),
        "calibration_path": clinical_model_calibration_path(paths),
        "report_path": clinical_model_report_path(paths),
        "roc_svg_path": clinical_model_roc_svg_path(paths),
        "pr_svg_path": clinical_model_pr_svg_path(paths),
        "calibration_svg_path": clinical_model_calibration_svg_path(paths),
        "risk_decile_svg_path": eir_risk_decile_path(paths),
    }


__all__ = [
    "OUTCOME_COLUMN",
    "build_eir_cohort",
    "build_eir_cohort_artifacts",
    "characterize_eir_artifacts",
    "characterize_eir_cohort",
    "eir_characterization_report_path",
    "eir_consort_counts_path",
    "eir_consort_counts_report_path",
    "eir_missingness_path",
    "eir_model_input_path",
    "eir_risk_decile_path",
    "eir_sparse_report_path",
    "eir_sparse_status_path",
    "eir_split_summary_path",
    "eir_table1_path",
    "run_eir_clinical_model",
]
