"""Aggregate feasibility audit for incident rhabdomyolysis prediction."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

import pandas as pd

from .config import ProjectConfig
from .eir import (
    BIGQUERY_ON_DEMAND_USD_PER_TIB,
    BYTES_PER_TIB,
    LOOKBACK_DAYS,
    SEPSIS_COFACTOR,
    SEPSIS_EXCLUSION_END_DAYS,
    SEPSIS_EXCLUSION_START_DAYS,
    SUPPORTIVE_EIR_COFATORS,
    SUPPORTIVE_EIR_END_DAYS,
    SUPPORTIVE_EIR_START_DAYS,
    TIME_AT_RISK_DAYS,
    TRAUMA_COFACTORS,
    TRAUMA_EXCLUSION_END_DAYS,
    TRAUMA_EXCLUSION_START_DAYS,
    _bytes_from_tib,
    _cofactor_condition_predicate,
    _estimate_bigquery_cost_usd,
    _lab_measurement_predicate,
    build_eir_cohort,
)
from .io_utils import dry_run_bigquery_query, query_bigquery_dataframe, write_dataframe, write_json, write_text
from .paths import ProjectPaths, build_output_paths, join_path
from .phenotype_sql import _match_predicate, _qualify_table
from .preflight import apply_runtime_defaults


INCIDENT_ANALYSIS_NAME = "incident_rhabdo_prediction_v1"


def incident_feasibility_output_dir(paths: ProjectPaths) -> str:
    return join_path(paths.run_root, "feasibility")


def incident_feasibility_counts_path(paths: ProjectPaths) -> str:
    return join_path(incident_feasibility_output_dir(paths), "feasibility_counts.tsv")


def incident_case_funnel_path(paths: ProjectPaths) -> str:
    return join_path(incident_feasibility_output_dir(paths), "case_funnel.tsv")


def incident_control_funnel_path(paths: ProjectPaths) -> str:
    return join_path(incident_feasibility_output_dir(paths), "control_funnel.tsv")


def incident_baseline_history_bins_path(paths: ProjectPaths) -> str:
    return join_path(incident_feasibility_output_dir(paths), "baseline_history_bins.tsv")


def incident_microarray_overlap_counts_path(paths: ProjectPaths) -> str:
    return join_path(incident_feasibility_output_dir(paths), "microarray_overlap_counts.tsv")


def incident_feasibility_report_path(paths: ProjectPaths) -> str:
    return join_path(incident_feasibility_output_dir(paths), "report.md")


def _copy_with_incident_output(config: ProjectConfig) -> ProjectConfig:
    analysis = replace(
        config.analysis,
        analysis_name=INCIDENT_ANALYSIS_NAME,
        matched_outcome_column="incident_nontrauma_rhabdo_case",
    )
    return replace(config, analysis=analysis)


def _feasibility_base_sql(config: ProjectConfig) -> str:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    person_table = _qualify_table(cdr, config.phenotype.tables.person_table)
    observation_table = _qualify_table(cdr, config.phenotype.tables.observation_table)
    condition_table = _qualify_table(cdr, config.phenotype.tables.condition_table)
    measurement_table = _qualify_table(cdr, config.phenotype.tables.measurement_table)
    concept_table = _qualify_table(cdr, config.phenotype.tables.concept_table)

    rhabdo_predicate = _match_predicate(
        concept_id_column="condition_concept_id",
        concept_name_column="condition_concept_name",
        concept_ids=config.phenotype.broad.condition_concept_ids,
        concept_terms=config.phenotype.broad.condition_terms,
    )
    ck_predicate = _lab_measurement_predicate(config, "ck")
    ck_threshold = float(config.phenotype.definite.measurement_min or 5000)
    trauma_predicate = _cofactor_condition_predicate(config, TRAUMA_COFACTORS)
    sepsis_predicate = _cofactor_condition_predicate(config, (SEPSIS_COFACTOR,))
    supportive_predicate = _cofactor_condition_predicate(config, SUPPORTIVE_EIR_COFATORS)

    return f"""
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
    co.{config.phenotype.condition_concept_column} AS condition_concept_id,
    LOWER(COALESCE(condition_concept.concept_name, '')) AS condition_concept_name,
    UPPER(COALESCE(source_concept.vocabulary_id, '')) AS source_vocabulary_id
  FROM `{condition_table}` co
  LEFT JOIN `{concept_table}` condition_concept
    ON co.{config.phenotype.condition_concept_column} = condition_concept.concept_id
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
measurement_records AS (
  SELECT
    CAST(m.{config.phenotype.person_id_column} AS STRING) AS person_id,
    DATE(m.{config.phenotype.measurement_date_column}) AS measurement_date,
    m.{config.phenotype.measurement_concept_column} AS measurement_concept_id,
    LOWER(COALESCE(measurement_concept.concept_name, '')) AS measurement_concept_name,
    SAFE_CAST(m.{config.phenotype.measurement_value_column} AS FLOAT64) AS measurement_value
  FROM `{measurement_table}` m
  LEFT JOIN `{concept_table}` measurement_concept
    ON m.{config.phenotype.measurement_concept_column} = measurement_concept.concept_id
  WHERE m.{config.phenotype.measurement_date_column} IS NOT NULL
),
measurement_rollup AS (
  SELECT
    person_id,
    MAX(measurement_date) AS last_measurement_date
  FROM measurement_records
  GROUP BY person_id
),
activity_records AS (
  SELECT person_id, condition_date AS activity_date
  FROM filtered_condition_records
  UNION DISTINCT
  SELECT person_id, measurement_date AS activity_date
  FROM measurement_records
),
baseline AS (
  SELECT
    CAST(p.person_id AS STRING) AS person_id,
    observation_windows.obs_start_date,
    observation_windows.obs_end_date,
    DATE_ADD(observation_windows.obs_start_date, INTERVAL {LOOKBACK_DAYS} DAY) AS washout_end_date,
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
  LEFT JOIN condition_rollup
    ON CAST(p.person_id AS STRING) = condition_rollup.person_id
  LEFT JOIN measurement_rollup
    ON CAST(p.person_id AS STRING) = measurement_rollup.person_id
  CROSS JOIN denominator_mode
),
person_time AS (
  SELECT
    baseline.*,
    DATE_ADD(baseline.baseline_date, INTERVAL 1 DAY) AS horizon_start_date,
    DATE_ADD(baseline.baseline_date, INTERVAL {TIME_AT_RISK_DAYS} DAY) AS horizon_end_date,
    IF(
      baseline.baseline_date IS NOT NULL
      AND baseline.obs_start_date IS NOT NULL
      AND baseline.obs_end_date IS NOT NULL
      AND baseline.omop_condition_record_dates >= 2
      AND baseline.washout_end_date <= baseline.baseline_date,
      1,
      0
    ) AS eligible_ehr_denominator
  FROM baseline
),
rhabdo_rollup AS (
  SELECT
    pt.person_id,
    IF(COUNTIF({rhabdo_predicate}) > 0, 1, 0) AS has_rhabdo_ever,
    MIN(IF({rhabdo_predicate}, condition_date, NULL)) AS first_rhabdo_ever_date,
    MIN(IF(({rhabdo_predicate}) AND condition_date < pt.baseline_date, condition_date, NULL)) AS prior_rhabdo_date,
    MIN(IF(({rhabdo_predicate}) AND condition_date > pt.baseline_date AND condition_date <= pt.horizon_end_date, condition_date, NULL)) AS first_rhabdo_date,
    IF(COUNTIF(({rhabdo_predicate}) AND condition_date > pt.baseline_date AND condition_date <= pt.horizon_end_date) > 0, 1, 0) AS rhabdo_during_horizon
  FROM person_time pt
  LEFT JOIN condition_records
    ON pt.person_id = condition_records.person_id
  GROUP BY pt.person_id
),
person_index AS (
  SELECT
    pt.*,
    COALESCE(rhabdo_rollup.has_rhabdo_ever, 0) AS has_rhabdo_ever,
    rhabdo_rollup.first_rhabdo_ever_date,
    rhabdo_rollup.prior_rhabdo_date,
    rhabdo_rollup.first_rhabdo_date,
    COALESCE(rhabdo_rollup.rhabdo_during_horizon, 0) AS rhabdo_during_horizon
  FROM person_time pt
  LEFT JOIN rhabdo_rollup USING (person_id)
),
ck_rollup AS (
  SELECT
    pt.person_id,
    IF(COUNTIF(({ck_predicate}) AND measurement_value >= {ck_threshold} AND measurement_date < pt.baseline_date) > 0, 1, 0) AS prior_high_ck,
    IF(COUNTIF(({ck_predicate}) AND measurement_value >= {ck_threshold} AND measurement_date >= pt.horizon_start_date AND measurement_date <= pt.horizon_end_date) > 0, 1, 0) AS high_ck_during_horizon,
    IF(COUNTIF(({ck_predicate}) AND measurement_date >= pt.horizon_start_date AND measurement_date <= pt.horizon_end_date) > 0, 1, 0) AS ck_tested_during_horizon,
    IF(COUNTIF(({ck_predicate}) AND measurement_date < pt.baseline_date) > 0, 1, 0) AS ck_tested_prebaseline,
    IF(COUNTIF(
      ({ck_predicate})
      AND measurement_value >= {ck_threshold}
      AND pt.first_rhabdo_date IS NOT NULL
      AND DATE_DIFF(measurement_date, pt.first_rhabdo_date, DAY)
        BETWEEN {config.phenotype.definite.measurement_window_start_days} AND {config.phenotype.definite.measurement_window_end_days}
    ) > 0, 1, 0) AS ck_confirmed_index
  FROM person_index pt
  LEFT JOIN measurement_records
    ON pt.person_id = measurement_records.person_id
  GROUP BY pt.person_id
),
cofactor_rollup AS (
  SELECT
    pt.person_id,
    IF(COUNTIF(({trauma_predicate}) AND condition_date BETWEEN DATE_ADD(pt.first_rhabdo_date, INTERVAL {TRAUMA_EXCLUSION_START_DAYS} DAY) AND DATE_ADD(pt.first_rhabdo_date, INTERVAL {TRAUMA_EXCLUSION_END_DAYS} DAY)) > 0, 1, 0) AS excluded_periindex_trauma,
    IF(COUNTIF(({sepsis_predicate}) AND condition_date BETWEEN DATE_ADD(pt.first_rhabdo_date, INTERVAL {SEPSIS_EXCLUSION_START_DAYS} DAY) AND DATE_ADD(pt.first_rhabdo_date, INTERVAL {SEPSIS_EXCLUSION_END_DAYS} DAY)) > 0, 1, 0) AS periindex_sepsis_flag,
    IF(COUNTIF(({supportive_predicate}) AND condition_date BETWEEN DATE_ADD(pt.first_rhabdo_date, INTERVAL {SUPPORTIVE_EIR_START_DAYS} DAY) AND DATE_ADD(pt.first_rhabdo_date, INTERVAL {SUPPORTIVE_EIR_END_DAYS} DAY)) > 0, 1, 0) AS supportive_exertion_heat_dehydration_code
  FROM person_index pt
  LEFT JOIN condition_records
    ON pt.person_id = condition_records.person_id
  GROUP BY pt.person_id
),
assembled AS (
  SELECT
    pt.*,
    COALESCE(ck_rollup.prior_high_ck, 0) AS prior_high_ck,
    COALESCE(ck_rollup.high_ck_during_horizon, 0) AS high_ck_during_horizon,
    COALESCE(ck_rollup.ck_tested_during_horizon, 0) AS ck_tested_during_horizon,
    COALESCE(ck_rollup.ck_tested_prebaseline, 0) AS ck_tested_prebaseline,
    COALESCE(ck_rollup.ck_confirmed_index, 0) AS ck_confirmed_index,
    COALESCE(cofactor_rollup.excluded_periindex_trauma, 0) AS excluded_periindex_trauma,
    COALESCE(cofactor_rollup.periindex_sepsis_flag, 0) AS periindex_sepsis_flag,
    COALESCE(cofactor_rollup.supportive_exertion_heat_dehydration_code, 0) AS supportive_exertion_heat_dehydration_code
  FROM person_index pt
  LEFT JOIN ck_rollup USING (person_id)
  LEFT JOIN cofactor_rollup USING (person_id)
),
participant_flags AS (
  SELECT
    assembled.*,
    NULLIF(
      GREATEST(
        COALESCE(last_condition_date, DATE '1900-01-01'),
        COALESCE(last_measurement_date, DATE '1900-01-01')
      ),
      DATE '1900-01-01'
    ) AS last_clinical_activity_date,
    IF(eligible_ehr_denominator = 1 AND prior_rhabdo_date IS NULL AND prior_high_ck = 0, 1, 0) AS incident_denominator,
    IF(
      obs_end_date >= horizon_end_date
      AND (
        NULLIF(
          GREATEST(
            COALESCE(last_condition_date, DATE '1900-01-01'),
            COALESCE(last_measurement_date, DATE '1900-01-01')
          ),
          DATE '1900-01-01'
        ) IS NULL
        OR NULLIF(
          GREATEST(
            COALESCE(last_condition_date, DATE '1900-01-01'),
            COALESCE(last_measurement_date, DATE '1900-01-01')
          ),
          DATE '1900-01-01'
        ) >= baseline_date
      ),
      1,
      0
    ) AS observed_followup_through_horizon
  FROM assembled
),
feasibility AS (
  SELECT
    *,
    IF(incident_denominator = 1 AND first_rhabdo_date IS NOT NULL, 1, 0) AS incident_rhabdo_case,
    IF(incident_denominator = 1 AND first_rhabdo_date IS NOT NULL AND excluded_periindex_trauma = 0, 1, 0) AS incident_nontrauma_rhabdo_case,
    IF(incident_denominator = 1 AND first_rhabdo_date IS NOT NULL AND excluded_periindex_trauma = 0 AND ck_confirmed_index = 1, 1, 0) AS ck_confirmed_nontrauma_case,
    IF(incident_denominator = 1 AND first_rhabdo_date IS NOT NULL AND excluded_periindex_trauma = 0 AND supportive_exertion_heat_dehydration_code = 1, 1, 0) AS supportive_nontrauma_case,
    IF(incident_denominator = 1 AND first_rhabdo_date IS NOT NULL AND excluded_periindex_trauma = 0 AND periindex_sepsis_flag = 1, 1, 0) AS sepsis_flagged_nontrauma_case,
    IF(incident_denominator = 1 AND rhabdo_during_horizon = 0 AND high_ck_during_horizon = 1, 1, 0) AS ck_only_no_rhabdo,
    IF(
      incident_denominator = 1
      AND observed_followup_through_horizon = 1
      AND rhabdo_during_horizon = 0
      AND high_ck_during_horizon = 0,
      1,
      0
    ) AS eligible_control,
    IF(
      incident_denominator = 1
      AND observed_followup_through_horizon = 1
      AND rhabdo_during_horizon = 0
      AND high_ck_during_horizon = 0
      AND ck_tested_during_horizon = 1,
      1,
      0
    ) AS ck_tested_control
  FROM participant_flags
)
""".strip()


def render_incident_feasibility_sql(config: ProjectConfig) -> str:
    base = _feasibility_base_sql(config)
    return f"""
{base}
SELECT 'case_funnel' AS section, 'All participants' AS metric, COUNT(*) AS n
FROM feasibility
UNION ALL
SELECT 'case_funnel', 'Any rhabdo ever', COUNTIF(has_rhabdo_ever = 1)
FROM feasibility
UNION ALL
SELECT 'case_funnel', 'First rhabdo after baseline within 2 years', COUNTIF(first_rhabdo_date IS NOT NULL)
FROM feasibility
UNION ALL
SELECT 'case_funnel', 'Incident rhabdo with 365d lookback and >=2 condition dates', COUNTIF(first_rhabdo_date IS NOT NULL AND eligible_ehr_denominator = 1)
FROM feasibility
UNION ALL
SELECT 'case_funnel', 'Incident rhabdo with no prior rhabdo or CK >=5000', COUNTIF(incident_rhabdo_case = 1)
FROM feasibility
UNION ALL
SELECT 'case_funnel', 'Non-trauma incident rhabdo cases', COUNTIF(incident_nontrauma_rhabdo_case = 1)
FROM feasibility
UNION ALL
SELECT 'case_funnel', 'CK-confirmed non-trauma incident rhabdo cases', COUNTIF(ck_confirmed_nontrauma_case = 1)
FROM feasibility
UNION ALL
SELECT 'case_funnel', 'Exertion/heat/dehydration-coded non-trauma cases', COUNTIF(supportive_nontrauma_case = 1)
FROM feasibility
UNION ALL
SELECT 'case_funnel', 'Peri-index sepsis-flagged non-trauma cases', COUNTIF(sepsis_flagged_nontrauma_case = 1)
FROM feasibility
UNION ALL
SELECT 'control_funnel', 'Incident denominator', COUNTIF(incident_denominator = 1)
FROM feasibility
UNION ALL
SELECT 'control_funnel', 'No rhabdo during 2-year horizon', COUNTIF(incident_denominator = 1 AND rhabdo_during_horizon = 0)
FROM feasibility
UNION ALL
SELECT 'control_funnel', 'No rhabdo and no CK >=5000 during horizon', COUNTIF(incident_denominator = 1 AND rhabdo_during_horizon = 0 AND high_ck_during_horizon = 0)
FROM feasibility
UNION ALL
SELECT 'control_funnel', 'Observed through 2-year horizon', COUNTIF(incident_denominator = 1 AND observed_followup_through_horizon = 1)
FROM feasibility
UNION ALL
SELECT 'control_funnel', 'Eligible controls', COUNTIF(eligible_control = 1)
FROM feasibility
UNION ALL
SELECT 'control_funnel', 'CK-tested eligible controls', COUNTIF(ck_tested_control = 1)
FROM feasibility
UNION ALL
SELECT 'feasibility_counts', 'Primary non-trauma incident cases', COUNTIF(incident_nontrauma_rhabdo_case = 1)
FROM feasibility
UNION ALL
SELECT 'feasibility_counts', 'CK-confirmed non-trauma cases', COUNTIF(ck_confirmed_nontrauma_case = 1)
FROM feasibility
UNION ALL
SELECT 'feasibility_counts', 'Eligible controls', COUNTIF(eligible_control = 1)
FROM feasibility
UNION ALL
SELECT 'feasibility_counts', 'Case-control rows', COUNTIF(incident_nontrauma_rhabdo_case = 1 OR eligible_control = 1)
FROM feasibility
UNION ALL
SELECT
  'baseline_history_bins',
  CASE
    WHEN first_rhabdo_ever_date IS NULL OR obs_start_date IS NULL THEN 'No rhabdo or no observation start'
    WHEN DATE_DIFF(first_rhabdo_ever_date, obs_start_date, DAY) < 365 THEN '<365d before first rhabdo'
    WHEN DATE_DIFF(first_rhabdo_ever_date, obs_start_date, DAY) < 730 THEN '365-729d before first rhabdo'
    WHEN DATE_DIFF(first_rhabdo_ever_date, obs_start_date, DAY) < 1095 THEN '730-1094d before first rhabdo'
    ELSE '>=1095d before first rhabdo'
  END AS metric,
  COUNTIF(has_rhabdo_ever = 1) AS n
FROM feasibility
WHERE has_rhabdo_ever = 1
GROUP BY metric
ORDER BY section, metric
""".strip()


def render_incident_feasibility_id_sql(config: ProjectConfig) -> str:
    base = _feasibility_base_sql(config)
    return f"""
{base}
SELECT
  person_id,
  incident_nontrauma_rhabdo_case,
  ck_confirmed_nontrauma_case,
  supportive_nontrauma_case,
  sepsis_flagged_nontrauma_case,
  eligible_control,
  ck_tested_control
FROM feasibility
WHERE incident_nontrauma_rhabdo_case = 1
  OR eligible_control = 1
""".strip()


def _section(frame: pd.DataFrame, name: str) -> pd.DataFrame:
    subset = frame[frame["section"].astype(str) == name].copy()
    return subset.drop(columns=["section"]).reset_index(drop=True)


def _md_table(frame: pd.DataFrame) -> str:
    try:
        return frame.to_markdown(index=False)
    except ImportError:
        return frame.to_string(index=False)


def _read_fam_ids(path: str) -> set[str]:
    fam = pd.read_csv(path, sep=r"\s+", header=None, usecols=[1], names=["IID"], dtype=str)
    return set(fam["IID"].dropna().astype(str))


def _numeric_column(frame: pd.DataFrame, column: str, default: int = 0) -> pd.Series:
    if column in frame.columns:
        values = frame[column]
    else:
        values = pd.Series(default, index=frame.index)
    return pd.to_numeric(values, errors="coerce").fillna(default)


def _microarray_overlap_counts(flags: pd.DataFrame, fam_path: str) -> pd.DataFrame:
    fam_ids = _read_fam_ids(fam_path)
    output = flags.copy()
    output["person_id"] = output["person_id"].astype(str)
    overlap = output[output["person_id"].isin(fam_ids)].copy()

    rows = [
        ("Microarray FAM IDs", len(fam_ids)),
        ("Feasibility case/control IDs queried", int(output["person_id"].nunique())),
        ("Primary non-trauma incident cases", int(pd.to_numeric(output["incident_nontrauma_rhabdo_case"], errors="coerce").fillna(0).sum())),
        ("Primary non-trauma incident cases with microarray", int(pd.to_numeric(overlap["incident_nontrauma_rhabdo_case"], errors="coerce").fillna(0).sum())),
        ("CK-confirmed non-trauma cases with microarray", int(pd.to_numeric(overlap["ck_confirmed_nontrauma_case"], errors="coerce").fillna(0).sum())),
        ("Exertion/heat/dehydration-coded cases with microarray", int(pd.to_numeric(overlap["supportive_nontrauma_case"], errors="coerce").fillna(0).sum())),
        ("Peri-index sepsis-flagged cases with microarray", int(pd.to_numeric(overlap["sepsis_flagged_nontrauma_case"], errors="coerce").fillna(0).sum())),
        ("Eligible controls with microarray", int(pd.to_numeric(overlap["eligible_control"], errors="coerce").fillna(0).sum())),
        ("CK-tested eligible controls with microarray", int(pd.to_numeric(overlap["ck_tested_control"], errors="coerce").fillna(0).sum())),
    ]
    return pd.DataFrame(rows, columns=["metric", "n"])


def _counts_from_local_cohort(cohort: pd.DataFrame) -> pd.DataFrame:
    frame = cohort.copy()
    for column in (
        "has_rhabdo_ever",
        "eligible_ehr_denominator",
        "incident_denominator",
        "rhabdo_during_horizon",
        "prior_high_ck",
        "high_ck_during_horizon",
        "excluded_periindex_trauma",
        "periindex_sepsis_flag",
        "supportive_exertion_heat_dehydration_code",
        "ck_confirmed_index",
        "observed_followup_through_horizon",
        "eligible_control",
        "ck_tested_control",
    ):
        if column not in frame.columns:
            frame[column] = 0
        frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0).astype(int)

    if "periindex_sepsis_flag" not in cohort.columns and "excluded_pre_or_same_day_sepsis" in cohort.columns:
        frame["periindex_sepsis_flag"] = pd.to_numeric(cohort["excluded_pre_or_same_day_sepsis"], errors="coerce").fillna(0).astype(int)
    if "ck_tested_control" not in cohort.columns and "eir_ck_tested_control" in cohort.columns:
        frame["ck_tested_control"] = pd.to_numeric(cohort["eir_ck_tested_control"], errors="coerce").fillna(0).astype(int)
    if "has_rhabdo_ever" not in cohort.columns:
        first_ever = frame.get("first_rhabdo_ever_date", frame.get("first_rhabdo_date"))
        frame["has_rhabdo_ever"] = pd.Series(first_ever).notna().astype(int)

    frame["first_rhabdo_date"] = pd.to_datetime(frame.get("first_rhabdo_date"), errors="coerce")
    frame["first_rhabdo_ever_date"] = pd.to_datetime(frame.get("first_rhabdo_ever_date", frame["first_rhabdo_date"]), errors="coerce")
    frame["obs_start_date"] = pd.to_datetime(frame.get("obs_start_date"), errors="coerce")
    frame["incident_rhabdo_case"] = ((frame["incident_denominator"] == 1) & frame["first_rhabdo_date"].notna()).astype(int)
    frame["incident_nontrauma_rhabdo_case"] = ((frame["incident_rhabdo_case"] == 1) & (frame["excluded_periindex_trauma"] == 0)).astype(int)
    frame["ck_confirmed_nontrauma_case"] = ((frame["incident_nontrauma_rhabdo_case"] == 1) & (frame["ck_confirmed_index"] == 1)).astype(int)
    frame["supportive_nontrauma_case"] = ((frame["incident_nontrauma_rhabdo_case"] == 1) & (frame["supportive_exertion_heat_dehydration_code"] == 1)).astype(int)
    frame["sepsis_flagged_nontrauma_case"] = ((frame["incident_nontrauma_rhabdo_case"] == 1) & (frame["periindex_sepsis_flag"] == 1)).astype(int)
    frame["eligible_control"] = pd.to_numeric(frame["eligible_control"], errors="coerce").fillna(0).astype(int)

    rows = [
        ("case_funnel", "All participants", len(frame)),
        ("case_funnel", "Any rhabdo ever", int(frame["has_rhabdo_ever"].sum())),
        ("case_funnel", "First rhabdo after baseline within 2 years", int(frame["first_rhabdo_date"].notna().sum())),
        ("case_funnel", "Incident rhabdo with 365d lookback and >=2 condition dates", int(((frame["first_rhabdo_date"].notna()) & (frame["eligible_ehr_denominator"] == 1)).sum())),
        ("case_funnel", "Incident rhabdo with no prior rhabdo or CK >=5000", int(frame["incident_rhabdo_case"].sum())),
        ("case_funnel", "Non-trauma incident rhabdo cases", int(frame["incident_nontrauma_rhabdo_case"].sum())),
        ("case_funnel", "CK-confirmed non-trauma incident rhabdo cases", int(frame["ck_confirmed_nontrauma_case"].sum())),
        ("case_funnel", "Exertion/heat/dehydration-coded non-trauma cases", int(frame["supportive_nontrauma_case"].sum())),
        ("case_funnel", "Peri-index sepsis-flagged non-trauma cases", int(frame["sepsis_flagged_nontrauma_case"].sum())),
        ("control_funnel", "Incident denominator", int(frame["incident_denominator"].sum())),
        ("control_funnel", "No rhabdo during 2-year horizon", int(((frame["incident_denominator"] == 1) & (frame["rhabdo_during_horizon"] == 0)).sum())),
        ("control_funnel", "No rhabdo and no CK >=5000 during horizon", int(((frame["incident_denominator"] == 1) & (frame["rhabdo_during_horizon"] == 0) & (frame["high_ck_during_horizon"] == 0)).sum())),
        ("control_funnel", "Observed through 2-year horizon", int(((frame["incident_denominator"] == 1) & (frame["observed_followup_through_horizon"] == 1)).sum())),
        ("control_funnel", "Eligible controls", int(frame["eligible_control"].sum())),
        ("control_funnel", "CK-tested eligible controls", int(frame["ck_tested_control"].sum())),
        ("feasibility_counts", "Primary non-trauma incident cases", int(frame["incident_nontrauma_rhabdo_case"].sum())),
        ("feasibility_counts", "CK-confirmed non-trauma cases", int(frame["ck_confirmed_nontrauma_case"].sum())),
        ("feasibility_counts", "Eligible controls", int(frame["eligible_control"].sum())),
        ("feasibility_counts", "Case-control rows", int(((frame["incident_nontrauma_rhabdo_case"] == 1) | (frame["eligible_control"] == 1)).sum())),
    ]

    has_rhabdo = frame[frame["has_rhabdo_ever"] == 1].copy()
    deltas = (has_rhabdo["first_rhabdo_ever_date"] - has_rhabdo["obs_start_date"]).dt.days
    bins = pd.cut(
        deltas,
        bins=[-10**9, 364, 729, 1094, 10**9],
        labels=[
            "<365d before first rhabdo",
            "365-729d before first rhabdo",
            "730-1094d before first rhabdo",
            ">=1095d before first rhabdo",
        ],
    )
    for label, count in bins.value_counts(sort=False).items():
        rows.append(("baseline_history_bins", str(label), int(count)))
    return pd.DataFrame(rows, columns=["section", "metric", "n"])


def _write_outputs(paths: ProjectPaths, combined: pd.DataFrame, microarray: pd.DataFrame | None = None) -> dict[str, pd.DataFrame]:
    outputs = {
        "feasibility_counts": _section(combined, "feasibility_counts"),
        "case_funnel": _section(combined, "case_funnel"),
        "control_funnel": _section(combined, "control_funnel"),
        "baseline_history_bins": _section(combined, "baseline_history_bins"),
    }
    write_dataframe(outputs["feasibility_counts"], incident_feasibility_counts_path(paths))
    write_dataframe(outputs["case_funnel"], incident_case_funnel_path(paths))
    write_dataframe(outputs["control_funnel"], incident_control_funnel_path(paths))
    write_dataframe(outputs["baseline_history_bins"], incident_baseline_history_bins_path(paths))
    if microarray is not None:
        outputs["microarray_overlap_counts"] = microarray
        write_dataframe(microarray, incident_microarray_overlap_counts_path(paths))
    lines = [
        "# Incident Rhabdomyolysis Feasibility Audit",
        "",
        f"- Lookback: {LOOKBACK_DAYS} days",
        f"- Time at risk: {TIME_AT_RISK_DAYS} days",
        "- Primary case for feasibility: incident non-traumatic rhabdomyolysis diagnosis.",
        "- Sepsis is flagged, not excluded.",
        "",
        "## Feasibility Counts",
        "",
        _md_table(outputs["feasibility_counts"]),
        "",
        "## Case Funnel",
        "",
        _md_table(outputs["case_funnel"]),
        "",
        "## Control Funnel",
        "",
        _md_table(outputs["control_funnel"]),
        "",
        "## Baseline History Bins",
        "",
        _md_table(outputs["baseline_history_bins"]),
    ]
    if microarray is not None:
        lines.extend(["", "## Microarray Overlap", "", _md_table(microarray)])
    write_text("\n".join(lines) + "\n", incident_feasibility_report_path(paths))
    return outputs


def estimate_incident_feasibility_artifacts(
    config: ProjectConfig,
    *,
    max_tib: float | None = None,
    write_sql_path: str | None = None,
) -> tuple[ProjectConfig, ProjectPaths, dict[str, Any]]:
    effective = apply_runtime_defaults(_copy_with_incident_output(config))
    paths = build_output_paths(effective)
    maximum_bytes_billed = _bytes_from_tib(max_tib)
    if effective.phenotype.tables.cohort_table:
        estimate = {
            "workflow": INCIDENT_ANALYSIS_NAME,
            "mode": "local",
            "total_bytes_processed": 0,
            "total_tib_processed": 0.0,
            "maximum_bytes_billed": maximum_bytes_billed,
            "would_exceed_maximum_bytes_billed": False,
            "estimated_query_cost_usd": 0.0,
            "message": "Local cohort tables are configured; BigQuery dry-run is not needed.",
        }
        return effective, paths, estimate
    sql = render_incident_feasibility_sql(effective)
    if write_sql_path:
        write_text(sql, write_sql_path)
    estimate = dry_run_bigquery_query(sql, maximum_bytes_billed=maximum_bytes_billed)
    total_bytes = int(estimate.get("total_bytes_processed", 0))
    estimate.update(
        {
            "workflow": INCIDENT_ANALYSIS_NAME,
            "mode": "bigquery",
            "sql_path": write_sql_path,
            "estimated_query_cost_usd": _estimate_bigquery_cost_usd(total_bytes),
            "on_demand_usd_per_tib": BIGQUERY_ON_DEMAND_USD_PER_TIB,
        }
    )
    return effective, paths, estimate


def run_incident_feasibility(
    config: ProjectConfig,
    *,
    max_tib: float | None = None,
    write_sql_path: str | None = None,
    microarray_fam: str | None = None,
) -> tuple[ProjectConfig, ProjectPaths, dict[str, pd.DataFrame]]:
    effective = apply_runtime_defaults(_copy_with_incident_output(config))
    paths = build_output_paths(effective)
    maximum_bytes_billed = _bytes_from_tib(max_tib)

    microarray_counts: pd.DataFrame | None = None
    if effective.phenotype.tables.cohort_table:
        cohort = build_eir_cohort(effective)
        combined = _counts_from_local_cohort(cohort)
        if microarray_fam:
            flags = cohort.copy()
            flags["incident_nontrauma_rhabdo_case"] = (
                _numeric_column(flags, "incident_denominator").eq(1)
                & pd.to_datetime(flags["first_rhabdo_date"], errors="coerce").notna()
                & _numeric_column(flags, "excluded_periindex_trauma").eq(0)
            ).astype(int)
            flags["ck_confirmed_nontrauma_case"] = (
                _numeric_column(flags, "incident_nontrauma_rhabdo_case").eq(1)
                & _numeric_column(flags, "ck_confirmed_index").eq(1)
            ).astype(int)
            flags["supportive_nontrauma_case"] = (
                _numeric_column(flags, "incident_nontrauma_rhabdo_case").eq(1)
                & _numeric_column(flags, "supportive_exertion_heat_dehydration_code").eq(1)
            ).astype(int)
            sepsis_column = "periindex_sepsis_flag" if "periindex_sepsis_flag" in flags.columns else "excluded_pre_or_same_day_sepsis"
            flags["sepsis_flagged_nontrauma_case"] = (
                _numeric_column(flags, "incident_nontrauma_rhabdo_case").eq(1)
                & _numeric_column(flags, sepsis_column).eq(1)
            ).astype(int)
            ck_tested_column = "eir_ck_tested_control" if "eir_ck_tested_control" in flags.columns else "ck_tested_control"
            flags["ck_tested_control"] = _numeric_column(flags, ck_tested_column).astype(int)
            microarray_counts = _microarray_overlap_counts(
                flags[
                    [
                        "person_id",
                        "incident_nontrauma_rhabdo_case",
                        "ck_confirmed_nontrauma_case",
                        "supportive_nontrauma_case",
                        "sepsis_flagged_nontrauma_case",
                        "eligible_control",
                        "ck_tested_control",
                    ]
                ],
                microarray_fam,
            )
        outputs = _write_outputs(paths, combined, microarray_counts)
        return effective, paths, outputs

    sql = render_incident_feasibility_sql(effective)
    if write_sql_path:
        write_text(sql, write_sql_path)
    combined = query_bigquery_dataframe(
        sql,
        maximum_bytes_billed=maximum_bytes_billed,
        progress_label="Incident rhabdo feasibility",
    )
    if microarray_fam:
        id_sql = render_incident_feasibility_id_sql(effective)
        flags = query_bigquery_dataframe(
            id_sql,
            maximum_bytes_billed=maximum_bytes_billed,
            progress_label="Incident rhabdo feasibility microarray overlap",
        )
        microarray_counts = _microarray_overlap_counts(flags, microarray_fam)
    outputs = _write_outputs(paths, combined, microarray_counts)
    metadata = {
        "workflow": INCIDENT_ANALYSIS_NAME,
        "lookback_days": LOOKBACK_DAYS,
        "time_at_risk_days": TIME_AT_RISK_DAYS,
        "max_tib": max_tib,
        "microarray_fam": microarray_fam,
    }
    write_json(metadata, join_path(incident_feasibility_output_dir(paths), "feasibility_qc.json"))
    return effective, paths, outputs


__all__ = [
    "estimate_incident_feasibility_artifacts",
    "incident_baseline_history_bins_path",
    "incident_case_funnel_path",
    "incident_control_funnel_path",
    "incident_feasibility_counts_path",
    "incident_feasibility_report_path",
    "incident_microarray_overlap_counts_path",
    "render_incident_feasibility_id_sql",
    "render_incident_feasibility_sql",
    "run_incident_feasibility",
]
