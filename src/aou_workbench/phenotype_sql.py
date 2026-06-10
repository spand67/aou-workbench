"""SQL rendering helpers for Workbench phenotype and covariate extraction."""

from __future__ import annotations

from .config import CaseTierRule, ClinicalCofactorRule, ProjectConfig


def _escape_term(term: str) -> str:
    return term.replace("\\", "\\\\").replace("'", "\\'")


def _qualify_table(cdr: str, reference: str) -> str:
    normalized = reference[5:] if reference.startswith("bq://") else reference
    if normalized.count(".") == 2 and "/" not in normalized:
        return normalized
    return f"{cdr}.{normalized}"


def _match_predicate(
    *,
    concept_id_column: str,
    concept_name_column: str,
    concept_ids: tuple[int, ...],
    concept_terms: tuple[str, ...],
) -> str:
    predicates: list[str] = []
    if concept_ids:
        predicates.append(f"{concept_id_column} IN ({', '.join(str(value) for value in concept_ids)})")
    for term in concept_terms:
        predicates.append(f"LOWER({concept_name_column}) LIKE '%{_escape_term(term.lower())}%'")
    return " OR ".join(predicates) if predicates else "FALSE"


def _empty_result_sql(*, columns: tuple[str, ...]) -> str:
    projection = ", ".join(columns)
    return f"SELECT {projection} FROM (SELECT 1 AS _unused) WHERE FALSE"


def render_case_tier_sql(config: ProjectConfig, tier: CaseTierRule) -> str:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    condition_table = _qualify_table(cdr, config.phenotype.tables.condition_table)
    measurement_table = _qualify_table(cdr, config.phenotype.tables.measurement_table)
    concept_table = _qualify_table(cdr, config.phenotype.tables.concept_table)
    condition_predicate = _match_predicate(
        concept_id_column=f"co.{config.phenotype.condition_concept_column}",
        concept_name_column="condition_concept.concept_name",
        concept_ids=tier.condition_concept_ids,
        concept_terms=tier.condition_terms,
    )
    measurement_predicate = _match_predicate(
        concept_id_column=f"m.{config.phenotype.measurement_concept_column}",
        concept_name_column="measurement_concept.concept_name",
        concept_ids=tier.measurement_concept_ids,
        concept_terms=tier.measurement_terms,
    )
    condition_sql = (
        f"""
SELECT
  CAST(co.{config.phenotype.person_id_column} AS STRING) AS person_id,
  MIN(DATE(co.{config.phenotype.condition_date_column})) AS condition_date
FROM `{condition_table}` co
LEFT JOIN `{concept_table}` condition_concept
  ON co.{config.phenotype.condition_concept_column} = condition_concept.concept_id
WHERE {condition_predicate}
GROUP BY person_id
""".strip()
        if condition_predicate != "FALSE"
        else _empty_result_sql(
            columns=("CAST(NULL AS STRING) AS person_id", "CAST(NULL AS DATE) AS condition_date")
        )
    )
    measurement_threshold = (
        f"AND SAFE_CAST(m.{config.phenotype.measurement_value_column} AS FLOAT64) >= {tier.measurement_min}"
        if tier.measurement_min is not None
        else ""
    )
    measurement_events_sql = (
        f"""
SELECT
  CAST(m.{config.phenotype.person_id_column} AS STRING) AS person_id,
  DATE(m.{config.phenotype.measurement_date_column}) AS measurement_date,
  SAFE_CAST(m.{config.phenotype.measurement_value_column} AS FLOAT64) AS measurement_value
FROM `{measurement_table}` m
LEFT JOIN `{concept_table}` measurement_concept
  ON m.{config.phenotype.measurement_concept_column} = measurement_concept.concept_id
WHERE {measurement_predicate}
  {measurement_threshold}
""".strip()
        if measurement_predicate != "FALSE"
        else _empty_result_sql(
            columns=(
                "CAST(NULL AS STRING) AS person_id",
                "CAST(NULL AS DATE) AS measurement_date",
                "CAST(NULL AS FLOAT64) AS measurement_value",
            )
        )
    )
    if tier.require_condition and tier.require_measurement:
        return f"""
WITH condition_events AS (
  {condition_sql}
),
measurement_events AS (
  {measurement_events_sql}
)
SELECT
  condition_events.person_id,
  condition_events.condition_date,
  MIN(measurement_events.measurement_date) AS measurement_date,
  MAX(measurement_events.measurement_value) AS measurement_value
FROM condition_events
JOIN measurement_events USING (person_id)
WHERE DATE_DIFF(measurement_events.measurement_date, condition_events.condition_date, DAY)
  BETWEEN {tier.measurement_window_start_days} AND {tier.measurement_window_end_days}
GROUP BY person_id, condition_date
""".strip()
    if tier.require_condition and not tier.require_measurement:
        return f"""
WITH condition_events AS (
  {condition_sql}
)
SELECT
  condition_events.person_id,
  condition_events.condition_date,
  CAST(NULL AS DATE) AS measurement_date,
  CAST(NULL AS FLOAT64) AS measurement_value
FROM condition_events
""".strip()
    if tier.require_measurement and not tier.require_condition:
        return f"""
WITH measurement_events AS (
  {measurement_events_sql}
)
SELECT
  measurement_events.person_id,
  CAST(NULL AS DATE) AS condition_date,
  MIN(measurement_events.measurement_date) AS measurement_date,
  MAX(measurement_events.measurement_value) AS measurement_value
FROM measurement_events
GROUP BY person_id
""".strip()
    return f"""
WITH condition_events AS (
  {condition_sql}
),
measurement_events AS (
  {measurement_events_sql}
)
SELECT
  COALESCE(condition_events.person_id, measurement_summary.person_id) AS person_id,
  condition_events.condition_date,
  measurement_summary.measurement_date,
  measurement_summary.measurement_value
FROM condition_events
FULL OUTER JOIN (
  SELECT
    person_id,
    MIN(measurement_date) AS measurement_date,
    MAX(measurement_value) AS measurement_value
  FROM measurement_events
  GROUP BY person_id
) measurement_summary USING (person_id)
""".strip()


def render_baseline_sql(config: ProjectConfig) -> str:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    person_table = _qualify_table(cdr, config.phenotype.tables.person_table)
    observation_table = _qualify_table(cdr, config.phenotype.tables.observation_table)
    condition_table = _qualify_table(cdr, config.phenotype.tables.condition_table)
    concept_table = _qualify_table(cdr, config.phenotype.tables.concept_table)
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
    UPPER(COALESCE(source_concept.vocabulary_id, '')) AS source_vocabulary_id
  FROM `{condition_table}` co
  LEFT JOIN `{concept_table}` source_concept
    ON co.condition_source_concept_id = source_concept.concept_id
  WHERE co.{config.phenotype.condition_date_column} IS NOT NULL
),
denominator_mode AS (
  SELECT
    IF(
      COUNTIF(STARTS_WITH(source_vocabulary_id, 'ICD')) > 0,
      'icd_source_condition_records',
      'all_condition_records'
    ) AS omop_condition_record_source
  FROM condition_records
),
denominator_counts AS (
  SELECT
    condition_records.person_id,
    COUNT(DISTINCT condition_records.condition_date) AS omop_condition_record_dates
  FROM condition_records
  CROSS JOIN denominator_mode
  WHERE denominator_mode.omop_condition_record_source = 'all_condition_records'
    OR STARTS_WITH(condition_records.source_vocabulary_id, 'ICD')
  GROUP BY person_id
)
SELECT
  CAST(p.person_id AS STRING) AS person_id,
  observation_windows.obs_start_date,
  observation_windows.obs_end_date,
  CAST(NULL AS DATE) AS baseline_index_date,
  SAFE_CAST(p.year_of_birth AS FLOAT64) AS year_of_birth,
  SAFE_CAST(EXTRACT(YEAR FROM CURRENT_DATE()) - p.year_of_birth AS FLOAT64) AS age_raw,
  g.concept_name AS gender_concept_name,
  CASE
    WHEN LOWER(g.concept_name) = 'female' THEN 'female'
    WHEN LOWER(g.concept_name) = 'male' THEN 'male'
    WHEN g.concept_name IS NULL
      OR LOWER(g.concept_name) IN ('', 'no matching concept', 'none', 'unknown', 'skip', 'pmi: skip', 'prefer not to answer', 'pmi: prefer not to answer')
      THEN 'missing'
    ELSE 'other_or_unknown'
  END AS sex_category,
  CASE
    WHEN LOWER(g.concept_name) = 'female' THEN 1.0
    WHEN LOWER(g.concept_name) = 'male' THEN 0.0
    ELSE NULL
  END AS is_female,
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
  COALESCE(denominator_counts.omop_condition_record_dates, 0) AS omop_condition_record_dates,
  denominator_mode.omop_condition_record_source,
  COALESCE(denominator_counts.omop_condition_record_dates, 0) >= 2 AS eligible_ehr_denominator
FROM `{person_table}` p
LEFT JOIN observation_windows
  ON CAST(p.person_id AS STRING) = observation_windows.person_id
LEFT JOIN `{concept_table}` g
  ON p.gender_concept_id = g.concept_id
LEFT JOIN denominator_counts
  ON CAST(p.person_id AS STRING) = denominator_counts.person_id
CROSS JOIN denominator_mode
""".strip()


def _render_cofactor_cte(
    config: ProjectConfig,
    rule: ClinicalCofactorRule,
    *,
    cdr: str,
    index: int,
) -> str:
    condition_table = _qualify_table(cdr, config.phenotype.tables.condition_table)
    concept_table = _qualify_table(cdr, config.phenotype.tables.concept_table)
    predicate = _match_predicate(
        concept_id_column=f"co.{config.phenotype.condition_concept_column}",
        concept_name_column="condition_concept.concept_name",
        concept_ids=rule.condition_concept_ids,
        concept_terms=rule.condition_terms,
    )
    if predicate == "FALSE":
        return f"""
cofactor_{index} AS (
  SELECT CAST(NULL AS STRING) AS person_id, 0 AS {rule.name}
  WHERE FALSE
)
""".strip()
    return f"""
cofactor_{index} AS (
  SELECT
    CAST(co.{config.phenotype.person_id_column} AS STRING) AS person_id,
    1 AS {rule.name}
  FROM `{condition_table}` co
  LEFT JOIN `{concept_table}` condition_concept
    ON co.{config.phenotype.condition_concept_column} = condition_concept.concept_id
  WHERE {predicate}
  GROUP BY person_id
)
""".strip()


def render_clinical_cofactors_sql(config: ProjectConfig) -> str:
    if not config.phenotype.clinical_cofactors:
        return _empty_result_sql(columns=("CAST(NULL AS STRING) AS person_id",))
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    person_table = _qualify_table(cdr, config.phenotype.tables.person_table)
    ctes = [
        _render_cofactor_cte(config, rule, cdr=cdr, index=index)
        for index, rule in enumerate(config.phenotype.clinical_cofactors, start=1)
    ]
    joins = "\n".join(
        f"LEFT JOIN cofactor_{index} USING (person_id)"
        for index, _rule in enumerate(config.phenotype.clinical_cofactors, start=1)
    )
    columns = ",\n  ".join(
        f"COALESCE(cofactor_{index}.{rule.name}, 0) AS {rule.name}"
        for index, rule in enumerate(config.phenotype.clinical_cofactors, start=1)
    )
    ctes_sql = ",\n".join(ctes)
    joins = "\n".join(
        f"LEFT JOIN cofactor_{index}\n  ON person.person_id = cofactor_{index}.person_id"
        for index, _rule in enumerate(config.phenotype.clinical_cofactors, start=1)
    )
    return f"""
WITH
person AS (
  SELECT CAST(person_id AS STRING) AS person_id
  FROM `{person_table}`
),
{ctes_sql}
SELECT
  person.person_id,
  {columns}
FROM person
{joins}
""".strip()


def _render_cofactor_events_select(
    config: ProjectConfig,
    rule: ClinicalCofactorRule,
    *,
    cdr: str,
) -> str:
    condition_table = _qualify_table(cdr, config.phenotype.tables.condition_table)
    concept_table = _qualify_table(cdr, config.phenotype.tables.concept_table)
    predicate = _match_predicate(
        concept_id_column=f"co.{config.phenotype.condition_concept_column}",
        concept_name_column="condition_concept.concept_name",
        concept_ids=rule.condition_concept_ids,
        concept_terms=rule.condition_terms,
    )
    if predicate == "FALSE":
        return _empty_result_sql(
            columns=(
                "CAST(NULL AS STRING) AS person_id",
                "CAST(NULL AS STRING) AS cofactor",
                "CAST(NULL AS DATE) AS condition_date",
            )
        )
    name = _escape_term(rule.name)
    return f"""
SELECT
  CAST(co.{config.phenotype.person_id_column} AS STRING) AS person_id,
  '{name}' AS cofactor,
  DATE(co.{config.phenotype.condition_date_column}) AS condition_date
FROM `{condition_table}` co
LEFT JOIN `{concept_table}` condition_concept
  ON co.{config.phenotype.condition_concept_column} = condition_concept.concept_id
WHERE ({predicate})
  AND co.{config.phenotype.condition_date_column} IS NOT NULL
GROUP BY person_id, cofactor, condition_date
""".strip()


def render_clinical_cofactor_events_sql(config: ProjectConfig) -> str:
    if not config.phenotype.clinical_cofactors:
        return _empty_result_sql(
            columns=(
                "CAST(NULL AS STRING) AS person_id",
                "CAST(NULL AS STRING) AS cofactor",
                "CAST(NULL AS DATE) AS condition_date",
            )
        )
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    return "\nUNION ALL\n".join(
        _render_cofactor_events_select(config, rule, cdr=cdr)
        for rule in config.phenotype.clinical_cofactors
    )


def render_covariate_sql(config: ProjectConfig) -> str:
    return render_baseline_sql(config)


__all__ = [
    "render_baseline_sql",
    "render_case_tier_sql",
    "render_clinical_cofactor_events_sql",
    "render_clinical_cofactors_sql",
    "render_covariate_sql",
]
