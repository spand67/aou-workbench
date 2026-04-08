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
    return f"SELECT {projection} WHERE FALSE"


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
    measurement_sql = (
        f"""
SELECT
  CAST(m.{config.phenotype.person_id_column} AS STRING) AS person_id,
  MIN(DATE(m.{config.phenotype.measurement_date_column})) AS measurement_date,
  MAX(SAFE_CAST(m.{config.phenotype.measurement_value_column} AS FLOAT64)) AS measurement_value
FROM `{measurement_table}` m
LEFT JOIN `{concept_table}` measurement_concept
  ON m.{config.phenotype.measurement_concept_column} = measurement_concept.concept_id
WHERE {measurement_predicate}
  {measurement_threshold}
GROUP BY person_id
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
    return f"""
WITH condition_events AS (
  {condition_sql}
),
measurement_events AS (
  {measurement_sql}
)
SELECT
  COALESCE(condition_events.person_id, measurement_events.person_id) AS person_id,
  condition_events.condition_date,
  measurement_events.measurement_date,
  measurement_events.measurement_value
FROM condition_events
FULL OUTER JOIN measurement_events USING (person_id)
""".strip()


def render_baseline_sql(config: ProjectConfig) -> str:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    person_table = _qualify_table(cdr, config.phenotype.tables.person_table)
    observation_table = _qualify_table(cdr, config.phenotype.tables.observation_table)
    concept_table = _qualify_table(cdr, config.phenotype.tables.concept_table)
    return f"""
WITH observation_windows AS (
  SELECT
    CAST(person_id AS STRING) AS person_id,
    MIN(DATE(observation_period_start_date)) AS obs_start_date,
    MAX(DATE(observation_period_end_date)) AS obs_end_date
  FROM `{observation_table}`
  GROUP BY person_id
)
SELECT
  CAST(p.person_id AS STRING) AS person_id,
  observation_windows.obs_start_date,
  observation_windows.obs_end_date,
  CAST(NULL AS DATE) AS baseline_index_date,
  SAFE_CAST(p.year_of_birth AS FLOAT64) AS year_of_birth,
  SAFE_CAST(EXTRACT(YEAR FROM CURRENT_DATE()) - p.year_of_birth AS FLOAT64) AS age_raw,
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
  CAST(NULL AS FLOAT64) AS pc10
FROM `{person_table}` p
LEFT JOIN observation_windows
  ON CAST(p.person_id AS STRING) = observation_windows.person_id
LEFT JOIN `{concept_table}` g
  ON p.gender_concept_id = g.concept_id
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


def render_covariate_sql(config: ProjectConfig) -> str:
    return render_baseline_sql(config)


__all__ = [
    "render_baseline_sql",
    "render_case_tier_sql",
    "render_clinical_cofactors_sql",
    "render_covariate_sql",
]
