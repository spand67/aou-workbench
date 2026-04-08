"""SQL rendering helpers for Workbench phenotype and covariate extraction."""

from __future__ import annotations

from .config import CaseTierRule, ProjectConfig


def _concept_sql(
    ids: tuple[int, ...],
    *,
    cdr: str,
    table: str,
    concept_column: str,
    date_column: str,
) -> str:
    if not ids:
        return "SELECT CAST(NULL AS INT64) AS person_id WHERE FALSE"
    joined = ", ".join(str(value) for value in ids)
    return (
        f"SELECT person_id, MIN({date_column}) AS event_date "
        f"FROM `{cdr}.{table}` "
        f"WHERE {concept_column} IN ({joined}) "
        f"GROUP BY person_id"
    )


def render_case_tier_sql(config: ProjectConfig, tier: CaseTierRule) -> str:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    condition_table = "condition_occurrence"
    measurement_table = "measurement"
    condition_sql = _concept_sql(
        tier.condition_concept_ids,
        cdr=cdr,
        table=condition_table,
        concept_column=config.phenotype.condition_concept_column,
        date_column=config.phenotype.condition_date_column,
    )
    measurement_ids = ", ".join(str(value) for value in tier.measurement_concept_ids) or "NULL"
    measurement_predicate = (
        f"{config.phenotype.measurement_concept_column} IN ({measurement_ids})"
        if tier.measurement_concept_ids
        else "FALSE"
    )
    measurement_threshold = (
        f"AND value_as_number >= {tier.measurement_min}"
        if tier.measurement_min is not None
        else ""
    )
    return f"""
WITH condition_events AS (
  {condition_sql}
),
measurement_events AS (
  SELECT
    person_id,
    MIN({config.phenotype.measurement_date_column}) AS measurement_date,
    MAX({config.phenotype.measurement_value_column}) AS max_value
  FROM `{cdr}.{measurement_table}`
  WHERE {measurement_predicate}
    {measurement_threshold}
  GROUP BY person_id
)
SELECT
  person_id,
  condition_events.event_date AS condition_date,
  measurement_events.measurement_date,
  measurement_events.max_value
FROM condition_events
FULL OUTER JOIN measurement_events USING (person_id)
""".strip()


def render_covariate_sql(config: ProjectConfig) -> str:
    cdr = config.workbench.workspace_cdr or "{{workspace_cdr}}"
    pcs = ",\n    ".join(
        f"SAFE_CAST(pc{index} AS FLOAT64) AS pc{index}" for index in range(1, 11)
    )
    return f"""
SELECT
  CAST(p.person_id AS STRING) AS person_id,
  SAFE_CAST(EXTRACT(YEAR FROM CURRENT_DATE()) - p.year_of_birth AS FLOAT64) AS age,
  CASE
    WHEN LOWER(g.concept_name) = 'female' THEN 1.0
    WHEN LOWER(g.concept_name) = 'male' THEN 0.0
    ELSE NULL
  END AS is_female,
  ancestry.ancestry_pred,
  {pcs}
FROM `{cdr}.person` p
LEFT JOIN `{cdr}.concept` g
  ON p.gender_concept_id = g.concept_id
LEFT JOIN `{cdr}.person_ext` ancestry
  ON p.person_id = ancestry.person_id
""".strip()


__all__ = ["render_case_tier_sql", "render_covariate_sql"]
