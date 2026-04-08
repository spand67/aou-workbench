WITH condition_events AS (
  SELECT
    person_id,
    MIN(condition_start_date) AS condition_date
  FROM `{{workspace_cdr}}.condition_occurrence`
  WHERE condition_concept_id IN ({{definite_condition_concept_ids}})
  GROUP BY person_id
),
measurement_events AS (
  SELECT
    person_id,
    MIN(measurement_date) AS measurement_date,
    MAX(value_as_number) AS max_value
  FROM `{{workspace_cdr}}.measurement`
  WHERE measurement_concept_id IN ({{ck_measurement_concept_ids}})
    AND value_as_number >= {{definite_ck_threshold}}
  GROUP BY person_id
)
SELECT
  condition_events.person_id,
  condition_date,
  measurement_date,
  max_value
FROM condition_events
INNER JOIN measurement_events USING (person_id)
