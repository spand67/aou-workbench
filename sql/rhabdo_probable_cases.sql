WITH probable_measurements AS (
  SELECT
    person_id,
    MIN(measurement_date) AS measurement_date,
    MAX(value_as_number) AS max_value
  FROM `{{workspace_cdr}}.measurement`
  WHERE measurement_concept_id IN ({{ck_measurement_concept_ids}})
    AND value_as_number >= {{probable_ck_threshold}}
  GROUP BY person_id
),
probable_conditions AS (
  SELECT
    person_id,
    MIN(condition_start_date) AS condition_date
  FROM `{{workspace_cdr}}.condition_occurrence`
  WHERE condition_concept_id IN ({{probable_condition_concept_ids}})
  GROUP BY person_id
)
SELECT
  COALESCE(probable_measurements.person_id, probable_conditions.person_id) AS person_id,
  probable_conditions.condition_date,
  probable_measurements.measurement_date,
  probable_measurements.max_value
FROM probable_measurements
FULL OUTER JOIN probable_conditions USING (person_id)
