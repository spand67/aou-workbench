SELECT
  CAST(p.person_id AS STRING) AS person_id,
  SAFE_CAST(EXTRACT(YEAR FROM CURRENT_DATE()) - p.year_of_birth AS FLOAT64) AS age,
  CASE
    WHEN LOWER(g.concept_name) = 'female' THEN 1.0
    WHEN LOWER(g.concept_name) = 'male' THEN 0.0
    ELSE NULL
  END AS is_female,
  ancestry.ancestry_pred,
  ancestry.pc1,
  ancestry.pc2,
  ancestry.pc3,
  ancestry.pc4,
  ancestry.pc5,
  ancestry.pc6,
  ancestry.pc7,
  ancestry.pc8,
  ancestry.pc9,
  ancestry.pc10
FROM `{{workspace_cdr}}.person` p
LEFT JOIN `{{workspace_cdr}}.concept` g
  ON p.gender_concept_id = g.concept_id
LEFT JOIN `{{workspace_cdr}}.person_ext` ancestry
  ON p.person_id = ancestry.person_id
