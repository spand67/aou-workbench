# Rhabdomyolysis Phenotype Notes

The default workflow uses a tiered case definition:

- `broad`: at least one rhabdomyolysis condition record among denominator-eligible participants.
- `definite`: broad rhabdomyolysis plus CK >5000 from 7 days before through 45 days after the first qualifying rhabdomyolysis condition date.
- `control`: denominator-eligible participants with no rhabdomyolysis condition record and no CK >5000 measurement. Missing CK is allowed.
- `indeterminate_ck_only`: CK >5000 without rhabdomyolysis diagnosis; excluded from primary case-control analyses.

The study denominator requires at least 2 distinct OMOP `condition_occurrence` dates. When ICD-derived source concepts are available, the denominator uses ICD-derived condition records; otherwise it falls back to all OMOP condition records and reports that fallback in cohort outputs.

The implementation is intentionally config-driven:

- Diagnosis concepts live in `configs/rhabdo/phenotype.yaml`.
- Measurement concept IDs, thresholds, and directional timing windows are editable per project.
- Observation windows, exclusion concepts, and clinical cofactors are all configurable.

Primary analyses use `broad` cases. `definite` cases are preserved in the built cohort for sensitivity analyses.

Sepsis, crush injury, and acute renal injury are retained as clinical cofactors rather than primary exclusions. The cohort builder writes each configured cofactor as an ever/never flag and as time-anchored flags relative to the case or inherited-control index date:

- `preindex_*`: any matching condition before the peri-index window.
- `periindex_*`: any matching condition from 7 days before through 45 days after index.
- `postindex_*`: any matching condition after the peri-index window.

Pre-index flags are the preferred baseline covariates for matched prediction and genetic models. Peri-index sepsis and renal injury should be interpreted as acute-event context, severity, or sensitivity-stratification variables because they may lie on the rhabdomyolysis causal/ascertainment pathway.
