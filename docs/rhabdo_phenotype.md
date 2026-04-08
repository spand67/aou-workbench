# Rhabdomyolysis Phenotype Notes

The default workflow uses a tiered case definition:

- `definite`: diagnosis evidence plus CK or other measurement evidence above the strict threshold within the configured lookback window.
- `probable`: broader diagnosis or measurement evidence using a lower measurement threshold.
- `control`: no case evidence and adequate observation time.

The implementation is intentionally config-driven:

- Diagnosis concepts live in `configs/rhabdo/phenotype.yaml`.
- Measurement concept IDs and thresholds are editable per project.
- Observation windows, exclusion concepts, and clinical cofactors are all configurable.

Primary analyses use `definite` cases. `probable` cases are preserved in the built cohort for sensitivity analyses.
