# aou-workbench

`aou-workbench` is the source-of-truth repository for reusable All of Us genomics analysis code, configs, SQL, and runbooks. The local Git repo is meant to be pushed to `spand67/aou-workbench`, then pulled into All of Us Researcher Workbench as the execution copy.

The first implemented workflow is a rhabdomyolysis program with:

1. A priori candidate variant analysis.
2. Pathogenic and likely pathogenic variant analysis in genes of interest.
3. Gene burden analysis.
4. Common-variant GWAS with covariate-adjusted summaries.

The default Workbench path is now cohort-first and BigQuery-native:

- `build-cohort` and `match-controls` read directly from the attached AoU CDR tables in Researcher Workbench.
- The checked-in starter config targets the non-`prep_` clinical dataset and the AoU genomics ancestry TSV in the attached controlled bucket.
- Stage 1 now uses direct Hail/VDS extraction for the exact a priori variant panel and writes a tiny derived genotype table for the matched cohort only.
- Hail is therefore part of the default Workbench path for Stage 1 exact variants and for broader MT/VDS analyses.

## Quickstart

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
python -m unittest discover -s tests
```

Install the optional Workbench extras inside a Hail-capable notebook environment:

```bash
python -m pip install -e ".[workbench]"
```

On All of Us Researcher Workbench, prefer this flow for the package itself:

```bash
python -m pip install --user --no-deps -e .
```

Install the separate Hail extra only on environments that do not already ship Hail and only when you truly need Hail-backed larger analyses:

```bash
python -m pip install -e ".[workbench,hail]"
```

If a Workbench session already pulled in `numpy>=2` and started failing with `pyarrow` or `_ARRAY_API` errors, repair the environment with:

```bash
python -m pip install --user --force-reinstall "numpy<2" "pandas<2.2" "scipy<1.12"
python -m pip install --user --no-deps --force-reinstall -e .
```

## Common commands

Run preflight checks:

```bash
aou-workbench preflight
```

Build the tiered rhabdomyolysis cohort from BigQuery:

```bash
aou-workbench build-cohort
aou-workbench match-controls
```

For a genomics-first rerun, restrict every cohort artifact directly to WGS-available participants:

```bash
aou-workbench build-cohort --require-wgs
aou-workbench match-controls --require-wgs
aou-workbench characterize-cohort --require-wgs
```

The `--require-wgs` commands filter directly in BigQuery with `cb_search_person.has_whole_genome_variant = 1`, then save a WGS-restricted built cohort, matched cohort, CONSORT, Table 1, split summaries, and clinical model input. Hail/ACAF MatrixTable overlap is checked later by the GWAS runner, but the cohort denominator starts from the CDR WGS flag. This is the preferred setup before any GWAS or PRS work. `prepare-wgs-manifest` remains available only as an optional audit/debug command.

Characterize the inclusive case-control cohort before genomic modeling:

```bash
aou-workbench characterize-cohort
```

This writes CONSORT counts, a matched clinical Table 1 with effect sizes and 95% confidence intervals, a matched-group train/test split summary, model eligibility counts, split-specific Table 1 outputs, a missingness summary, a clinical model input table, and sepsis/renal-injury summaries under `cohort/` in the configured run root. Table 1 reports sex as explicit female, male, other/unknown, and missing categories while retaining `is_female` as the binary matching/modeling indicator; ancestry missingness is reported as an explicit category. Clinical cofactors are retained as ever/never flags and split into `remote_preindex_*` (>30 days before index), `near_preindex_*` (8-30 days before index), `preindex_*`, `periindex_*`, and `postindex_*` columns using the `[-7, +45]` day peri-index window. The `case_cofactor_prior_timing.tsv` output summarizes the nearest sepsis or renal-injury condition record on or before the first rhabdomyolysis index date among cases, using same-day, 1-7, 8-30, 31-90, 91-365, and >365 day bins; `case_cofactor_prior_timing_histogram.svg` plots the same binned distribution. Use remote pre-index columns as baseline covariates for matched models; treat near-pre-index, peri-index, and post-index sepsis/renal injury as acute-event characterization or sensitivity-stratification variables. The split is assigned at the matched-group level with an 80/20 train/test holdout and 5 training-only CV folds, so final model performance can be evaluated once on the holdout while clinical model selection and later PRS tuning happen inside the training folds. The primary clinical modeling flag excludes pre/peri-index crush injury at the row and matched-case group level, while sensitivity flags retain trauma or additionally remove peri-index sepsis groups. The command also refreshes `rhabdo_dashboard.html`, a static dashboard with formatted cards, figures, and tables; `rhabdo_summary_report.md` remains available as a lightweight Markdown digest.

Run the first clinical-only prediction model:

```bash
aou-workbench run-clinical-model
```

This trains an L2-regularized logistic model on the primary non-traumatic training set and evaluates it once on the held-out test set. The model uses a deliberately simple baseline feature set: age, observation depth, condition-record depth, sex category, and ancestry category. Sepsis and renal injury are retained for cohort characterization and sensitivity interpretation, but are not used as clinical model predictors. Metrics, coefficients, predictions, calibration, ROC, and precision-recall outputs are written under `clinical/model/`.

Run the first reduced-marker Hail GWAS pilot:

```bash
aou-workbench run-hail-pilot-gwas \
  --chromosomes 22 \
  --min-maf 0.05 \
  --min-mac 20 \
  --min-call-rate 0.98 \
  --hwe-p-control 1e-6 \
  --analysis-split train \
  --eligibility-flag primary_model_eligible \
  --label acaf_chr22_maf05_train_qc
```

This uses the AoU ACAF threshold split MatrixTable, not the full variant database. It restricts to the training split and primary model-eligible rows, then computes variant QC inside that analysis sample before association testing. The pilot keeps autosomal biallelic SNPs with MAF >= 0.05, minor allele count >= 20, call rate >= 0.98, and control-only Hardy-Weinberg equilibrium p >= 1e-6. Outputs are isolated under `stage4/hail_pilot/<label>/` and include GWAS results, lead hits, QC JSON, sequential variant QC counts, Manhattan and QQ plots, and a markdown report.

Rebuild the cross-analysis dashboard and Markdown digest from existing outputs without rerunning queries or models:

```bash
aou-workbench report
```

Prepare and run Stage 1 with the direct WGS VDS workflow:

```bash
aou-workbench preflight
aou-workbench prepare-stage1
aou-workbench run-stage1
```

`prepare-stage1` expects:

- `gsutil`
- `hail`
- access to the AoU WGS VDS path configured in `configs/workbench.yaml`
- requester-pays access through your AoU workspace project

If a Workbench session already ships Hail, reinstall the repo without dependency upgrades:

```bash
python -m pip install --user --no-deps --force-reinstall -e .
```

If Hail is unavailable, preflight will warn and Stage 1 will fail fast with a clear setup message instead of falling back to a slower ad hoc extraction path.

Run the full pipeline after configuring stage-specific derived tables:

```bash
aou-workbench run-all
```

Profile pre-index clinical data availability for the primary case tier:

```bash
aou-workbench profile-preindex-cases --windows 365,1095,all --top-n 25
```

This writes `availability_summary.tsv`, `top_conditions.tsv`, `top_measurements.tsv`, `biomarker_availability.tsv`, generated SQL, and `report.md` under `clinical/preindex_case_profile/` in the configured run root. Treat the entire folder as controlled-tier output: the aggregate tables and report summarize case history, and the generated SQL embeds case IDs and index dates.

Write the consolidated markdown report:

```bash
aou-workbench report
```

## Project layout

- `src/aou_workbench/`: reusable package code.
- `configs/`: Workbench, phenotype, cohort, panel, and analysis configs.
- `docs/`: phenotype definitions, runbooks, and Git workflow notes.
- `notebooks/`: lightweight notebook entrypoints for Researcher Workbench.
- `sql/`: BigQuery SQL templates for phenotype and covariate extraction.
- `tests/`: unit and synthetic integration coverage.
- `outputs/`: local outputs, ignored by Git.

## Design notes

- Assume execution happens inside All of Us Researcher Workbench.
- Default cohort inputs come from `person`, `observation_period`, `condition_occurrence`, `measurement`, and `person_ext` in the attached CDR.
- Use direct WGS VDS extraction for tiny exact-site panels when it avoids broader callset scanning.
- Reserve broader Hail MT/VDS workflows for rare-variant and genome-scale analyses where Hail's distributed model is worth the added environment complexity.
- Keep row-level outputs in controlled-tier storage only; do not commit them.
- Treat `main` as the branch Workbench pulls for execution.
- Keep notebooks thin and delegate business logic to the package. Avoid leaving tracked notebooks open while pulling repo updates because Jupyter can autosave metadata changes.
