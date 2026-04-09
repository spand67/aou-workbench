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
- Stage 1 uses the AoU smaller-callset VCF workflow and writes a tiny derived genotype table for the matched cohort only.
- Larger analyses can still use Hail, but Hail is treated as the path for broader MT/VDS work rather than as a fallback for Stage 1.

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

Prepare and run Stage 1 with the AoU smaller-callset workflow:

```bash
aou-workbench preflight
aou-workbench prepare-stage1
aou-workbench run-stage1
```

`prepare-stage1` expects:

- `gsutil`
- `bcftools`
- requester-pays access through your AoU workspace project

If `bcftools` is unavailable, preflight will warn and Stage 1 will fail fast with a clear setup message instead of falling back to an unstable Hail extraction path.

Run the full pipeline after configuring stage-specific derived tables:

```bash
aou-workbench run-all
```

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
- Use smaller All of Us callsets whenever possible.
- Prefer VCF/PLINK-style exact extraction for tiny targeted analyses like Stage 1.
- Reserve Hail MT/VDS workflows for broader rare-variant and genome-scale analyses where Hail's distributed model is worth the added environment complexity.
- Keep row-level outputs in controlled-tier storage only; do not commit them.
- Treat `main` as the branch Workbench pulls for execution.
- Keep notebooks thin and delegate business logic to the package. Avoid leaving tracked notebooks open while pulling repo updates because Jupyter can autosave metadata changes.
