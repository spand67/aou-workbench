# Workbench Bootstrap

Use this repository as the Git source of truth and keep your All of Us Researcher Workbench copy as a disposable execution checkout.

## Recommended flow

1. Create or pull `spand67/aou-workbench` into your Workbench persistent disk.
2. Run `python -m pip install --user --no-deps -e .`.
3. Set `configs/workbench.yaml` to the attached non-`prep_` CDR if your workspace differs from the checked-in default.
4. Run `aou-workbench preflight` before any stage execution.
5. Start with `aou-workbench build-cohort` and `aou-workbench match-controls`.
6. Enable Stage 1-4 in `configs/rhabdo/analysis.yaml` only after their derived variant or genotype tables are configured.
7. Commit reusable code, configs, docs, and notebooks here. Do not commit row-level outputs.

## If a session breaks after installing Hail

Some Workbench images ship with `pyarrow` extensions compiled against NumPy 1.x. If pip upgrades NumPy to 2.x, commands can fail with `_ARRAY_API not found`.

Repair the session with:

```bash
python -m pip install --user --force-reinstall "numpy<2" "pandas<2.2" "scipy<1.12"
python -m pip install --user --no-deps --force-reinstall -e .
```

## Bootstrap notebook

See `notebooks/01_workbench_bootstrap.ipynb` for an executable setup notebook that clones or pulls the repo, installs the package, and runs preflight.

If `git pull` is blocked by notebook changes, close the notebook tabs in JupyterLab and stash them before pulling:

```bash
git stash push -m "temp notebook edits" notebooks/01_workbench_bootstrap.ipynb notebooks/02_rhabdo_end_to_end.ipynb
git pull --ff-only origin main
```
