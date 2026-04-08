# Workbench Bootstrap

Use this repository as the Git source of truth and keep your All of Us Researcher Workbench copy as a disposable execution checkout.

## Recommended flow

1. Create or pull `spand67/aou-workbench` into your Workbench persistent disk.
2. Run `python -m pip install -e ".[workbench]"`.
3. Edit only configs and notebooks in Workbench if the change should be committed back here.
4. Run `aou-workbench preflight` before any stage execution.
5. Execute `aou-workbench run-all` or the stage-specific CLI commands.
6. Commit reusable code, configs, docs, and notebooks here. Do not commit row-level outputs.

## Bootstrap notebook

See `notebooks/01_workbench_bootstrap.ipynb` for an executable setup notebook that clones or pulls the repo, installs the package, and runs preflight.
