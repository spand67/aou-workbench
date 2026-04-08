# Git Workflow

`aou-workbench` is meant to become the canonical GitHub repository at `spand67/aou-workbench`.

## Working model

- Do development locally in this repository.
- Push reusable changes to `main`.
- Pull `main` inside All of Us Researcher Workbench before running analyses.
- Keep output data in controlled-tier storage and out of Git.

## Suggested bootstrap once the remote exists

```bash
git remote add origin git@github.com:spand67/aou-workbench.git
git add .
git commit -m "Initial AoU workbench scaffold"
git push -u origin main
```
