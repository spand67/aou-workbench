from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.config import load_project_config
from aou_workbench.matching import match_case_controls
from aou_workbench.paths import build_output_paths
from aou_workbench.regenie import prepare_regenie_inputs
from tests.support import build_demo_project_tree


class RegeniePreparationTests(unittest.TestCase):
    def test_prepare_regenie_inputs_writes_expected_files(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        cohort_df = build_rhabdo_cohort(config)
        matched_df = match_case_controls(cohort_df, config)
        output_paths = build_output_paths(config)

        outputs = prepare_regenie_inputs(config, matched_df, output_paths)

        for path in outputs.values():
            self.assertTrue(Path(path).exists(), path)

        matched_manifest = pd.read_csv(outputs["matched_samples"], sep="\t")
        gwas_manifest = pd.read_csv(outputs["gwas_samples"], sep="\t")
        keep = pd.read_csv(outputs["keep"], sep="\t")
        pheno = pd.read_csv(outputs["phenotypes"], sep="\t")
        covar = pd.read_csv(outputs["covariates"], sep="\t")
        covar_numeric = pd.read_csv(outputs["covariates_numeric"], sep="\t")
        commands = Path(outputs["commands"]).read_text(encoding="utf-8")

        self.assertGreaterEqual(len(matched_manifest), len(gwas_manifest))
        self.assertIn("FID", keep.columns)
        self.assertIn("IID", pheno.columns)
        self.assertIn("rhabdo_case", pheno.columns)
        self.assertIn("ancestry_pred", covar.columns)
        self.assertNotIn("ancestry_pred", covar_numeric.columns)
        self.assertIn("--step 1", commands)
        self.assertIn("--step 2", commands)
        self.assertIn("--catCovarList \"ancestry_pred\"", commands)
        self.assertIn("STEP2_BED_PREFIX_TEMPLATE", commands)
        self.assertIn("--bed", commands)


if __name__ == "__main__":
    unittest.main()
