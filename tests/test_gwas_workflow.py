from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.config import load_project_config
from aou_workbench.gwas_workflow import prepare_terminal_gwas_workspace
from aou_workbench.matching import match_case_controls
from aou_workbench.paths import build_output_paths
from tests.support import build_demo_project_tree


class GwasWorkflowPreparationTests(unittest.TestCase):
    def test_prepare_terminal_gwas_workspace_writes_scripts_and_manifests(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        config = replace(
            config,
            workbench=replace(
                config.workbench,
                workspace_bucket="gs://unit-test-bucket",
                max_unrelated_path=paths["max_unrelated"],
            ),
        )
        cohort_df = build_rhabdo_cohort(config)
        matched_df = match_case_controls(cohort_df, config)
        output_paths = build_output_paths(config)

        outputs = prepare_terminal_gwas_workspace(config, matched_df, output_paths)

        for key in (
            "gwas_workspace",
            "matched_manifest",
            "restricted_manifest",
            "rewrite_script",
            "prune_script",
            "step1_script",
            "step2_script",
            "merge_script",
            "sync_script",
            "dsub_step1",
            "dsub_step2",
            "readme",
        ):
            self.assertTrue(Path(outputs[key]).exists(), key)

        restricted = pd.read_csv(outputs["restricted_manifest"], sep="\t")
        self.assertFalse(restricted.empty)
        self.assertNotIn("11", set(restricted["IID"].astype(str)))
        self.assertNotIn("12", set(restricted["IID"].astype(str)))

        rewrite_script = Path(outputs["rewrite_script"]).read_text(encoding="utf-8")
        step2_script = Path(outputs["step2_script"]).read_text(encoding="utf-8")
        dsub_step2 = Path(outputs["dsub_step2"]).read_text(encoding="utf-8")

        self.assertIn("--max-alleles 2", rewrite_script)
        self.assertIn("--maf 0.01", rewrite_script)
        self.assertIn("--pred \"$OUTDIR/step1_pred.list\"", step2_script)
        self.assertIn("us-central1", dsub_step2)
        self.assertIn("gs://unit-test-bucket", dsub_step2)


if __name__ == "__main__":
    unittest.main()
