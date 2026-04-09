from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.config import load_project_config
from aou_workbench.paths import build_output_paths
from aou_workbench.stage2_plp_panel import run_stage2_plp_panel
from aou_workbench.stage2_prepare import stage2_sample_manifest_path
from tests.support import build_demo_project_tree


class Stage2PlpPanelTests(unittest.TestCase):
    def test_stage2_runs_case_vs_rest_comparisons_using_manifest_filtered_samples(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        cohort_df = build_rhabdo_cohort(config)
        pd.DataFrame({"person_id": cohort_df["person_id"].astype(str)}).to_csv(
            stage2_sample_manifest_path(paths["stage2_table"]),
            sep="\t",
            index=False,
        )

        variant_df, gene_df, person_df = run_stage2_plp_panel(config, cohort_df, build_output_paths(config))

        self.assertFalse(variant_df.empty)
        self.assertFalse(gene_df.empty)
        self.assertFalse(person_df.empty)
        self.assertIn("omop_rhabdo_vs_non_rhabdo", set(gene_df["comparison"]))
        self.assertIn("omop_rhabdo_plus_ck_vs_non_rhabdo", set(gene_df["comparison"]))
        self.assertTrue(set(gene_df["gene"]).issubset(set(config.panel.genes_of_interest)))

    def test_stage2_requires_sample_manifest(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        cohort_df = build_rhabdo_cohort(config)
        manifest_path = Path(stage2_sample_manifest_path(paths["stage2_table"]))
        if manifest_path.exists():
            manifest_path.unlink()

        with self.assertRaises(RuntimeError):
            run_stage2_plp_panel(config, cohort_df, build_output_paths(config))


if __name__ == "__main__":
    unittest.main()
