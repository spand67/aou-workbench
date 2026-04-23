from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.config import load_project_config
from aou_workbench.cohort_summary import (
    cohort_summary_report_path,
    cohort_summary_table_path,
    summarize_clinical_demographics,
)
from aou_workbench.pipeline import build_cohort_artifacts, match_controls_artifacts
from aou_workbench.stage1_prepare import stage1_sample_manifest_path
from tests.support import build_demo_project_tree


class CohortSummaryTests(unittest.TestCase):
    def test_summary_writes_wgs_restricted_unmatched_and_matched_groups(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        effective, output_paths, cohort_df = build_cohort_artifacts(config)
        _, _, matched_df = match_controls_artifacts(effective, cohort_df)

        present = pd.DataFrame({"person_id": cohort_df["person_id"].astype(str)})
        present.to_csv(
            stage1_sample_manifest_path(effective.analysis.stage1.variant_table),
            sep="\t",
            index=False,
        )

        summary = summarize_clinical_demographics(effective, cohort_df, matched_df, output_paths)

        self.assertFalse(summary.empty)
        self.assertIn("unmatched_rhabdo_wgs", summary.columns)
        self.assertIn("unmatched_non_rhabdo_wgs", summary.columns)
        self.assertIn("matched_rhabdo_wgs", summary.columns)
        self.assertIn("matched_controls_wgs", summary.columns)
        self.assertTrue(Path(cohort_summary_table_path(output_paths)).exists())
        self.assertTrue(Path(cohort_summary_report_path(output_paths)).exists())

        n_row = summary[summary["variable"] == "N"].iloc[0]
        self.assertEqual(n_row["unmatched_rhabdo_wgs"], "4")
        self.assertEqual(n_row["matched_rhabdo_wgs"], "3")
        self.assertEqual(n_row["matched_controls_wgs"], "6")


if __name__ == "__main__":
    unittest.main()
