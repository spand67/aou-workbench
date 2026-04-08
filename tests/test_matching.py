from __future__ import annotations

import unittest

from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.config import load_project_config
from aou_workbench.matching import match_case_controls
from tests.support import build_demo_project_tree


class MatchingTests(unittest.TestCase):
    def test_matching_builds_primary_case_control_sets(self) -> None:
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
        self.assertGreaterEqual(int((cohort_df["case_tier"] == "definite").sum()), 3)
        self.assertGreater(int((matched_df["analysis_case"] == 1).sum()), 0)
        self.assertGreater(int((matched_df["analysis_case"] == 0).sum()), 0)
        self.assertIn("match_group_id", matched_df.columns)
