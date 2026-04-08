from __future__ import annotations

import unittest

from aou_workbench.config import load_project_config
from aou_workbench.paths import build_output_paths
from tests.support import build_demo_project_tree


class PathPlanningTests(unittest.TestCase):
    def test_build_output_paths_uses_analysis_slug(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        output_paths = build_output_paths(config)
        self.assertIn("unit-test-rhabdo", output_paths.run_root)
        self.assertTrue(output_paths.stage4_manhattan_svg.endswith("manhattan.svg"))
