from __future__ import annotations

import unittest

from aou_workbench.config import load_project_config
from tests.support import build_demo_project_tree


class ConfigLoadingTests(unittest.TestCase):
    def test_load_project_config_reads_nested_yaml(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        self.assertEqual(config.analysis.analysis_name, "unit_test_rhabdo")
        self.assertEqual(len(config.panel.a_priori_variants), 6)
        self.assertEqual(config.phenotype.definite.measurement_min, 5000)
        self.assertTrue(config.config_hash)
