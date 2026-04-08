from __future__ import annotations

import unittest

from aou_workbench.config import load_project_config
from aou_workbench.phenotype_sql import render_case_tier_sql, render_covariate_sql
from tests.support import build_demo_project_tree


class SqlRenderingTests(unittest.TestCase):
    def test_rendered_sql_uses_configured_thresholds(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        definite_sql = render_case_tier_sql(config, config.phenotype.definite)
        covariate_sql = render_covariate_sql(config)
        self.assertIn("5000", definite_sql)
        self.assertIn("condition_occurrence", definite_sql)
        self.assertIn("person", covariate_sql)
