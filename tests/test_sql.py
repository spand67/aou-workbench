from __future__ import annotations

import unittest

from aou_workbench.config import load_project_config
from aou_workbench.phenotype_sql import (
    _empty_result_sql,
    render_case_tier_sql,
    render_clinical_cofactor_events_sql,
    render_covariate_sql,
)
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
        self.assertIn("condition_concept_id IN (100)", definite_sql)
        self.assertIn("BETWEEN -7 AND 45", definite_sql)
        self.assertIn("person", covariate_sql)
        self.assertIn("omop_condition_record_dates", covariate_sql)
        self.assertIn("condition_source_concept_id", covariate_sql)
        self.assertIn("eligible_ehr_denominator", covariate_sql)
        self.assertIn("gender_concept_name", covariate_sql)
        self.assertIn("sex_category", covariate_sql)
        cofactor_sql = render_clinical_cofactor_events_sql(config)
        self.assertIn("'sepsis' AS cofactor", cofactor_sql)
        self.assertIn("DATE(co.condition_start_date) AS condition_date", cofactor_sql)
        self.assertIn("UNION ALL", cofactor_sql)

    def test_empty_result_sql_is_valid_for_condition_only_tiers(self) -> None:
        sql = _empty_result_sql(columns=("CAST(NULL AS STRING) AS person_id",))
        self.assertIn("FROM (SELECT 1 AS _unused)", sql)
        self.assertIn("WHERE FALSE", sql)
