from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import unittest
from unittest.mock import patch

import pandas as pd

from aou_workbench import preindex_profile
from aou_workbench.config import load_project_config
from aou_workbench.pipeline import build_cohort_artifacts
from aou_workbench.preindex_profile import (
    preindex_biomarker_path,
    preindex_condition_top_path,
    preindex_measurement_top_path,
    preindex_report_path,
    preindex_summary_path,
    profile_preindex_case_data,
)
from tests.support import build_demo_project_tree


class PreindexProfileTests(unittest.TestCase):
    def test_profile_writes_case_preindex_availability_outputs(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        effective, output_paths, cohort_df = build_cohort_artifacts(config)

        outputs = profile_preindex_case_data(
            effective,
            cohort_df,
            output_paths,
            windows=["365", "all"],
            top_n=10,
        )

        self.assertFalse(outputs["summary"].empty)
        self.assertFalse(outputs["condition_top"].empty)
        self.assertFalse(outputs["measurement_top"].empty)
        self.assertFalse(outputs["biomarker"].empty)

        for path in (
            preindex_summary_path(output_paths),
            preindex_condition_top_path(output_paths),
            preindex_measurement_top_path(output_paths),
            preindex_biomarker_path(output_paths),
            preindex_report_path(output_paths),
        ):
            self.assertTrue(Path(path).exists())

        condition_names = set(outputs["condition_top"]["concept_name"])
        self.assertIn("Type 2 diabetes mellitus", condition_names)
        self.assertNotIn("Rhabdomyolysis", condition_names)

        creatinine = outputs["biomarker"][
            (outputs["biomarker"]["biomarker"] == "creatinine")
            & (outputs["biomarker"]["window"] == "365d")
        ].iloc[0]
        self.assertEqual(creatinine["n_cases_with_measurement"], 2)

    def test_bigquery_profile_renders_sql_outputs(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        effective, output_paths, cohort_df = build_cohort_artifacts(config)
        bq_tables = replace(
            effective.phenotype.tables,
            cohort_table=None,
            condition_table="condition_occurrence",
            measurement_table="measurement",
            concept_table="concept",
        )
        bq_config = replace(effective, phenotype=replace(effective.phenotype, tables=bq_tables))
        sql_calls: list[str] = []
        query_results = iter(
            [
                pd.DataFrame(
                    [
                        {
                            "domain": "condition",
                            "window": "365d",
                            "n_cases": 3,
                            "n_cases_with_any": 2,
                            "pct_cases_with_any": 66.67,
                            "total_events": 2,
                            "median_events_per_case_with_any": 1.0,
                        }
                    ]
                ),
                pd.DataFrame(
                    [
                        {
                            "domain": "measurement",
                            "window": "365d",
                            "n_cases": 3,
                            "n_cases_with_any": 2,
                            "pct_cases_with_any": 66.67,
                            "total_events": 2,
                            "median_events_per_case_with_any": 1.0,
                        }
                    ]
                ),
                pd.DataFrame(
                    [
                        {
                            "domain": "condition",
                            "window": "365d",
                            "concept_id": "201",
                            "concept_name": "BigQuery hypertension",
                            "n_cases": 2,
                            "total_events": 2,
                            "median_days_before_index": 120,
                        }
                    ]
                ),
                pd.DataFrame(
                    [
                        {
                            "domain": "measurement",
                            "window": "365d",
                            "concept_id": "901",
                            "concept_name": "Creatinine",
                            "n_cases": 2,
                            "total_events": 2,
                            "median_days_before_index": 30,
                        }
                    ]
                ),
                pd.DataFrame(
                    [
                        {
                            "biomarker": "creatinine",
                            "window": "365d",
                            "n_cases": 3,
                            "n_cases_with_measurement": 2,
                            "pct_cases_with_measurement": 66.67,
                            "total_measurements": 2,
                            "median_latest_value": 1.0,
                        }
                    ]
                ),
            ]
        )

        def fake_query(sql: str) -> pd.DataFrame:
            sql_calls.append(sql)
            return next(query_results)

        with patch("aou_workbench.preindex_profile.query_bigquery_dataframe", side_effect=fake_query):
            outputs = profile_preindex_case_data(
                bq_config,
                cohort_df,
                output_paths,
                windows=["365", "all"],
                top_n=5,
            )

        self.assertEqual(len(sql_calls), 5)
        self.assertIn("observation_history", set(outputs["summary"]["domain"]))
        self.assertIn("condition", set(outputs["summary"]["domain"]))
        self.assertIn("measurement", set(outputs["summary"]["domain"]))

        sql_root = Path(preindex_profile.preindex_case_profile_root(output_paths)) / "sql"
        condition_sql = (sql_root / "top_conditions.sql").read_text(encoding="utf-8")
        condition_summary_sql = (sql_root / "condition_summary.sql").read_text(encoding="utf-8")
        biomarker_sql = (sql_root / "biomarker_availability.sql").read_text(encoding="utf-8")
        self.assertIn("`test-project.C2024Q3R9.condition_occurrence`", condition_sql)
        self.assertIn("STRUCT('365d' AS window_label, 365 AS days)", condition_sql)
        self.assertIn("window_label AS `window`", condition_summary_sql)
        self.assertNotIn(" AS window,", condition_sql)
        self.assertIn("ROW_NUMBER() OVER", condition_sql)
        self.assertIn("`test-project.C2024Q3R9.measurement`", biomarker_sql)
        self.assertIn("biomarker_terms", biomarker_sql)
        self.assertIn("creatine_kinase", biomarker_sql)

        report = Path(preindex_report_path(output_paths)).read_text(encoding="utf-8")
        self.assertIn("BigQuery hypertension", report)


if __name__ == "__main__":
    unittest.main()
