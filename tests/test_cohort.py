from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.cohort import build_rhabdo_cohort, cohort_qc_summary
from aou_workbench.config import load_project_config
from tests.support import build_demo_project_tree


class CohortPhenotypeDefinitionTests(unittest.TestCase):
    def _config(self):
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        return paths, config

    def test_broad_definite_control_and_ck_only_definitions(self) -> None:
        _, config = self._config()
        cohort = build_rhabdo_cohort(config).set_index("person_id")

        self.assertEqual(cohort.loc["1", "case_tier"], "definite")
        self.assertEqual(cohort.loc["2", "case_tier"], "definite")
        self.assertEqual(cohort.loc["3", "case_tier"], "definite")
        self.assertEqual(cohort.loc["4", "case_tier"], "broad")
        self.assertEqual(cohort.loc["5", "case_tier"], "control")
        self.assertEqual(cohort.loc["15", "case_tier"], "indeterminate_ck_only")

        self.assertEqual(int(cohort.loc["4", "broad_rhabdo_case"]), 1)
        self.assertEqual(int(cohort.loc["4", "definite_rhabdo_case"]), 0)
        self.assertTrue(bool(cohort.loc["5", "eligible_control"]))
        self.assertFalse(bool(cohort.loc["15", "eligible_control"]))

    def test_definite_ck_window_is_directional(self) -> None:
        paths, config = self._config()
        measurements = pd.read_csv(paths["measurement_table"], sep="\t")
        measurements.loc[measurements["person_id"] == 2, "measurement_date"] = "2022-02-01"
        measurements.to_csv(paths["measurement_table"], sep="\t", index=False)

        cohort = build_rhabdo_cohort(config).set_index("person_id")

        self.assertEqual(cohort.loc["2", "case_tier"], "broad")
        self.assertEqual(int(cohort.loc["2", "definite_rhabdo_case"]), 0)

    def test_denominator_requires_two_distinct_condition_dates(self) -> None:
        paths, config = self._config()
        conditions = pd.read_csv(paths["condition_table"], sep="\t")
        conditions = conditions[
            ~(
                (conditions["person_id"] == 4)
                & (conditions["condition_concept_name"] == "Myalgia")
            )
        ]
        conditions.to_csv(paths["condition_table"], sep="\t", index=False)

        cohort = build_rhabdo_cohort(config).set_index("person_id")

        self.assertEqual(cohort.loc["4", "case_tier"], "excluded_denominator")
        self.assertEqual(int(cohort.loc["4", "omop_condition_record_dates"]), 1)
        self.assertFalse(bool(cohort.loc["4", "eligible_ehr_denominator"]))

    def test_legacy_probable_config_maps_to_broad(self) -> None:
        paths, _ = self._config()
        cohort_path = Path(paths["cohort"])
        text = cohort_path.read_text(encoding="utf-8")
        text = text.replace("primary_case_tier: broad", "primary_case_tier: probable")
        text = text.replace("- definite", "- probable")
        cohort_path.write_text(text, encoding="utf-8")

        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        self.assertEqual(config.cohort.primary_case_tier, "broad")
        self.assertEqual(config.cohort.sensitivity_case_tiers, ("broad",))

    def test_qc_summary_reports_denominator_and_exclusions(self) -> None:
        _, config = self._config()
        cohort = build_rhabdo_cohort(config)
        qc = cohort_qc_summary(cohort)

        self.assertEqual(qc["broad_rhabdo_cases"], 4)
        self.assertEqual(qc["definite_rhabdo_cases"], 3)
        self.assertEqual(qc["indeterminate_ck_only"], 1)
        self.assertGreater(qc["eligible_ehr_denominator"], 0)


if __name__ == "__main__":
    unittest.main()
