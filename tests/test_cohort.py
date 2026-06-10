from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.cohort import apply_time_anchored_clinical_cofactors, build_rhabdo_cohort, cohort_qc_summary
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

    def test_sex_category_preserves_missing_and_other_values(self) -> None:
        paths, config = self._config()
        cohort_input = pd.read_csv(paths["cohort_table"], sep="\t")
        cohort_input.loc[cohort_input["person_id"] == 14, "sex"] = ""
        cohort_input.loc[cohort_input["person_id"] == 15, "sex"] = "not male or female"
        cohort_input.to_csv(paths["cohort_table"], sep="\t", index=False)

        cohort = build_rhabdo_cohort(config).set_index("person_id")

        self.assertEqual(cohort.loc["1", "sex_category"], "female")
        self.assertEqual(cohort.loc["2", "sex_category"], "male")
        self.assertEqual(cohort.loc["14", "sex_category"], "missing")
        self.assertEqual(cohort.loc["15", "sex_category"], "other_or_unknown")
        self.assertTrue(pd.isna(cohort.loc["14", "is_female"]))
        self.assertTrue(pd.isna(cohort.loc["15", "is_female"]))

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

    def test_time_anchored_cofactors_split_preindex_and_periindex_events(self) -> None:
        _, config = self._config()
        cohort = build_rhabdo_cohort(config).set_index("person_id")

        self.assertEqual(int(cohort.loc["1", "sepsis"]), 1)
        self.assertEqual(int(cohort.loc["1", "periindex_sepsis"]), 1)
        self.assertEqual(int(cohort.loc["1", "preindex_sepsis"]), 0)
        self.assertEqual(int(cohort.loc["2", "renal_injury"]), 1)
        self.assertEqual(int(cohort.loc["2", "periindex_renal_injury"]), 1)
        self.assertEqual(int(cohort.loc["2", "preindex_renal_injury"]), 0)
        self.assertEqual(int(cohort.loc["3", "preindex_sepsis"]), 1)
        self.assertEqual(int(cohort.loc["3", "preindex_renal_injury"]), 1)

    def test_matched_controls_are_time_anchored_to_inherited_case_index_date(self) -> None:
        _, config = self._config()
        matched = pd.DataFrame(
            [
                {"person_id": "5", "index_date": "2022-01-10", "analysis_case": 0},
                {"person_id": "5", "index_date": "2021-01-01", "analysis_case": 0},
            ]
        )
        events = pd.DataFrame(
            [
                {"person_id": "5", "cofactor": "sepsis", "condition_date": "2022-01-11"},
                {"person_id": "5", "cofactor": "renal_injury", "condition_date": "2020-06-01"},
            ]
        )

        anchored = apply_time_anchored_clinical_cofactors(config, matched, events)

        self.assertEqual(int(anchored.loc[0, "periindex_sepsis"]), 1)
        self.assertEqual(int(anchored.loc[0, "preindex_renal_injury"]), 1)
        self.assertEqual(int(anchored.loc[1, "postindex_sepsis"]), 1)
        self.assertEqual(int(anchored.loc[1, "preindex_renal_injury"]), 1)

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
