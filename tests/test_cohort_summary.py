from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.config import load_project_config
from aou_workbench.cohort_summary import (
    add_model_eligibility_flags,
    add_model_splits,
    build_case_cofactor_prior_timing_summary,
    build_matched_table1,
    case_cofactor_prior_timing_path,
    characterize_case_control_cohort,
    clinical_model_input_path,
    clinical_characterization_report_path,
    consort_counts_path,
    critical_illness_summary_path,
    cohort_summary_report_path,
    cohort_summary_table_path,
    matched_table1_path,
    model_eligibility_summary_path,
    missingness_summary_path,
    model_split_summary_path,
    split_table1_path,
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
        self.assertIn("unmatched_smd", summary.columns)
        self.assertIn("matched_p_value", summary.columns)
        self.assertTrue(Path(cohort_summary_table_path(output_paths)).exists())
        self.assertTrue(Path(cohort_summary_report_path(output_paths)).exists())

        n_row = summary[summary["variable"] == "N"].iloc[0]
        self.assertEqual(n_row["unmatched_rhabdo_wgs"], "4")
        self.assertEqual(n_row["matched_rhabdo_wgs"], "3")
        self.assertEqual(n_row["matched_controls_wgs"], "6")

        age_row = summary[summary["variable"] == "Age at index, mean (SD)"].iloc[0]
        self.assertNotEqual(age_row["unmatched_smd"], "")
        self.assertNotEqual(age_row["matched_p_value"], "")

    def test_characterization_writes_consort_table1_and_critical_illness_outputs(self) -> None:
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

        outputs = characterize_case_control_cohort(effective, cohort_df, matched_df, output_paths)

        self.assertFalse(outputs["consort"].empty)
        self.assertFalse(outputs["table1"].empty)
        self.assertFalse(outputs["critical_illness"].empty)
        self.assertTrue(Path(consort_counts_path(output_paths)).exists())
        self.assertTrue(Path(matched_table1_path(output_paths)).exists())
        self.assertTrue(Path(critical_illness_summary_path(output_paths)).exists())
        self.assertTrue(Path(case_cofactor_prior_timing_path(output_paths)).exists())
        self.assertTrue(Path(split_table1_path(output_paths)).exists())
        self.assertTrue(Path(model_split_summary_path(output_paths)).exists())
        self.assertTrue(Path(model_eligibility_summary_path(output_paths)).exists())
        self.assertTrue(Path(missingness_summary_path(output_paths)).exists())
        self.assertTrue(Path(clinical_model_input_path(output_paths)).exists())
        self.assertTrue(Path(clinical_characterization_report_path(output_paths)).exists())
        self.assertIn("periindex_sepsis", set(outputs["critical_illness"]["variable"]))
        self.assertFalse(outputs["case_cofactor_prior_timing"].empty)
        self.assertIn("ci_95", outputs["table1"].columns)
        self.assertIn("missing", set(outputs["table1"]["variable"].str.extract(r"Ancestry: (.*), n", expand=False).dropna()))
        self.assertIn("train", set(outputs["split_summary"]["group"]))
        self.assertIn("test", set(outputs["split_summary"]["group"]))
        self.assertIn("primary_nontraumatic", set(outputs["eligibility"]["eligibility_set"]))
        self.assertIn("age_at_index", set(outputs["missingness"]["variable"]))
        self.assertIn("analysis_split", outputs["clinical_model_input"].columns)
        self.assertIn("primary_model_eligible", outputs["clinical_model_input"].columns)

    def test_case_cofactor_prior_timing_bins_nearest_prior_events(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        effective, _, cohort_df = build_cohort_artifacts(config)
        _, _, matched_df = match_controls_artifacts(effective, cohort_df)

        timing = build_case_cofactor_prior_timing_summary(effective, cohort_df, matched_df)
        built = timing[timing["group"] == "built_broad_cases"]
        renal_near = built[
            (built["cofactor"] == "renal_injury")
            & (built["window"] == "1_to_7_days_before")
        ].iloc[0]
        sepsis_remote = built[
            (built["cofactor"] == "sepsis")
            & (built["window"] == "91_to_365_days_before")
        ].iloc[0]
        sepsis_any = built[
            (built["cofactor"] == "sepsis")
            & (built["window"] == "any_prior")
        ].iloc[0]

        self.assertEqual(int(renal_near["n_cases_in_window"]), 1)
        self.assertEqual(int(sepsis_remote["n_cases_in_window"]), 1)
        self.assertEqual(float(sepsis_any["median_days_before_index"]), 104.0)

    def test_model_split_keeps_matched_groups_together(self) -> None:
        matched = pd.DataFrame(
            [
                {"person_id": "case1", "analysis_case": 1, "match_group_id": "match-1"},
                {"person_id": "ctrl1", "analysis_case": 0, "match_group_id": "match-1"},
                {"person_id": "case2", "analysis_case": 1, "match_group_id": "match-2"},
                {"person_id": "ctrl2", "analysis_case": 0, "match_group_id": "match-2"},
                {"person_id": "case3", "analysis_case": 1, "match_group_id": "match-3"},
                {"person_id": "ctrl3", "analysis_case": 0, "match_group_id": "match-3"},
            ]
        )

        split = add_model_splits(matched, train_fraction=0.67, cv_folds=2)

        self.assertEqual(set(split["analysis_split"]), {"train", "test"})
        self.assertEqual(split.groupby("match_group_id")["analysis_split"].nunique().max(), 1)
        self.assertTrue(split[split["analysis_split"] == "test"]["cv_fold"].isna().all())
        self.assertFalse(split[split["analysis_split"] == "train"]["cv_fold"].isna().any())

    def test_model_eligibility_flags_exclude_traumatic_and_periindex_sepsis_sets(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        matched = pd.DataFrame(
            [
                {
                    "person_id": "case-trauma",
                    "analysis_case": 1,
                    "match_group_id": "match-trauma",
                    "preindex_crush_injury": 1,
                    "periindex_crush_injury": 0,
                    "periindex_sepsis": 0,
                    "age_at_index": 50,
                    "observation_days": 1000,
                    "omop_condition_record_dates": 3,
                },
                {
                    "person_id": "ctrl-trauma",
                    "analysis_case": 0,
                    "match_group_id": "match-trauma",
                    "preindex_crush_injury": 0,
                    "periindex_crush_injury": 0,
                    "periindex_sepsis": 0,
                    "age_at_index": 50,
                    "observation_days": 1000,
                    "omop_condition_record_dates": 3,
                },
                {
                    "person_id": "case-sepsis",
                    "analysis_case": 1,
                    "match_group_id": "match-sepsis",
                    "preindex_crush_injury": 0,
                    "periindex_crush_injury": 0,
                    "periindex_sepsis": 1,
                    "age_at_index": 51,
                    "observation_days": 1000,
                    "omop_condition_record_dates": 3,
                },
                {
                    "person_id": "ctrl-sepsis",
                    "analysis_case": 0,
                    "match_group_id": "match-sepsis",
                    "preindex_crush_injury": 0,
                    "periindex_crush_injury": 0,
                    "periindex_sepsis": 0,
                    "age_at_index": 51,
                    "observation_days": 1000,
                    "omop_condition_record_dates": 3,
                },
            ]
        )

        flagged = add_model_eligibility_flags(config, matched)

        trauma_group = flagged[flagged["match_group_id"] == "match-trauma"]
        sepsis_group = flagged[flagged["match_group_id"] == "match-sepsis"]
        self.assertTrue((trauma_group["primary_model_eligible"] == 0).all())
        self.assertTrue((sepsis_group["primary_model_eligible"] == 1).all())
        self.assertTrue((sepsis_group["sensitivity_no_periindex_sepsis_eligible"] == 0).all())
        self.assertTrue((flagged["sensitivity_include_trauma_eligible"] == 1).all())

    def test_matched_table1_reports_sex_categories_including_missing(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        matched = pd.DataFrame(
            [
                {"person_id": "1", "analysis_case": 1, "age_at_index": 50, "sex_category": "female"},
                {"person_id": "2", "analysis_case": 1, "age_at_index": 51, "sex_category": "other_or_unknown"},
                {"person_id": "3", "analysis_case": 0, "age_at_index": 50, "sex_category": "male"},
                {"person_id": "4", "analysis_case": 0, "age_at_index": 51, "sex_category": "missing"},
            ]
        )

        table1 = build_matched_table1(config, matched)

        self.assertNotIn("Female sex, n (%)", set(table1["variable"]))
        self.assertIn("Sex category distribution, overall", set(table1["variable"]))
        self.assertIn("Sex: other/unknown, n (%)", set(table1["variable"]))
        missing_row = table1[table1["variable"] == "Sex: missing, n (%)"].iloc[0]
        self.assertEqual(missing_row["matched_cases"], "0 (0.0%)")
        self.assertEqual(missing_row["matched_controls"], "1 (50.0%)")


if __name__ == "__main__":
    unittest.main()
