from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.config import load_project_config
from aou_workbench.matching import match_case_controls
from aou_workbench.stage1_prepare import stage1_sample_manifest_path
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

    def test_matching_can_use_definite_as_primary_case_tier(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        config = replace(
            config,
            cohort=replace(config.cohort, primary_case_tier="definite", sensitivity_case_tiers=("broad",)),
        )
        cohort_df = build_rhabdo_cohort(config)
        matched_df = match_case_controls(cohort_df, config)
        case_rows = matched_df[matched_df["analysis_case"] == 1]

        self.assertGreater(int(case_rows.shape[0]), 0)
        self.assertTrue(case_rows["case_tier"].eq("definite").all())
        self.assertTrue(case_rows["definite_rhabdo_case"].fillna(0).astype(int).eq(1).all())

    def test_matching_accepts_reloaded_tsv_with_string_dates(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        cohort_df = build_rhabdo_cohort(config)
        reload_path = Path(paths["output_dir"]) / "cohort_reload.tsv"
        reload_path.parent.mkdir(parents=True, exist_ok=True)
        cohort_df.to_csv(reload_path, sep="\t", index=False)
        reloaded = pd.read_csv(reload_path, sep="\t")

        matched_df = match_case_controls(reloaded, config)

        self.assertGreater(int((matched_df["analysis_case"] == 1).sum()), 0)
        self.assertGreater(int((matched_df["analysis_case"] == 0).sum()), 0)

    def test_matching_respects_age_sex_and_ancestry_criteria(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        config = replace(
            config,
            cohort=replace(
                config.cohort,
                exact_match_columns=("is_female", "ancestry_pred"),
                age_tolerance_years=8,
            ),
        )
        cohort_df = build_rhabdo_cohort(config)
        matched_df = match_case_controls(cohort_df, config)

        cases = matched_df[matched_df["analysis_case"] == 1].set_index("match_group_id")
        controls = matched_df[matched_df["analysis_case"] == 0]
        self.assertFalse(controls.empty)
        for control in controls.itertuples(index=False):
            case = cases.loc[control.match_group_id]
            self.assertEqual(control.is_female, case.is_female)
            self.assertEqual(control.ancestry_pred, case.ancestry_pred)
            self.assertLessEqual(abs(float(control.age_at_index) - float(case.age_at_index)), 8.0)

    def test_matching_applies_wgs_and_max_unrelated_restriction_before_matching(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        pd.DataFrame({"person_id": [str(person_id) for person_id in range(1, 13)]}).to_csv(
            stage1_sample_manifest_path(paths["stage1_table"]),
            sep="\t",
            index=False,
        )
        config = replace(
            config,
            workbench=replace(config.workbench, max_unrelated_path=paths["max_unrelated"]),
        )
        cohort_df = build_rhabdo_cohort(config)

        matched_df = match_case_controls(cohort_df, config)

        allowed = {str(person_id) for person_id in range(1, 11)}
        self.assertTrue(set(matched_df["person_id"].astype(str)).issubset(allowed))
