from __future__ import annotations

import unittest

import pandas as pd

from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.config import load_project_config
from aou_workbench.paths import build_output_paths
from aou_workbench.stage1_prepare import stage1_sample_manifest_path
from aou_workbench.stage1_prior_variants import (
    _restrict_to_stage1_wgs_samples,
    _stage1_comparisons,
    _variant_exposure,
    run_stage1_prior_variants,
)
from tests.support import build_demo_project_tree


class Stage1PriorVariantTests(unittest.TestCase):
    def test_stage1_comparisons_use_non_rhabdo_controls_for_both_levels(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        cohort_df = build_rhabdo_cohort(config)

        comparisons = _stage1_comparisons(cohort_df, config)

        comparison_names = [item["comparison"] for item in comparisons]
        self.assertIn("omop_rhabdo_vs_non_rhabdo", comparison_names)
        self.assertIn("omop_rhabdo_plus_ck_vs_non_rhabdo", comparison_names)
        omop_any = comparisons[0]["sample_df"]
        omop_ck = comparisons[1]["sample_df"]
        self.assertEqual(int((omop_any["rhabdo_case"] == 1).sum()), 4)
        self.assertEqual(int((omop_ck["rhabdo_primary_case"] == 1).sum()), 3)
        self.assertEqual(int((omop_ck["rhabdo_case"] == 1).sum()), 3)
        self.assertIn("omop_rhabdo_vs_non_rhabdo_ancestry_afr", {item["comparison"] for item in comparisons})
        self.assertIn("omop_rhabdo_vs_non_rhabdo_ancestry_eur", {item["comparison"] for item in comparisons})

    def test_run_stage1_prior_variants_writes_two_comparisons_from_built_cohort(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        cohort_df = build_rhabdo_cohort(config)
        output_paths = build_output_paths(config)
        manifest_path = stage1_sample_manifest_path(config.analysis.stage1.variant_table)
        pd.DataFrame({"person_id": cohort_df["person_id"].astype(str)}).to_csv(manifest_path, sep="\t", index=False)

        result = run_stage1_prior_variants(config, cohort_df, output_paths)

        self.assertIn("omop_rhabdo_vs_non_rhabdo", set(result["comparison"]))
        self.assertIn("omop_rhabdo_plus_ck_vs_non_rhabdo", set(result["comparison"]))
        hbb = result[result["label"] == "HBB sickle trait"].set_index("comparison")
        self.assertEqual(int(hbb.loc["omop_rhabdo_vs_non_rhabdo", "case_carriers"]), 2)
        self.assertEqual(int(hbb.loc["omop_rhabdo_vs_non_rhabdo", "control_carriers"]), 1)
        self.assertEqual(int(hbb.loc["omop_rhabdo_plus_ck_vs_non_rhabdo", "case_carriers"]), 1)
        self.assertEqual(int(hbb.loc["omop_rhabdo_plus_ck_vs_non_rhabdo", "control_carriers"]), 1)

    def test_stage1_prior_variants_restricts_denominators_to_wgs_manifest(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        cohort_df = build_rhabdo_cohort(config)
        manifest_path = stage1_sample_manifest_path(config.analysis.stage1.variant_table)
        pd.DataFrame({"person_id": ["1", "2", "3", "4", "5", "9"]}).to_csv(manifest_path, sep="\t", index=False)

        restricted = _restrict_to_stage1_wgs_samples(cohort_df, config)

        self.assertEqual(set(restricted["person_id"].astype(str)), {"1", "2", "3", "4", "5", "9"})
        self.assertEqual(int((restricted["rhabdo_case"] == 1).sum()), 4)
        self.assertEqual(int((restricted["rhabdo_case"] == 0).sum()), 2)

    def test_stage1_prior_variants_requires_wgs_manifest(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        cohort_df = build_rhabdo_cohort(config)

        with self.assertRaises(RuntimeError):
            _restrict_to_stage1_wgs_samples(cohort_df, config)

    def test_variant_exposure_uses_homozygous_alt_model_when_requested(self) -> None:
        subset = pd.DataFrame(
            [
                {"person_id": "1", "dosage": 1.0},
                {"person_id": "2", "dosage": 2.0},
                {"person_id": "3", "dosage": 0.0},
            ]
        )

        exposure, threshold = _variant_exposure(subset, exact_test_model="hom_alt_vs_rest")

        self.assertEqual(threshold, 1.0)
        self.assertEqual(exposure.to_dict(), {"1": 0.0, "2": 1.0, "3": 0.0})

    def test_variant_exposure_preserves_dosage_for_carrier_model(self) -> None:
        subset = pd.DataFrame(
            [
                {"person_id": "1", "dosage": 1.0},
                {"person_id": "1", "dosage": 2.0},
                {"person_id": "2", "dosage": 1.0},
            ]
        )

        exposure, threshold = _variant_exposure(subset, exact_test_model="carrier_vs_noncarrier")

        self.assertEqual(threshold, 1.0)
        self.assertEqual(exposure.to_dict(), {"1": 2.0, "2": 1.0})


if __name__ == "__main__":
    unittest.main()
