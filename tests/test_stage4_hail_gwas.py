from __future__ import annotations

from dataclasses import replace
import unittest

import pandas as pd

from aou_workbench.cli import _build_parser
from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.config import load_project_config
from aou_workbench.matching import match_case_controls
from aou_workbench.paths import build_output_paths
from aou_workbench.stage1_prepare import stage1_sample_manifest_path
from aou_workbench.stage4_hail_gwas import (
    _hail_pilot_sample_frame,
    _hail_sample_frame,
    _normalize_autosomal_chromosomes,
    _normalize_chromosomes,
    hail_pilot_results_path,
)
from tests.support import build_demo_project_tree


class Stage4HailGwasTests(unittest.TestCase):
    def test_normalize_chromosomes_defaults_and_deduplicates(self) -> None:
        self.assertEqual(_normalize_chromosomes(None), [str(chrom) for chrom in range(1, 23)])
        self.assertEqual(_normalize_chromosomes(["chr1", "1", "19", "chr19", "22"]), ["1", "19", "22"])

    def test_normalize_autosomal_chromosomes_rejects_non_autosomes(self) -> None:
        self.assertEqual(_normalize_autosomal_chromosomes(["chr22", "1"]), ["22", "1"])
        with self.assertRaises(ValueError):
            _normalize_autosomal_chromosomes(["22", "X"])

    def test_hail_sample_frame_builds_complete_case_numeric_covariates(self) -> None:
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
        matched_df = match_case_controls(cohort_df, config).copy()
        matched_df.loc[matched_df.index[0], "pc1"] = None
        matched_df.loc[matched_df.index[1], "ancestry_pred"] = ""

        sample_df, hail_covariates, raw_covariates, dropped_covariates = _hail_sample_frame(matched_df, config)

        self.assertLess(len(sample_df), matched_df["person_id"].astype(str).nunique())
        self.assertIn("age_at_index", hail_covariates)
        self.assertIn("pc1", hail_covariates)
        self.assertIn("ancestry_pred", raw_covariates)
        self.assertNotIn("ancestry_pred", hail_covariates)
        self.assertNotIn("ancestry_pred", dropped_covariates)
        self.assertFalse(sample_df[hail_covariates + [config.analysis.matched_outcome_column]].isna().any().any())
        self.assertTrue(set(sample_df["person_id"].astype(str)).issubset({str(person_id) for person_id in range(1, 11)}))

    def test_hail_pilot_sample_frame_filters_train_eligible_complete_cases(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        cohort_df = build_rhabdo_cohort(config)
        matched_df = match_case_controls(cohort_df, config).copy()
        matched_df["analysis_split"] = "train"
        matched_df["primary_model_eligible"] = 1

        test_person = str(matched_df.loc[matched_df.index[0], "person_id"])
        ineligible_person = str(matched_df.loc[matched_df.index[1], "person_id"])
        incomplete_person = str(matched_df.loc[matched_df.index[2], "person_id"])
        matched_df.loc[matched_df.index[0], "analysis_split"] = "test"
        matched_df.loc[matched_df.index[1], "primary_model_eligible"] = 0
        matched_df.loc[matched_df.index[2], "pc1"] = None

        sample_df, hail_covariates, raw_covariates, dropped_covariates, counts, wgs_manifest_used = (
            _hail_pilot_sample_frame(matched_df, config)
        )
        observed = set(sample_df["person_id"].astype(str))

        self.assertNotIn(test_person, observed)
        self.assertNotIn(ineligible_person, observed)
        self.assertNotIn(incomplete_person, observed)
        self.assertIn("pc1", hail_covariates)
        self.assertEqual(raw_covariates, ["age_at_index", "is_female", "pc1", "pc2", "pc3", "pc4", "pc5"])
        self.assertIsInstance(dropped_covariates, list)
        self.assertTrue(wgs_manifest_used)
        self.assertLess(counts["after_complete_case_participants"], counts["matched_input_participants"])

    def test_hail_pilot_paths_are_labeled_and_do_not_overwrite_stage4(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        output_paths = build_output_paths(config)
        pilot_path = hail_pilot_results_path(output_paths, "acaf_chr22_maf05_train_qc")

        self.assertIn("/stage4/hail_pilot/acaf-chr22-maf05-train-qc/gwas_results.tsv", pilot_path)
        self.assertNotEqual(pilot_path, output_paths.stage4_full_results_tsv)

    def test_hail_pilot_cli_defaults_match_guidance_qc(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["run-hail-pilot-gwas"])

        self.assertEqual(args.command, "run-hail-pilot-gwas")
        self.assertEqual(args.chromosomes, "22")
        self.assertEqual(args.min_maf, 0.05)
        self.assertEqual(args.min_mac, 20)
        self.assertEqual(args.min_call_rate, 0.98)
        self.assertEqual(args.hwe_p_control, 1e-6)
        self.assertEqual(args.analysis_split, "train")
        self.assertEqual(args.eligibility_flag, "primary_model_eligible")


if __name__ == "__main__":
    unittest.main()
