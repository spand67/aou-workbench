from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.cli import _build_parser, _load_hail_pilot_matched_input
from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.cohort_summary import characterize_case_control_cohort, clinical_model_input_path
from aou_workbench.config import load_project_config
from aou_workbench.matching import match_case_controls
from aou_workbench.paths import build_output_paths
from aou_workbench.pipeline import build_cohort_artifacts, match_controls_artifacts
from aou_workbench.stage1_prepare import stage1_sample_manifest_path
from aou_workbench.stage4_hail_gwas import (
    PASSING_GENOTYPE_FILTER_VALUES,
    _hail_pilot_sample_frame,
    _hail_sample_frame,
    _matrix_table_gt,
    _normalize_autosomal_chromosomes,
    _normalize_chromosomes,
    _pilot_default_target_partitions,
    _pilot_restricts_to_wgs_manifest,
    _variant_qc_summary_rows,
    hail_pilot_qc_pass_mt_uri,
    hail_pilot_default_label,
    hail_pilot_results_ht_uri,
    hail_pilot_results_path,
)
from tests.support import build_demo_project_tree


class Stage4HailGwasTests(unittest.TestCase):
    def test_matrix_table_gt_masks_filtered_gt_entries(self) -> None:
        class FakeExpression:
            def __init__(self, text: str):
                self.text = text

            def __or__(self, other):
                return FakeExpression(f"({self.text}) OR ({other.text})")

        class FakeLiteral:
            def __init__(self, values):
                self.values = values

            def contains(self, value):
                return FakeExpression(f"{value.text} IN {sorted(self.values)}")

        class FakeHail:
            def __init__(self):
                self.condition = None

            def literal(self, values):
                return FakeLiteral(values)

            def is_missing(self, value):
                return FakeExpression(f"is_missing({value.text})")

            def or_missing(self, condition, value):
                self.condition = condition
                return ("or_missing", condition, value)

        class FakeMatrixTable:
            entry = {"GT": object(), "FT": object()}
            GT = FakeExpression("GT")
            FT = FakeExpression("FT")

        fake_hail = FakeHail()
        result = _matrix_table_gt(FakeMatrixTable(), fake_hail)

        self.assertEqual(result[0], "or_missing")
        self.assertIs(result[2], FakeMatrixTable.GT)
        self.assertIn("is_missing(FT)", fake_hail.condition.text)
        self.assertIn("PASS", fake_hail.condition.text)
        self.assertIn(".", PASSING_GENOTYPE_FILTER_VALUES)

    def test_normalize_chromosomes_defaults_and_deduplicates(self) -> None:
        self.assertEqual(_normalize_chromosomes(None), [str(chrom) for chrom in range(1, 23)])
        self.assertEqual(_normalize_chromosomes(["chr1", "1", "19", "chr19", "22"]), ["1", "19", "22"])

    def test_normalize_autosomal_chromosomes_rejects_non_autosomes(self) -> None:
        self.assertEqual(_normalize_autosomal_chromosomes(["chr22", "1"]), ["22", "1"])
        with self.assertRaises(ValueError):
            _normalize_autosomal_chromosomes(["22", "X"])

    def test_hail_pilot_genotype_source_defaults(self) -> None:
        self.assertEqual(hail_pilot_default_label("acaf", ["22"]), "acaf_chr22_maf05_train_qc")
        self.assertEqual(
            hail_pilot_default_label("microarray", ["22"], min_maf=0.05),
            "microarray_chr22_maf05_train_qc",
        )
        self.assertTrue(_pilot_restricts_to_wgs_manifest("acaf"))
        self.assertFalse(_pilot_restricts_to_wgs_manifest("microarray"))
        self.assertEqual(_pilot_default_target_partitions("acaf", ["22"]), 0)
        self.assertEqual(_pilot_default_target_partitions("microarray", ["22"]), 64)
        self.assertEqual(_pilot_default_target_partitions("microarray", ["21", "22"]), 128)

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

        bad_case = matched_df[matched_df["analysis_case"] == 1].iloc[0].copy()
        bad_case["person_id"] = "bad-broad-case"
        bad_case["broad_rhabdo_case"] = 0
        bad_control = matched_df[matched_df["analysis_case"] == 0].iloc[0].copy()
        bad_control["person_id"] = "bad-control"
        bad_control["case_tier"] = "indeterminate_ck_only"
        bad_control["eligible_control"] = 0
        matched_df = pd.concat([matched_df, pd.DataFrame([bad_case, bad_control])], ignore_index=True)

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
        self.assertNotIn("bad-broad-case", observed)
        self.assertNotIn("bad-control", observed)
        self.assertIn("pc1", hail_covariates)
        self.assertEqual(raw_covariates, ["age_at_index", "is_female", "pc1", "pc2", "pc3", "pc4", "pc5"])
        self.assertNotIn("observation_days", raw_covariates)
        self.assertNotIn("omop_condition_record_dates", raw_covariates)
        self.assertNotIn("sepsis", raw_covariates)
        self.assertNotIn("renal_injury", raw_covariates)
        self.assertNotIn("crush_injury", raw_covariates)
        self.assertIsInstance(dropped_covariates, list)
        self.assertTrue(wgs_manifest_used)
        self.assertLess(counts["after_case_control_definition_participants"], counts["after_eligibility_participants"])
        self.assertLess(counts["after_complete_case_participants"], counts["matched_input_participants"])
        self.assertGreater(counts["after_complete_case_cases"], 0)
        self.assertGreater(counts["after_complete_case_controls"], 0)

    def test_hail_pilot_sample_frame_can_skip_wgs_manifest_for_microarray(self) -> None:
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
        pd.DataFrame({"person_id": [str(matched_df["person_id"].iloc[0])] }).to_csv(
            stage1_sample_manifest_path(paths["stage1_table"]),
            sep="\t",
            index=False,
        )

        _, _, _, _, restricted_counts, restricted_manifest_used = _hail_pilot_sample_frame(
            matched_df,
            config,
            restrict_to_wgs_manifest=True,
        )
        _, _, _, _, unrestricted_counts, unrestricted_manifest_used = _hail_pilot_sample_frame(
            matched_df,
            config,
            restrict_to_wgs_manifest=False,
        )

        self.assertTrue(restricted_manifest_used)
        self.assertFalse(unrestricted_manifest_used)
        self.assertLess(
            restricted_counts["after_wgs_manifest_participants"],
            unrestricted_counts["after_wgs_manifest_participants"],
        )
        self.assertEqual(
            unrestricted_counts["after_wgs_manifest_participants"],
            unrestricted_counts["after_case_control_definition_participants"],
        )

    def test_hail_pilot_loader_refreshes_stale_clinical_model_input(self) -> None:
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
        characterize_case_control_cohort(effective, cohort_df, matched_df, output_paths)

        input_path = Path(clinical_model_input_path(output_paths))
        stale = pd.read_csv(input_path, sep="\t").drop(
            columns=[
                "eligible_control",
                "eligible_ehr_denominator",
                "broad_rhabdo_case",
                "definite_rhabdo_case",
                "high_ck_without_rhabdo",
            ],
            errors="ignore",
        )
        stale.to_csv(input_path, sep="\t", index=False)

        _, _, refreshed = _load_hail_pilot_matched_input(config, eligibility_flag="primary_model_eligible")

        for column in (
            "eligible_control",
            "eligible_ehr_denominator",
            "broad_rhabdo_case",
            "definite_rhabdo_case",
            "high_ck_without_rhabdo",
        ):
            self.assertIn(column, refreshed.columns)

    def test_variant_qc_summary_rows_preserve_sequential_counts(self) -> None:
        rows = _variant_qc_summary_rows(
            chromosomes=["22"],
            initial_rows=100,
            biallelic_rows=90,
            maf_rows=80,
            mac_rows=75,
            call_rate_rows=70,
            hwe_rows=60,
            min_maf=0.05,
            min_mac=20,
            min_call_rate=0.98,
            hwe_p_control=1e-6,
        )
        summary = pd.DataFrame(rows)

        self.assertEqual(
            summary["filter"].tolist(),
            [
                "interval_and_sample_restriction",
                "autosomal_biallelic_snp",
                "maf",
                "minor_allele_count",
                "call_rate",
                "control_hwe_p",
            ],
        )
        self.assertEqual(summary["rows_before"].tolist(), [100, 100, 90, 80, 75, 70])
        self.assertEqual(summary["rows_after"].tolist(), [100, 90, 80, 75, 70, 60])
        self.assertEqual(summary["rows_removed"].tolist(), [0, 10, 10, 5, 5, 10])

    def test_variant_qc_summary_rows_can_report_hwe_without_filtering(self) -> None:
        rows = _variant_qc_summary_rows(
            chromosomes=["19"],
            initial_rows=100,
            biallelic_rows=90,
            maf_rows=80,
            mac_rows=75,
            call_rate_rows=70,
            hwe_rows=50,
            final_rows=70,
            min_maf=0.01,
            min_mac=20,
            min_call_rate=0.98,
            hwe_p_control=1e-6,
            hwe_filter_mode="report-only",
        )
        summary = pd.DataFrame(rows)

        self.assertEqual(summary["filter"].tolist()[-2:], ["control_hwe_p_report_only", "final_qc"])
        self.assertEqual(summary["rows_after"].tolist()[-2:], [50, 70])
        self.assertEqual(summary["rows_removed"].tolist()[-2:], [20, 0])

    def test_variant_qc_summary_rows_allow_skipped_denominator_counts(self) -> None:
        rows = _variant_qc_summary_rows(
            chromosomes=["19"],
            initial_rows=None,
            biallelic_rows=None,
            maf_rows=1000,
            mac_rows=900,
            call_rate_rows=850,
            hwe_rows=800,
            min_maf=0.01,
            min_mac=20,
            min_call_rate=0.98,
            hwe_p_control=1e-6,
        )
        summary = pd.DataFrame(rows)

        self.assertTrue(summary.loc[summary["filter"].eq("interval_and_sample_restriction"), "rows_after"].isna().all())
        self.assertTrue(summary.loc[summary["filter"].eq("autosomal_biallelic_snp"), "rows_after"].isna().all())
        self.assertEqual(
            int(summary.loc[summary["filter"].eq("minor_allele_count"), "rows_removed"].iloc[0]),
            100,
        )

    def test_variant_qc_summary_rows_allow_skipped_qc_counts(self) -> None:
        rows = _variant_qc_summary_rows(
            chromosomes=["5"],
            initial_rows=None,
            biallelic_rows=None,
            maf_rows=None,
            mac_rows=None,
            call_rate_rows=None,
            hwe_rows=None,
            final_rows=None,
            min_maf=0.01,
            min_mac=20,
            min_call_rate=0.98,
            hwe_p_control=1e-6,
            hwe_filter_mode="report-only",
        )
        summary = pd.DataFrame(rows)

        self.assertEqual(summary["filter"].tolist()[-2:], ["control_hwe_p_report_only", "final_qc"])
        self.assertTrue(summary[["rows_before", "rows_after", "rows_removed"]].isna().all().all())

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

    def test_hail_pilot_gcs_paths_use_workspace_bucket(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        config = replace(config, workbench=replace(config.workbench, workspace_bucket="gs://example-bucket"))

        self.assertEqual(
            hail_pilot_results_ht_uri(config, "acaf_chr19_ryr1_maf01"),
            "gs://example-bucket/aou-workbench/unit-test-rhabdo/stage4/hail_pilot/acaf-chr19-ryr1-maf01/gwas_results.ht",
        )
        self.assertEqual(
            hail_pilot_qc_pass_mt_uri(config, "acaf_chr19_ryr1_maf01"),
            "gs://example-bucket/aou-workbench/unit-test-rhabdo/stage4/hail_pilot/acaf-chr19-ryr1-maf01/qc_pass_genotypes.mt",
        )

    def test_hail_pilot_cli_defaults_match_guidance_qc(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["run-hail-pilot-gwas"])

        self.assertEqual(args.command, "run-hail-pilot-gwas")
        self.assertEqual(args.chromosomes, "22")
        self.assertEqual(args.min_maf, 0.05)
        self.assertEqual(args.min_mac, 20)
        self.assertEqual(args.min_call_rate, 0.98)
        self.assertEqual(args.hwe_p_control, 1e-6)
        self.assertEqual(args.hwe_filter_mode, "filter")
        self.assertEqual(args.analysis_split, "train")
        self.assertEqual(args.eligibility_flag, "primary_model_eligible")
        self.assertEqual(args.genotype_source, "acaf")
        self.assertIsNone(args.label)
        self.assertIsNone(args.target_partitions)
        self.assertFalse(args.write_qc_mt)
        self.assertFalse(args.export_hail_results_tsv)
        self.assertFalse(args.skip_variant_row_counts)
        self.assertFalse(args.skip_qc_counts)
        self.assertEqual(args.results_preview_n, 100000)

        microarray_args = parser.parse_args(
            ["run-hail-pilot-gwas", "--genotype-source", "microarray", "--target-partitions", "96"]
        )
        self.assertEqual(microarray_args.genotype_source, "microarray")
        self.assertEqual(microarray_args.target_partitions, 96)

        hwe_args = parser.parse_args(
            [
                "run-hail-pilot-gwas",
                "--hwe-filter-mode",
                "report-only",
                "--write-qc-mt",
                "--export-hail-results-tsv",
                "--skip-variant-row-counts",
                "--skip-qc-counts",
                "--results-preview-n",
                "250000",
            ]
        )
        self.assertEqual(hwe_args.hwe_filter_mode, "report-only")
        self.assertTrue(hwe_args.write_qc_mt)
        self.assertTrue(hwe_args.export_hail_results_tsv)
        self.assertTrue(hwe_args.skip_variant_row_counts)
        self.assertTrue(hwe_args.skip_qc_counts)
        self.assertEqual(hwe_args.results_preview_n, 250000)

        no_preview_args = parser.parse_args(["run-hail-pilot-gwas", "--results-preview-n", "0"])
        self.assertEqual(no_preview_args.results_preview_n, 0)


if __name__ == "__main__":
    unittest.main()
