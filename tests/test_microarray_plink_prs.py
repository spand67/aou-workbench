from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.cli import _build_parser
from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.config import load_project_config
from aou_workbench.matching import match_case_controls
from aou_workbench.microarray_plink_prs import (
    _valid_gwas_weights,
    _write_clumped_weights,
    _write_prs_sample_keep,
    build_microarray_prs_commands,
    microarray_prs_default_label,
    microarray_prs_keep_path,
    microarray_prs_output_dir,
    microarray_prs_range_path,
    parse_prs_score_files,
    parse_thresholds,
    prs_metrics,
)
from aou_workbench.paths import build_output_paths
from tests.support import build_demo_project_tree


class MicroarrayPlinkPrsTests(unittest.TestCase):
    def _config_and_matched(self):
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
        test_groups = matched_df["match_group_id"].drop_duplicates().astype(str).head(2)
        matched_df.loc[matched_df["match_group_id"].astype(str).isin(set(test_groups)), "analysis_split"] = "test"
        return config, matched_df

    def test_prs_sample_keep_scores_test_split_only(self) -> None:
        config, matched_df = self._config_and_matched()
        output_paths = build_output_paths(config)
        gwas_label = "microarray_plink_autosomes_maf05_train_qc"
        prs_label = "test_clumped_threshold_grid"
        fam_df = pd.DataFrame({"IID": matched_df["person_id"].astype(str).tolist()})

        sample, counts = _write_prs_sample_keep(
            config,
            matched_df,
            fam_df,
            output_paths,
            gwas_label=gwas_label,
            prs_label=prs_label,
            score_split="test",
            eligibility_flag="primary_model_eligible",
        )

        expected = matched_df[
            (matched_df["analysis_split"] == "test")
            & (pd.to_numeric(matched_df["primary_model_eligible"], errors="coerce") == 1)
        ]
        self.assertFalse(sample.empty)
        self.assertEqual(set(sample["person_id"].astype(str)), set(expected["person_id"].astype(str)))
        self.assertEqual(counts["after_microarray_fam_overlap_participants"], expected["person_id"].nunique())
        self.assertEqual(counts["after_microarray_fam_overlap_cases"], int((expected["analysis_case"] == 1).sum()))
        self.assertEqual(counts["after_microarray_fam_overlap_controls"], int((expected["analysis_case"] == 0).sum()))
        keep = pd.read_csv(
            microarray_prs_keep_path(output_paths, gwas_label, prs_label, "test"),
            sep=r"\s+",
            header=None,
            names=["FID", "IID"],
            dtype=str,
        )
        self.assertEqual(set(keep["IID"]), set(expected["person_id"].astype(str)))

    def test_gwas_weights_drop_malformed_and_duplicate_rows(self) -> None:
        gwas = pd.DataFrame(
            [
                {"variant_id": "v1", "effect_allele": "A", "beta": 0.1, "regression_p": 0.01},
                {"variant_id": "v2", "effect_allele": "G", "beta": 0.2, "regression_p": 0.02},
                {"variant_id": "v2", "effect_allele": "G", "beta": 0.3, "regression_p": 0.03},
                {"variant_id": "v3", "effect_allele": "", "beta": 0.4, "regression_p": 0.04},
                {"variant_id": "v4", "effect_allele": "T", "beta": None, "regression_p": 0.05},
                {"variant_id": "v5", "effect_allele": "C", "beta": 0.5, "regression_p": None},
            ]
        )

        weights = _valid_gwas_weights(gwas)

        self.assertEqual(weights["ID"].tolist(), ["v1"])
        self.assertEqual(weights[["ID", "A1", "BETA", "P"]].columns.tolist(), ["ID", "A1", "BETA", "P"])

    def test_clumped_weights_and_ranges_are_written(self) -> None:
        config, _ = self._config_and_matched()
        output_paths = build_output_paths(config)
        gwas_label = "gwas"
        prs_label = "prs"
        weights = pd.DataFrame(
            [
                {"ID": "v1", "A1": "A", "BETA": 0.1, "P": 0.01},
                {"ID": "v2", "A1": "G", "BETA": 0.2, "P": 0.02},
            ]
        )

        selected = _write_clumped_weights(weights, {"v2"}, str(Path(microarray_prs_output_dir(output_paths, gwas_label, prs_label)) / "weights.tsv"))

        self.assertEqual(selected["ID"].tolist(), ["v2"])
        self.assertEqual(parse_thresholds("5e-8,0.01,1"), [5e-08, 0.01, 1.0])

    def test_prs_commands_build_clump_and_qscore_range(self) -> None:
        config, _ = self._config_and_matched()
        output_paths = build_output_paths(config)
        commands = build_microarray_prs_commands(
            plink2_bin="plink2",
            plink_prefix="/tmp/arrays",
            paths=output_paths,
            gwas_label="gwas",
            prs_label="prs",
            clump_split="train",
            score_split="test",
            clump_input_path="/tmp/clump.tsv",
            clump_r2=0.1,
            clump_kb=250,
            clump_p1=1.0,
            clump_p2=1.0,
            threads=8,
            memory_mb=16000,
        )
        clump_cmd, score_cmd = commands

        self.assertIn("--clump", clump_cmd)
        self.assertIn("--clump-id-field", clump_cmd)
        self.assertIn("ID", clump_cmd)
        self.assertIn("train_keep.tsv", " ".join(clump_cmd))
        self.assertIn("--clump-r2", clump_cmd)
        self.assertIn("0.1", clump_cmd)
        self.assertIn("--score", score_cmd)
        self.assertIn("--extract", score_cmd)
        self.assertIn("prs_variant_ids.txt", " ".join(score_cmd))
        self.assertIn("--q-score-range", score_cmd)
        self.assertIn("cols=+scoresums", score_cmd)
        self.assertIn("list-variants", score_cmd)
        self.assertIn("--threads", score_cmd)
        self.assertIn("8", score_cmd)

    def test_parse_prs_score_files_combines_threshold_outputs(self) -> None:
        config, _ = self._config_and_matched()
        output_paths = build_output_paths(config)
        outdir = Path(microarray_prs_output_dir(output_paths, "gwas", "prs"))
        outdir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            [
                {"#FID": "0", "IID": "1", "SCORE1_SUM": "0.2"},
                {"#FID": "0", "IID": "2", "SCORE1_SUM": "-0.1"},
            ]
        ).to_csv(outdir / "plink_prs.p0_01.sscore", sep="\t", index=False)
        pd.DataFrame(
            [
                {"#FID": "0", "IID": "1", "SCORE1_SUM": "0.3"},
                {"#FID": "0", "IID": "2", "SCORE1_SUM": "-0.2"},
            ]
        ).to_csv(outdir / "plink_prs.p1.sscore", sep="\t", index=False)

        scores = parse_prs_score_files(output_paths, "gwas", "prs")

        self.assertEqual(scores.shape[0], 4)
        self.assertEqual(set(scores["threshold_label"]), {"p0_01", "p1"})
        self.assertEqual(set(scores["person_id"]), {"1", "2"})

    def test_parse_prs_score_files_accepts_plink_without_score_named_column(self) -> None:
        config, _ = self._config_and_matched()
        output_paths = build_output_paths(config)
        outdir = Path(microarray_prs_output_dir(output_paths, "gwas_fallback", "prs"))
        outdir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            [
                {"#FID": "0", "IID": "1", "ALLELE_CT": "10", "DENOM": "10", "PRS_TOTAL": "0.2"},
                {"#FID": "0", "IID": "2", "ALLELE_CT": "10", "DENOM": "10", "PRS_TOTAL": "-0.1"},
            ]
        ).to_csv(outdir / "plink_prs.p0_01.sscore", sep="\t", index=False)

        scores = parse_prs_score_files(output_paths, "gwas_fallback", "prs")

        self.assertEqual(scores["prs_score"].tolist(), [0.2, -0.1])

    def test_prs_metrics_compute_auc_and_or_per_sd(self) -> None:
        scores = pd.DataFrame(
            {
                "person_id": ["1", "2", "3", "4"],
                "threshold_label": ["p1", "p1", "p1", "p1"],
                "prs_score": [0.8, 0.7, -0.2, -0.1],
            }
        )
        sample = pd.DataFrame({"person_id": ["1", "2", "3", "4"], "analysis_case": [1, 1, 0, 0]})
        weights = pd.DataFrame({"ID": ["v1", "v2"], "P": [0.01, 0.2]})

        metrics = prs_metrics(scores, sample, weights=weights)

        self.assertEqual(metrics.shape[0], 1)
        self.assertAlmostEqual(metrics.loc[0, "roc_auc"], 1.0)
        self.assertEqual(metrics.loc[0, "n_participants"], 4)
        self.assertEqual(metrics.loc[0, "n_variants_scored"], 2)
        self.assertGreater(metrics.loc[0, "or_per_prs_sd"], 1.0)

    def test_prs_cli_defaults(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["run-microarray-plink-prs", "--gwas-label", "microarray_plink_autosomes_maf05_train_qc"])

        self.assertEqual(args.command, "run-microarray-plink-prs")
        self.assertEqual(args.gwas_label, "microarray_plink_autosomes_maf05_train_qc")
        self.assertEqual(args.score_split, "test")
        self.assertEqual(args.clump_r2, 0.1)
        self.assertEqual(args.clump_kb, 250)
        self.assertEqual(args.clump_p1, 0.01)
        self.assertEqual(args.clump_p2, 0.01)
        self.assertEqual(microarray_prs_default_label("test"), "test_clumped_threshold_grid")
        self.assertEqual(args.thresholds, "0.01")


if __name__ == "__main__":
    unittest.main()
