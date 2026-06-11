from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.cli import _build_parser
from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.config import DEFAULT_MICROARRAY_PLINK_PREFIX, load_project_config
from aou_workbench.matching import match_case_controls
from aou_workbench.microarray_plink_gwas import (
    build_microarray_plink_commands,
    microarray_plink_controls_keep_path,
    microarray_plink_default_label,
    microarray_plink_keep_path,
    microarray_plink_output_dir,
    parse_plink_glm_results,
    write_microarray_plink_sample_files,
)
from aou_workbench.paths import build_output_paths
from tests.support import build_demo_project_tree


class MicroarrayPlinkGwasTests(unittest.TestCase):
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
        return config, matched_df

    def test_default_config_includes_microarray_plink_prefix(self) -> None:
        config, _ = self._config_and_matched()

        self.assertEqual(config.workbench.microarray_plink_prefix, DEFAULT_MICROARRAY_PLINK_PREFIX)

    def test_microarray_plink_sample_files_intersect_fam_and_keep_controls(self) -> None:
        config, matched_df = self._config_and_matched()
        output_paths = build_output_paths(config)
        label = "plink_test"
        kept_people = matched_df["person_id"].astype(str).head(6).tolist()
        fam_df = pd.DataFrame(
            {
                "FID": ["0"] * len(kept_people),
                "IID": kept_people,
                "father": ["0"] * len(kept_people),
                "mother": ["0"] * len(kept_people),
                "sex": ["0"] * len(kept_people),
                "phenotype": ["NA"] * len(kept_people),
            }
        )

        sample_df, covariates, raw_covariates, dropped_covariates, counts = write_microarray_plink_sample_files(
            config,
            matched_df,
            fam_df,
            output_paths,
            label=label,
        )

        keep_path = Path(microarray_plink_keep_path(output_paths, label))
        controls_keep_path = Path(microarray_plink_controls_keep_path(output_paths, label))
        keep = pd.read_csv(keep_path, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
        controls_keep = pd.read_csv(controls_keep_path, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)

        self.assertEqual(set(sample_df["person_id"].astype(str)), set(kept_people))
        self.assertTrue((keep["FID"] == "0").all())
        self.assertTrue(set(controls_keep["IID"]).issubset(set(keep["IID"])))
        self.assertEqual(covariates, ["age_at_index", "is_female", "pc1", "pc2", "pc3", "pc4", "pc5"])
        self.assertEqual(raw_covariates, covariates)
        self.assertEqual(dropped_covariates, [])
        self.assertEqual(counts["after_microarray_fam_overlap_participants"], len(kept_people))

    def test_microarray_plink_commands_use_qc_and_control_hwe(self) -> None:
        config, _ = self._config_and_matched()
        output_paths = build_output_paths(config)
        label = "plink_test"
        commands = build_microarray_plink_commands(
            plink2_bin="plink2",
            plink_prefix="/tmp/arrays",
            paths=output_paths,
            label=label,
            chromosomes=["22"],
            min_maf=0.05,
            min_mac=20,
            min_call_rate=0.98,
            hwe_p_control=1e-6,
            outcome_column="analysis_case",
            covariates=["age_at_index", "is_female", "pc1", "pc2", "pc3", "pc4", "pc5"],
            threads=8,
            memory_mb=32000,
        )
        analysis_qc_cmd, control_hwe_cmd, assoc_cmd = commands

        self.assertIn("--snps-only", analysis_qc_cmd)
        self.assertIn("--maf", analysis_qc_cmd)
        self.assertIn("0.05", analysis_qc_cmd)
        self.assertIn("--mac", analysis_qc_cmd)
        self.assertIn("20", analysis_qc_cmd)
        self.assertIn("--geno", analysis_qc_cmd)
        self.assertIn("0.02", analysis_qc_cmd)
        self.assertIn("--keep", control_hwe_cmd)
        self.assertIn("plink_controls_keep.tsv", " ".join(control_hwe_cmd))
        self.assertIn("--hwe", control_hwe_cmd)
        self.assertIn("1e-06", control_hwe_cmd)
        self.assertIn("--glm", assoc_cmd)
        self.assertIn("firth-fallback", assoc_cmd)
        self.assertIn("--covar-variance-standardize", assoc_cmd)
        self.assertIn("--threads", assoc_cmd)
        self.assertIn("8", assoc_cmd)

    def test_parse_plink_glm_results_normalizes_logistic_columns(self) -> None:
        tmp = Path(microarray_plink_output_dir(build_output_paths(self._config_and_matched()[0]), "parse_test"))
        tmp.mkdir(parents=True, exist_ok=True)
        result_path = tmp / "microarray_rhabdo.analysis_case.glm.logistic.hybrid"
        pd.DataFrame(
            [
                {
                    "#CHROM": "22",
                    "POS": "101",
                    "ID": "22:101:A:G",
                    "REF": "A",
                    "ALT": "G",
                    "A1": "G",
                    "TEST": "ADD",
                    "OBS_CT": "100",
                    "A1_FREQ": "0.12",
                    "OR": "1.5",
                    "LOG(OR)_SE": "0.2",
                    "Z_STAT": "2.0",
                    "P": "0.0455",
                    "L95": "1.01",
                    "U95": "2.22",
                },
                {
                    "#CHROM": "22",
                    "POS": "101",
                    "ID": "22:101:A:G",
                    "REF": "A",
                    "ALT": "G",
                    "A1": "G",
                    "TEST": "DOMDEV",
                    "OBS_CT": "100",
                    "A1_FREQ": "0.12",
                    "OR": "1.1",
                    "LOG(OR)_SE": "0.2",
                    "Z_STAT": "0.5",
                    "P": "0.61",
                    "L95": "0.75",
                    "U95": "1.61",
                },
            ]
        ).to_csv(result_path, sep="\t", index=False)

        parsed = parse_plink_glm_results(str(result_path))

        self.assertEqual(parsed.shape[0], 1)
        self.assertEqual(parsed.loc[0, "variant_id"], "22:101:A:G")
        self.assertEqual(parsed.loc[0, "chromosome"], "22")
        self.assertAlmostEqual(parsed.loc[0, "odds_ratio"], 1.5)
        self.assertAlmostEqual(parsed.loc[0, "beta"], 0.4054651081)
        self.assertAlmostEqual(parsed.loc[0, "regression_p"], 0.0455)

    def test_microarray_plink_cli_defaults(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["run-microarray-plink-gwas"])

        self.assertEqual(args.command, "run-microarray-plink-gwas")
        self.assertEqual(args.chromosomes, "22")
        self.assertEqual(args.min_maf, 0.05)
        self.assertEqual(args.min_mac, 20)
        self.assertEqual(args.min_call_rate, 0.98)
        self.assertEqual(args.hwe_p_control, 1e-6)
        self.assertEqual(args.analysis_split, "train")
        self.assertEqual(args.eligibility_flag, "primary_model_eligible")
        self.assertIsNone(args.copy_plink_to)
        self.assertIsNone(args.plink_prefix)
        self.assertIsNone(args.threads)
        self.assertEqual(args.plink2_bin, "plink2")
        self.assertEqual(microarray_plink_default_label(["22"]), "microarray_plink_chr22_maf05_train_qc")


if __name__ == "__main__":
    unittest.main()
