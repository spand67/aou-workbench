from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.cli import _build_parser
from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.config import load_project_config
from aou_workbench.matching import match_case_controls
from aou_workbench.microarray_bigsnpr import (
    BIGSNPR_DEFAULT_ALPHAS,
    build_microarray_bigsnpr_plink_commands,
    microarray_bigsnpr_default_label,
    microarray_bigsnpr_metadata_path,
    microarray_bigsnpr_script_path,
    write_bigsnpr_r_script,
    write_microarray_bigsnpr_inputs,
)
from aou_workbench.paths import build_output_paths
from tests.support import build_demo_project_tree


class MicroarrayBigsnprTests(unittest.TestCase):
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
        matched_df.loc[matched_df.index[::3], "analysis_split"] = "test"
        matched_df["primary_model_eligible"] = 1
        matched_df["train_cv_fold"] = "train_cv_fold_1"
        matched_df.loc[matched_df["analysis_split"] == "test", "train_cv_fold"] = ""
        return config, matched_df

    def test_bigsnpr_inputs_write_train_test_and_control_keep_files(self) -> None:
        config, matched_df = self._config_and_matched()
        output_paths = build_output_paths(config)
        label = "bigsnpr_test"
        fam_df = pd.DataFrame(
            {
                "FID": ["0"] * matched_df.shape[0],
                "IID": matched_df["person_id"].astype(str).tolist(),
                "father": ["0"] * matched_df.shape[0],
                "mother": ["0"] * matched_df.shape[0],
                "sex": ["0"] * matched_df.shape[0],
                "phenotype": ["NA"] * matched_df.shape[0],
            }
        )

        metadata, covariates, dropped_covariates, counts = write_microarray_bigsnpr_inputs(
            config,
            matched_df,
            fam_df,
            output_paths,
            label=label,
        )

        train_ids = set(metadata.loc[metadata["analysis_split"] == "train", "IID"])
        test_ids = set(metadata.loc[metadata["analysis_split"] == "test", "IID"])
        self.assertTrue(train_ids)
        self.assertTrue(test_ids)
        self.assertFalse(train_ids.intersection(test_ids))
        self.assertTrue(set(covariates).issubset({"age_at_index", "is_female", "pc1", "pc2", "pc3", "pc4", "pc5"}))
        self.assertIn("age_at_index", covariates)
        self.assertIn("is_female", dropped_covariates)
        self.assertGreater(counts["after_microarray_fam_overlap_train_participants"], 0)
        self.assertGreater(counts["after_microarray_fam_overlap_test_participants"], 0)
        self.assertTrue(Path(microarray_bigsnpr_metadata_path(output_paths, label)).exists())

    def test_bigsnpr_plink_commands_build_qc_and_subset(self) -> None:
        config, _ = self._config_and_matched()
        output_paths = build_output_paths(config)
        label = "bigsnpr_test"
        commands = build_microarray_bigsnpr_plink_commands(
            plink2_bin="plink2",
            plink_prefix="/tmp/arrays",
            paths=output_paths,
            label=label,
            chromosomes=["22"],
            min_maf=0.05,
            min_mac=20,
            min_call_rate=0.98,
            hwe_p_control=1e-6,
            threads=8,
            memory_mb=32000,
        )
        analysis_qc_cmd, control_hwe_cmd, make_bed_cmd = commands

        self.assertIn("--snps-only", analysis_qc_cmd)
        self.assertIn("--maf", analysis_qc_cmd)
        self.assertIn("0.05", analysis_qc_cmd)
        self.assertIn("--hwe", control_hwe_cmd)
        self.assertIn("train_controls_keep.tsv", " ".join(control_hwe_cmd))
        self.assertIn("--make-bed", make_bed_cmd)
        self.assertIn("train_test_keep.tsv", " ".join(make_bed_cmd))
        self.assertIn("bigsnpr_array_subset", " ".join(make_bed_cmd))

    def test_bigsnpr_r_script_uses_sparse_logistic_model(self) -> None:
        config, _ = self._config_and_matched()
        output_paths = build_output_paths(config)
        label = "bigsnpr_script_test"

        script_path = write_bigsnpr_r_script(
            output_paths,
            label=label,
            outcome_column="analysis_case",
            covariates=["age_at_index", "is_female", "pc1"],
            alphas=BIGSNPR_DEFAULT_ALPHAS,
            folds=5,
            nlambda=100,
            dfmax=50000,
            ncores=4,
        )
        script = Path(script_path).read_text(encoding="utf-8")

        self.assertIn("bigsnpr::snp_readBed", script)
        self.assertIn("bigsnpr::snp_fastImputeSimple", script)
        self.assertIn('method = "mean2"', script)
        self.assertIn("bigstatsr::big_spLogReg", script)
        self.assertIn("pf.covar = rep(0, ncol(covar_train))", script)
        self.assertIn("predict(fit, G, ind.row = test_ind", script)
        self.assertIn("selected_variants.tsv", script)
        self.assertTrue(Path(microarray_bigsnpr_script_path(output_paths, label)).exists())

    def test_bigsnpr_cli_defaults(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["run-microarray-bigsnpr-model"])

        self.assertEqual(args.command, "run-microarray-bigsnpr-model")
        self.assertEqual(args.chromosomes, "22")
        self.assertEqual(args.min_maf, 0.05)
        self.assertEqual(args.min_mac, 20)
        self.assertEqual(args.min_call_rate, 0.98)
        self.assertEqual(args.hwe_p_control, 1e-6)
        self.assertEqual(args.eligibility_flag, "primary_model_eligible")
        self.assertFalse(args.prepare_only)
        self.assertEqual(microarray_bigsnpr_default_label(["22"]), "microarray_bigsnpr_chr22_maf05_train_qc")


if __name__ == "__main__":
    unittest.main()
