from __future__ import annotations

from dataclasses import replace
import unittest

import pandas as pd

from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.config import load_project_config
from aou_workbench.matching import match_case_controls
from aou_workbench.stage1_prepare import stage1_sample_manifest_path
from aou_workbench.stage4_hail_gwas import _hail_sample_frame, _normalize_chromosomes
from tests.support import build_demo_project_tree


class Stage4HailGwasTests(unittest.TestCase):
    def test_normalize_chromosomes_defaults_and_deduplicates(self) -> None:
        self.assertEqual(_normalize_chromosomes(None), [str(chrom) for chrom in range(1, 23)])
        self.assertEqual(_normalize_chromosomes(["chr1", "1", "19", "chr19", "22"]), ["1", "19", "22"])

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

        sample_df, hail_covariates, raw_covariates = _hail_sample_frame(matched_df, config)

        self.assertLess(len(sample_df), matched_df["person_id"].astype(str).nunique())
        self.assertIn("age_at_index", hail_covariates)
        self.assertIn("pc1", hail_covariates)
        self.assertIn("ancestry_pred", raw_covariates)
        self.assertNotIn("ancestry_pred", hail_covariates)
        self.assertTrue(any(column.startswith("ancestry_") for column in hail_covariates))
        self.assertFalse(sample_df[hail_covariates + [config.analysis.matched_outcome_column]].isna().any().any())
        self.assertTrue(set(sample_df["person_id"].astype(str)).issubset({str(person_id) for person_id in range(1, 11)}))


if __name__ == "__main__":
    unittest.main()
