from __future__ import annotations

from pathlib import Path
import unittest
from unittest import mock

import pandas as pd

from aou_workbench.config import load_project_config
from aou_workbench.pipeline import run_all
from aou_workbench.stage1_prepare import stage1_sample_manifest_path
from aou_workbench.stage2_prepare import stage2_sample_manifest_path
from tests.support import build_demo_project_tree


class PipelineIntegrationTests(unittest.TestCase):
    def test_run_all_writes_expected_outputs(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        cohort_df = pd.read_csv(paths["cohort_table"], sep="\t")
        pd.DataFrame({"person_id": cohort_df["person_id"].astype(str)}).to_csv(
            stage1_sample_manifest_path(paths["stage1_table"]),
            sep="\t",
            index=False,
        )
        pd.DataFrame({"person_id": cohort_df["person_id"].astype(str)}).to_csv(
            stage2_sample_manifest_path(paths["stage2_table"]),
            sep="\t",
            index=False,
        )
        with (
            mock.patch("aou_workbench.pipeline.prepare_stage1_variant_table") as mock_prepare_stage1,
            mock.patch("aou_workbench.pipeline.prepare_stage2_variant_table") as mock_prepare_stage2,
            mock.patch("aou_workbench.pipeline.prepare_stage4_acaf_subset") as mock_prepare_stage4,
        ):
            mock_prepare_stage1.return_value = pd.read_csv(paths["stage1_table"], sep="\t")
            mock_prepare_stage2.return_value = pd.read_csv(paths["stage2_table"], sep="\t")
            mock_prepare_stage4.return_value = {}
            output_paths = run_all(config, skip_preflight=True)
        mock_prepare_stage1.assert_called_once()
        mock_prepare_stage2.assert_called_once()
        mock_prepare_stage4.assert_called_once()
        expected = [
            output_paths.built_cohort_tsv,
            output_paths.matched_cohort_tsv,
            output_paths.stage1_results_tsv,
            output_paths.stage2_gene_tsv,
            output_paths.stage3_results_tsv,
            output_paths.stage4_full_results_tsv,
            output_paths.stage4_manhattan_svg,
            output_paths.final_report_md,
        ]
        for path in expected:
            self.assertTrue(Path(path).exists(), path)
