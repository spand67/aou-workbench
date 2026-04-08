from __future__ import annotations

from pathlib import Path
import unittest

from aou_workbench.config import load_project_config
from aou_workbench.pipeline import run_all
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
        output_paths = run_all(config, skip_preflight=True)
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
