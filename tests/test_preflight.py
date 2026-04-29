from __future__ import annotations

from dataclasses import replace
import unittest
from unittest import mock

from aou_workbench.config import load_project_config
from aou_workbench.preflight import PreflightCheck, RuntimeDefaults, run_preflight_checks
from tests.support import build_demo_project_tree


class PreflightTests(unittest.TestCase):
    def test_stage1_preflight_checks_vds_and_planned_output_without_requiring_existing_table(self) -> None:
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
            workbench=replace(config.workbench, max_unrelated_path=paths["max_unrelated"]),
            analysis=replace(
                config.analysis,
                run_stage1=False,
                stage1=replace(config.analysis.stage1, variant_table="/tmp/new-stage1-output.tsv"),
            ),
        )

        runtime = RuntimeDefaults(
            workspace_bucket="gs://test-workspace-bucket",
            workspace_cdr="test-project.C2024Q3R9",
            requester_pays_project="test-project",
        )

        def named_pass(name: str) -> PreflightCheck:
            return PreflightCheck(name=name, status="PASS", message="ok")

        with (
            mock.patch("aou_workbench.preflight.discover_runtime_defaults", return_value=runtime),
            mock.patch(
                "aou_workbench.preflight._check_local_or_gcs_path",
                side_effect=lambda _path, name, requester_pays_project=None: named_pass(name),
            ),
            mock.patch(
                "aou_workbench.preflight._bigquery_table_check",
                side_effect=lambda cdr, table_reference, name: named_pass(name),
            ),
            mock.patch("aou_workbench.preflight._bigquery_check", return_value=named_pass("bigquery:person")),
            mock.patch(
                "aou_workbench.preflight._tool_check",
                side_effect=lambda command, name, required_for, override_env=None: named_pass(name),
            ),
            mock.patch("aou_workbench.preflight._hail_check", return_value=named_pass("hail:available")),
        ):
            checks = run_preflight_checks(config)

        names = [check.name for check in checks]
        self.assertIn("input:wgs_vds", names)
        self.assertIn("output:stage1", names)
        self.assertIn("input:clinvar_mt", names)
        self.assertIn("output:stage2", names)
        self.assertIn("input:acaf_mt", names)
        self.assertIn("input:max_unrelated", names)
        self.assertIn("output:stage4", names)
        self.assertNotIn("input:stage1", names)
        self.assertNotIn("input:stage2", names)


if __name__ == "__main__":
    unittest.main()
