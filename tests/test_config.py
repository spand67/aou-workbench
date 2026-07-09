from __future__ import annotations

from pathlib import Path
import unittest

from aou_workbench.config import DEFAULT_MICROARRAY_MT_PATH, DEFAULT_MICROARRAY_PLINK_PREFIX, load_project_config
from aou_workbench.paths import build_output_paths
from tests.support import build_demo_project_tree


class ConfigLoadingTests(unittest.TestCase):
    def test_load_project_config_reads_nested_yaml(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        self.assertEqual(config.analysis.analysis_name, "unit_test_rhabdo")
        self.assertEqual(len(config.panel.a_priori_variants), 6)
        self.assertEqual(config.phenotype.definite.measurement_min, 5000)
        self.assertEqual(config.phenotype.definite.measurement_window_start_days, -7)
        self.assertEqual(config.phenotype.definite.measurement_window_end_days, 45)
        self.assertEqual(config.phenotype.broad.name, "broad")
        self.assertEqual(config.analysis.stage2.max_af, 0.001)
        self.assertEqual(config.analysis.stage3.max_af, 0.001)
        self.assertEqual(config.cohort.primary_case_tier, "broad")
        self.assertEqual(config.workbench.microarray_mt_path, DEFAULT_MICROARRAY_MT_PATH)
        self.assertEqual(config.workbench.microarray_plink_prefix, DEFAULT_MICROARRAY_PLINK_PREFIX)
        self.assertTrue(config.config_hash)

    def test_rhabdo_stage4_config_uses_primary_gwas_covariates(self) -> None:
        repo_root = Path(__file__).resolve().parents[1]
        config = load_project_config(
            workbench_path=str(repo_root / "configs" / "workbench.yaml"),
            phenotype_path=str(repo_root / "configs" / "rhabdo" / "phenotype.yaml"),
            cohort_path=str(repo_root / "configs" / "rhabdo" / "cohort.yaml"),
            panel_path=str(repo_root / "configs" / "rhabdo" / "panel.yaml"),
            analysis_path=str(repo_root / "configs" / "rhabdo" / "analysis.yaml"),
        )

        self.assertEqual(
            config.analysis.stage4.covariates,
            ("age_at_index", "is_female", "pc1", "pc2", "pc3", "pc4", "pc5"),
        )
        self.assertNotIn("preindex_crush_injury", config.analysis.stage4.covariates)
        self.assertNotIn("preindex_sepsis", config.analysis.stage4.covariates)
        self.assertNotIn("preindex_renal_injury", config.analysis.stage4.covariates)

    def test_ck_confirmed_config_uses_definite_primary_and_separate_output_root(self) -> None:
        repo_root = Path(__file__).resolve().parents[1]
        config = load_project_config(
            workbench_path=str(repo_root / "configs" / "workbench.yaml"),
            phenotype_path=str(repo_root / "configs" / "rhabdo" / "phenotype.yaml"),
            cohort_path=str(repo_root / "configs" / "rhabdo_ck_confirmed" / "cohort.yaml"),
            panel_path=str(repo_root / "configs" / "rhabdo" / "panel.yaml"),
            analysis_path=str(repo_root / "configs" / "rhabdo_ck_confirmed" / "analysis.yaml"),
        )

        self.assertEqual(config.analysis.analysis_name, "rhabdo_ck_confirmed_v1")
        self.assertEqual(config.cohort.primary_case_tier, "definite")
        self.assertEqual(config.cohort.sensitivity_case_tiers, ("broad",))
        self.assertEqual(config.phenotype.definite.measurement_min, 5000)
        self.assertEqual(config.phenotype.definite.measurement_window_start_days, -7)
        self.assertEqual(config.phenotype.definite.measurement_window_end_days, 45)
        self.assertIn("rhabdo-ck-confirmed-v1", build_output_paths(config).run_root)
        for stage in (config.analysis.stage1, config.analysis.stage2, config.analysis.stage3, config.analysis.stage4):
            self.assertEqual(stage.covariates, ("age_at_index", "is_female", "pc1", "pc2", "pc3", "pc4", "pc5"))
            self.assertNotIn("preindex_sepsis", stage.covariates)
            self.assertNotIn("preindex_renal_injury", stage.covariates)
            self.assertNotIn("preindex_crush_injury", stage.covariates)
