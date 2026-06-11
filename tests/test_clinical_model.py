from __future__ import annotations

from pathlib import Path
import unittest

from aou_workbench.clinical_model import (
    clinical_model_coefficients_path,
    clinical_model_metrics_path,
    clinical_model_predictions_path,
    clinical_model_report_path,
    run_clinical_model,
)
from aou_workbench.cohort_summary import characterize_case_control_cohort
from aou_workbench.config import load_project_config
from aou_workbench.pipeline import build_cohort_artifacts, match_controls_artifacts
from tests.support import build_demo_project_tree


class ClinicalModelTests(unittest.TestCase):
    def test_clinical_model_writes_metrics_predictions_and_coefficients(self) -> None:
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

        outputs = run_clinical_model(effective, output_paths)

        metrics = outputs["metrics"]
        coefficients = outputs["coefficients"]
        predictions = outputs["predictions"]
        self.assertEqual(set(metrics["evaluation_set"]), {"train", "test"})
        self.assertIn("roc_auc", metrics.columns)
        self.assertFalse(coefficients.empty)
        self.assertFalse(predictions.empty)
        self.assertTrue(Path(clinical_model_metrics_path(output_paths)).exists())
        self.assertTrue(Path(clinical_model_coefficients_path(output_paths)).exists())
        self.assertTrue(Path(clinical_model_predictions_path(output_paths)).exists())
        self.assertTrue(Path(clinical_model_report_path(output_paths)).exists())

        feature_text = "\n".join(coefficients["feature"].astype(str).tolist())
        report_text = Path(clinical_model_report_path(output_paths)).read_text()
        self.assertNotIn("sepsis", feature_text)
        self.assertNotIn("renal", feature_text)
        self.assertNotIn("remote_preindex", feature_text)
        self.assertNotIn("near_preindex", feature_text)
        self.assertNotIn("periindex", feature_text)
        self.assertNotIn("postindex", feature_text)
        self.assertNotIn("observation_days", feature_text)
        self.assertNotIn("condition_record", feature_text)
        self.assertIn("Primary predictors: age, sex category, ancestry category", report_text)


if __name__ == "__main__":
    unittest.main()
