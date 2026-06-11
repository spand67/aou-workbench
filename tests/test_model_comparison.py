from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.clinical_model import run_clinical_model
from aou_workbench.cohort_summary import characterize_case_control_cohort
from aou_workbench.config import load_project_config
from aou_workbench.io_utils import write_dataframe
from aou_workbench.microarray_plink_prs import microarray_prs_scores_path
from aou_workbench.model_comparison import (
    model_comparison_metrics_path,
    model_comparison_predictions_path,
    model_comparison_report_path,
    run_heldout_model_comparison,
)
from aou_workbench.pipeline import build_cohort_artifacts, match_controls_artifacts
from tests.support import build_demo_project_tree


class ModelComparisonTests(unittest.TestCase):
    def test_model_comparison_writes_clinical_and_prs_metrics(self) -> None:
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
        clinical_outputs = run_clinical_model(effective, output_paths)
        test_predictions = clinical_outputs["predictions"]
        test_predictions = test_predictions[test_predictions["analysis_split"] == "test"].copy()
        outcome = effective.analysis.matched_outcome_column
        prs_scores = pd.DataFrame(
            {
                "person_id": test_predictions["person_id"].astype(str),
                "threshold_label": "p0_01",
                "prs_score": pd.to_numeric(test_predictions[outcome], errors="coerce") * 0.5,
                "analysis_case": pd.to_numeric(test_predictions[outcome], errors="coerce"),
            }
        )
        write_dataframe(prs_scores, microarray_prs_scores_path(output_paths, "gwas", "prs"))

        outputs = run_heldout_model_comparison(
            effective,
            output_paths,
            gwas_label="gwas",
            prs_label="prs",
            label="comparison",
        )

        metrics = outputs["metrics"]
        self.assertEqual(set(metrics["model"]), {"clinical_only", "prs_only"})
        self.assertEqual(set(metrics["threshold_label"]), {"p0_01"})
        self.assertTrue(Path(model_comparison_metrics_path(output_paths, "comparison")).exists())
        self.assertTrue(Path(model_comparison_predictions_path(output_paths, "comparison")).exists())
        self.assertTrue(Path(model_comparison_report_path(output_paths, "comparison")).exists())
        report = Path(model_comparison_report_path(output_paths, "comparison")).read_text()
        self.assertIn("No clinical+PRS combiner is trained here", report)


if __name__ == "__main__":
    unittest.main()
