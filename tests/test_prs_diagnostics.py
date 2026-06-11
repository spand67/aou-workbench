from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.clinical_prs_model import clinical_prs_model_predictions_path
from aou_workbench.cohort_summary import characterize_case_control_cohort, clinical_model_input_path
from aou_workbench.config import load_project_config
from aou_workbench.io_utils import read_table, write_dataframe
from aou_workbench.microarray_plink_prs import microarray_prs_scores_path
from aou_workbench.pipeline import build_cohort_artifacts, match_controls_artifacts
from aou_workbench.prs_diagnostics import (
    prs_diagnostics_ancestry_path,
    prs_diagnostics_bootstrap_ci_path,
    prs_diagnostics_calibration_path,
    prs_diagnostics_cofactor_path,
    prs_diagnostics_deciles_path,
    prs_diagnostics_definite_path,
    prs_diagnostics_overall_metrics_path,
    prs_diagnostics_report_path,
    run_prs_diagnostics,
)
from tests.support import build_demo_project_tree


class PrsDiagnosticsTests(unittest.TestCase):
    def test_prs_diagnostics_writes_existing_output_summaries(self) -> None:
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
        model_input = read_table(clinical_model_input_path(output_paths))
        outcome = effective.analysis.matched_outcome_column
        test = model_input[model_input["analysis_split"] == "test"].copy()
        train = model_input[model_input["analysis_split"] == "train"].copy()

        write_dataframe(
            pd.DataFrame(
                {
                    "person_id": test["person_id"].astype(str),
                    "threshold_label": "p0_01",
                    "prs_score": pd.to_numeric(test[outcome], errors="coerce") * 0.5 + 0.1,
                    "analysis_case": pd.to_numeric(test[outcome], errors="coerce"),
                }
            ),
            microarray_prs_scores_path(output_paths, "gwas", "prs"),
        )
        predictions = pd.concat([train, test], ignore_index=True)
        predictions = pd.DataFrame(
            {
                "person_id": predictions["person_id"].astype(str),
                outcome: pd.to_numeric(predictions[outcome], errors="coerce"),
                "analysis_split": predictions["analysis_split"].astype(str),
                "threshold_label": "p0_01",
                "predicted_probability": pd.to_numeric(predictions[outcome], errors="coerce") * 0.6 + 0.2,
            }
        )
        write_dataframe(predictions, clinical_prs_model_predictions_path(output_paths, "clinical-prs"))

        outputs = run_prs_diagnostics(
            effective,
            output_paths,
            gwas_label="gwas",
            prs_label="prs",
            clinical_prs_label="clinical-prs",
            label="diag",
            bootstrap_iterations=5,
            seed=1,
        )

        self.assertFalse(outputs["overall"].empty)
        self.assertFalse(outputs["deciles"].empty)
        self.assertTrue(Path(prs_diagnostics_overall_metrics_path(output_paths, "diag")).exists())
        self.assertTrue(Path(prs_diagnostics_bootstrap_ci_path(output_paths, "diag")).exists())
        self.assertTrue(Path(prs_diagnostics_deciles_path(output_paths, "diag")).exists())
        self.assertTrue(Path(prs_diagnostics_ancestry_path(output_paths, "diag")).exists())
        self.assertTrue(Path(prs_diagnostics_definite_path(output_paths, "diag")).exists())
        self.assertTrue(Path(prs_diagnostics_cofactor_path(output_paths, "diag")).exists())
        self.assertTrue(Path(prs_diagnostics_calibration_path(output_paths, "diag")).exists())
        self.assertTrue(Path(prs_diagnostics_report_path(output_paths, "diag")).exists())


if __name__ == "__main__":
    unittest.main()
