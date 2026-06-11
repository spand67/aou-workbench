from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.clinical_prs_model import (
    clinical_prs_model_coefficients_path,
    clinical_prs_model_dir,
    clinical_prs_model_metrics_path,
    clinical_prs_model_predictions_path,
    clinical_prs_model_report_path,
    run_clinical_prs_model,
)
from aou_workbench.cohort_summary import characterize_case_control_cohort, clinical_model_input_path
from aou_workbench.config import load_project_config
from aou_workbench.io_utils import read_table, write_dataframe
from aou_workbench.microarray_plink_prs import microarray_prs_scores_path, microarray_prs_weights_path
from aou_workbench.pipeline import build_cohort_artifacts, match_controls_artifacts
from tests.support import build_demo_project_tree


class ClinicalPrsModelTests(unittest.TestCase):
    def test_clinical_prs_model_uses_existing_train_scores_and_writes_outputs(self) -> None:
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

        plink_prefix = Path(paths["root"]) / "plink" / "arrays"
        plink_prefix.parent.mkdir(parents=True, exist_ok=True)
        for suffix in ("bed", "bim"):
            (plink_prefix.parent / f"arrays.{suffix}").write_text("", encoding="utf-8")
        fam = pd.DataFrame({"FID": "0", "IID": model_input["person_id"].astype(str), "father": "0", "mother": "0", "sex": "0", "phenotype": "NA"})
        fam.to_csv(plink_prefix.parent / "arrays.fam", sep=" ", index=False, header=False)

        write_dataframe(
            pd.DataFrame([{"ID": "v1", "A1": "A", "BETA": 0.2, "P": 0.001}]),
            microarray_prs_weights_path(output_paths, "gwas", "test-clumped-p001"),
        )
        test_rows = model_input[model_input["analysis_split"] == "test"].copy()
        train_rows = model_input[model_input["analysis_split"] == "train"].copy()
        outcome = effective.analysis.matched_outcome_column
        write_dataframe(
            pd.DataFrame(
                {
                    "person_id": test_rows["person_id"].astype(str),
                    "threshold_label": "p0_01",
                    "prs_score": pd.to_numeric(test_rows[outcome], errors="coerce") * 0.4 + 0.1,
                    "analysis_case": pd.to_numeric(test_rows[outcome], errors="coerce"),
                }
            ),
            microarray_prs_scores_path(output_paths, "gwas", "test-clumped-p001"),
        )
        outdir = Path(clinical_prs_model_dir(output_paths, "clinical-prs"))
        outdir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {
                "#FID": "0",
                "IID": train_rows["person_id"].astype(str),
                "ALLELE_CT": 2,
                "SCORE1_SUM": pd.to_numeric(train_rows[outcome], errors="coerce") * 0.4 + 0.1,
            }
        ).to_csv(outdir / "train_prs.sscore", sep="\t", index=False)

        outputs = run_clinical_prs_model(
            effective,
            output_paths,
            gwas_label="gwas",
            prs_label="test-clumped-p001",
            plink_prefix=str(plink_prefix),
            plink2_bin="true",
            label="clinical-prs",
        )

        metrics = outputs["metrics"]
        coefficients = outputs["coefficients"]
        self.assertEqual(set(metrics["evaluation_set"]), {"train", "test"})
        self.assertIn("prs_score_per_sd", coefficients["feature"].tolist())
        self.assertTrue(Path(clinical_prs_model_metrics_path(output_paths, "clinical-prs")).exists())
        self.assertTrue(Path(clinical_prs_model_coefficients_path(output_paths, "clinical-prs")).exists())
        self.assertTrue(Path(clinical_prs_model_predictions_path(output_paths, "clinical-prs")).exists())
        self.assertTrue(Path(clinical_prs_model_report_path(output_paths, "clinical-prs")).exists())


if __name__ == "__main__":
    unittest.main()
