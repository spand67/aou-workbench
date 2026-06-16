from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.clinical_model import clinical_model_metrics_path
from aou_workbench.clinical_prs_model import (
    clinical_prs_model_coefficients_path,
    clinical_prs_model_dir,
    clinical_prs_model_metrics_path,
)
from aou_workbench.cli import _build_parser
from aou_workbench.cohort_summary import characterize_case_control_cohort, clinical_model_input_path
from aou_workbench.config import load_project_config
from aou_workbench.io_utils import read_table, write_dataframe
from aou_workbench.microarray_plink_gwas import (
    microarray_plink_manhattan_path,
    microarray_plink_qc_path,
    microarray_plink_results_path,
    microarray_plink_variant_qc_summary_path,
)
from aou_workbench.microarray_plink_prs import microarray_prs_metrics_path, microarray_prs_weights_path
from aou_workbench.paths import build_output_paths
from aou_workbench.pipeline import build_cohort_artifacts, match_controls_artifacts
from aou_workbench.presentation_dashboard import (
    _cofactor_timing_summary,
    presentation_consort_svg_path,
    presentation_timing_svg_path,
    presentation_timing_table_path,
    presentation_train_prs_summary_path,
    presentation_train_prs_svg_path,
    render_presentation_dashboard,
)
from tests.support import build_demo_project_tree


def _load_demo_config(paths: dict[str, str]):
    return load_project_config(
        workbench_path=paths["workbench"],
        phenotype_path=paths["phenotype"],
        cohort_path=paths["cohort"],
        panel_path=paths["panel"],
        analysis_path=paths["analysis"],
    )


def _build_characterized_demo():
    paths = build_demo_project_tree()
    config = _load_demo_config(paths)
    effective, output_paths, cohort_df = build_cohort_artifacts(config)
    _, _, matched_df = match_controls_artifacts(effective, cohort_df)
    characterize_case_control_cohort(effective, cohort_df, matched_df, output_paths)
    return paths, effective, output_paths


def _write_representative_dashboard_inputs(output_paths, *, gwas_label: str, prs_label: str, clinical_prs_label: str) -> None:
    model_input = read_table(clinical_model_input_path(output_paths))
    train = model_input[model_input["analysis_split"].astype(str) == "train"].copy()
    outcome = "analysis_case"

    write_dataframe(
        pd.DataFrame(
            [
                {"evaluation_set": "train", "n": 8, "cases": 2, "controls": 6, "roc_auc": 0.61, "average_precision": 0.30, "brier_score": 0.15, "log_loss": 0.48},
                {"evaluation_set": "test", "n": 4, "cases": 1, "controls": 3, "roc_auc": 0.55, "average_precision": 0.28, "brier_score": 0.17, "log_loss": 0.51},
            ]
        ),
        clinical_model_metrics_path(output_paths),
    )
    write_dataframe(
        pd.DataFrame(
            [
                {
                    "threshold_label": "p0_01",
                    "p_threshold": 0.01,
                    "n_participants": 4,
                    "n_cases": 1,
                    "n_controls": 3,
                    "n_variants_scored": 2,
                    "roc_auc": 0.57,
                    "average_precision": 0.29,
                }
            ]
        ),
        microarray_prs_metrics_path(output_paths, gwas_label, prs_label),
    )
    write_dataframe(
        pd.DataFrame(
            [
                {"evaluation_set": "train", "threshold_label": "p0_01", "n": 8, "cases": 2, "controls": 6, "roc_auc": 0.74, "average_precision": 0.42, "brier_score": 0.12, "log_loss": 0.41},
                {"evaluation_set": "test", "threshold_label": "p0_01", "n": 4, "cases": 1, "controls": 3, "roc_auc": 0.60, "average_precision": 0.33, "brier_score": 0.16, "log_loss": 0.49},
            ]
        ),
        clinical_prs_model_metrics_path(output_paths, clinical_prs_label),
    )
    write_dataframe(
        pd.DataFrame(
            [
                {"feature": "intercept", "source_column": "", "kind": "intercept", "reference": "", "beta": -1.0, "odds_ratio": 0.37},
                {"feature": "prs_score_per_sd", "source_column": "prs_score", "kind": "continuous", "reference": "", "beta": 0.8, "odds_ratio": 2.23},
                {"feature": "age_at_index_per_sd", "source_column": "age_at_index", "kind": "continuous", "reference": "", "beta": 0.2, "odds_ratio": 1.22},
            ]
        ),
        clinical_prs_model_coefficients_path(output_paths, clinical_prs_label),
    )

    Path(clinical_prs_model_dir(output_paths, clinical_prs_label)).mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "#FID": "0",
            "IID": train["person_id"].astype(str),
            "ALLELE_CT": 2,
            "SCORE1_SUM": pd.to_numeric(train[outcome], errors="coerce") * 0.4 + 0.1,
        }
    ).to_csv(Path(clinical_prs_model_dir(output_paths, clinical_prs_label)) / "train_prs.sscore", sep="\t", index=False)

    write_dataframe(
        pd.DataFrame(
            [
                {"variant_id": "rs1", "chromosome": "1", "position": 100, "effect_allele": "A", "beta": 0.2, "regression_p": 0.001},
                {"variant_id": "rs2", "chromosome": "22", "position": 200, "effect_allele": "G", "beta": -0.1, "regression_p": 0.02},
            ]
        ),
        microarray_plink_results_path(output_paths, gwas_label),
    )
    write_dataframe(
        pd.DataFrame(
            [
                {"filter": "loaded", "rows_before": 100, "rows_after": 100, "rows_removed": 0},
                {"filter": "final", "rows_before": 100, "rows_after": 2, "rows_removed": 98},
            ]
        ),
        microarray_plink_variant_qc_summary_path(output_paths, gwas_label),
    )
    write_dataframe(
        pd.DataFrame([{"ID": "rs1", "A1": "A", "BETA": 0.2, "P": 0.001}, {"ID": "rs2", "A1": "G", "BETA": -0.1, "P": 0.02}]),
        microarray_prs_weights_path(output_paths, gwas_label, prs_label),
    )
    Path(microarray_plink_manhattan_path(output_paths, gwas_label)).parent.mkdir(parents=True, exist_ok=True)
    Path(microarray_plink_manhattan_path(output_paths, gwas_label)).write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="300" height="120"><text x="10" y="30">Manhattan</text></svg>',
        encoding="utf-8",
    )
    Path(microarray_plink_qc_path(output_paths, gwas_label)).write_text(
        json.dumps(
            {
                "chromosomes_tested": ["1", "22"],
                "sample_counts": {
                    "after_microarray_fam_overlap_cases": 2,
                    "after_microarray_fam_overlap_controls": 6,
                },
                "n_variants_tested": 2,
                "min_maf": 0.05,
                "min_call_rate": 0.98,
            }
        ),
        encoding="utf-8",
    )


class PresentationDashboardTests(unittest.TestCase):
    def test_cli_parser_has_presentation_dashboard_defaults(self) -> None:
        args = _build_parser().parse_args(["presentation-dashboard"])
        self.assertEqual(args.command, "presentation-dashboard")
        self.assertEqual(args.gwas_label, "microarray_plink_autosomes_maf05_train_qc")
        self.assertEqual(args.prs_label, "test-clumped-p001")
        self.assertEqual(args.clinical_prs_label, "clinical_prs_p001")

    def test_presentation_dashboard_renders_demo_outputs(self) -> None:
        _, config, output_paths = _build_characterized_demo()
        _write_representative_dashboard_inputs(output_paths, gwas_label="gwas", prs_label="prs", clinical_prs_label="clinical-prs")

        dashboard = render_presentation_dashboard(
            config,
            output_paths,
            gwas_label="gwas",
            prs_label="prs",
            clinical_prs_label="clinical-prs",
            diagnostics_label="diag",
        )

        text = Path(dashboard).read_text(encoding="utf-8")
        self.assertIn("Rhabdomyolysis Analysis Dashboard", text)
        self.assertIn("CONSORT Flow", text)
        self.assertIn("Case And Control Definitions", text)
        self.assertIn("Training Set Table 1", text)
        self.assertIn("Testing Set Table 1", text)
        self.assertIn("Sepsis And Renal Injury Timing", text)
        self.assertIn("GWAS Parameters", text)
        self.assertIn("Training PRS Distribution", text)
        self.assertIn("Clinical + PRS Model Drivers", text)
        self.assertTrue(Path(presentation_consort_svg_path(output_paths)).exists())
        self.assertTrue(Path(presentation_timing_svg_path(output_paths)).exists())
        self.assertTrue(Path(presentation_timing_table_path(output_paths)).exists())
        self.assertTrue(Path(presentation_train_prs_svg_path(output_paths)).exists())
        self.assertTrue(Path(presentation_train_prs_summary_path(output_paths)).exists())

    def test_presentation_dashboard_labels_definite_primary_as_ck_confirmed(self) -> None:
        _, config, output_paths = _build_characterized_demo()
        config = replace(
            config,
            cohort=replace(config.cohort, primary_case_tier="definite", sensitivity_case_tiers=("broad",)),
            analysis=replace(config.analysis, analysis_name="rhabdo_ck_confirmed_v1"),
        )
        _write_representative_dashboard_inputs(output_paths, gwas_label="gwas", prs_label="prs", clinical_prs_label="clinical-prs")

        dashboard = render_presentation_dashboard(
            config,
            output_paths,
            gwas_label="gwas",
            prs_label="prs",
            clinical_prs_label="clinical-prs",
            diagnostics_label="diag",
        )

        text = Path(dashboard).read_text(encoding="utf-8")
        consort_svg = Path(presentation_consort_svg_path(output_paths)).read_text(encoding="utf-8")
        self.assertIn("CK-Confirmed Rhabdomyolysis Susceptibility Dashboard", text)
        self.assertIn("CK-confirmed cases", text)
        self.assertIn("Primary case tier: <code>definite</code>", text)
        self.assertIn("sepsis and renal injury are characterized but not excluded", text)
        self.assertIn("CK-confirmed rhabdomyolysis cases", consort_svg)

    def test_timing_summary_uses_participants_and_control_index_dates(self) -> None:
        paths = build_demo_project_tree()
        config = _load_demo_config(paths)
        output_paths = build_output_paths(config)
        condition_path = Path(paths["condition_table"])
        condition = pd.read_csv(condition_path, sep="\t")
        condition = pd.concat(
            [
                condition,
                pd.DataFrame(
                    [
                        {"person_id": 100, "condition_concept_id": 204, "condition_concept_name": "Sepsis", "condition_start_date": "2021-12-15"},
                        {"person_id": 100, "condition_concept_id": 205, "condition_concept_name": "Acute kidney injury", "condition_start_date": "2022-01-10"},
                        {"person_id": 200, "condition_concept_id": 204, "condition_concept_name": "Sepsis", "condition_start_date": "2022-01-10"},
                        {"person_id": 200, "condition_concept_id": 205, "condition_concept_name": "Acute kidney injury", "condition_start_date": "2022-04-20"},
                    ]
                ),
            ],
            ignore_index=True,
        )
        condition.to_csv(condition_path, sep="\t", index=False)
        model_input = pd.DataFrame(
            [
                {"person_id": "100", "analysis_case": 1, "index_date": "2022-01-10"},
                {"person_id": "200", "analysis_case": 0, "index_date": "2022-01-10"},
            ]
        )

        summary = _cofactor_timing_summary(config, output_paths, model_input)

        case_sepsis = summary[
            (summary["case_control"] == "cases")
            & (summary["cofactor"] == "sepsis")
            & (summary["timing_bin"] == "-30 to -8")
        ]
        control_sepsis = summary[
            (summary["case_control"] == "controls")
            & (summary["cofactor"] == "sepsis")
            & (summary["timing_bin"] == "index day")
        ]
        control_renal = summary[
            (summary["case_control"] == "controls")
            & (summary["cofactor"] == "renal_injury")
            & (summary["timing_bin"] == "+91 to +365")
        ]
        self.assertEqual(int(case_sepsis["n_participants"].iloc[0]), 1)
        self.assertEqual(int(control_sepsis["n_participants"].iloc[0]), 1)
        self.assertEqual(int(control_renal["n_participants"].iloc[0]), 1)

    def test_missing_optional_outputs_render_placeholders(self) -> None:
        _, config, output_paths = _build_characterized_demo()
        dashboard = render_presentation_dashboard(
            config,
            output_paths,
            gwas_label="missing-gwas",
            prs_label="missing-prs",
            clinical_prs_label="missing-clinical-prs",
            diagnostics_label="missing-diag",
        )
        text = Path(dashboard).read_text(encoding="utf-8")
        self.assertIn("Expected file not found", text)
        self.assertIn("No rows available", text)


if __name__ == "__main__":
    unittest.main()
