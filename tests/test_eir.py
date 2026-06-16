from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
import yaml

from aou_workbench.cli import _build_parser
from aou_workbench.clinical_model import (
    clinical_model_coefficients_path,
    clinical_model_metrics_path,
    clinical_model_predictions_path,
    clinical_model_report_path,
)
from aou_workbench.config import load_project_config
from aou_workbench.eir import (
    build_eir_cohort,
    build_eir_cohort_artifacts,
    characterize_eir_artifacts,
    estimate_eir_cohort_artifacts,
    eir_consort_counts_path,
    eir_missingness_path,
    eir_model_input_path,
    eir_table1_path,
    run_eir_clinical_model,
)
from aou_workbench.incident_feasibility import (
    estimate_incident_feasibility_artifacts,
    incident_case_funnel_path,
    incident_feasibility_counts_path,
    incident_feasibility_report_path,
    incident_microarray_overlap_counts_path,
    render_incident_feasibility_sql,
    run_incident_feasibility,
)


def _write_table(path: Path, rows: list[dict]) -> str:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return str(path)


def _write_yaml(path: Path, payload: dict) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return str(path)


def build_eir_demo_project_tree() -> dict[str, str]:
    root = Path(tempfile.mkdtemp(prefix="aou-eir-test-"))
    inputs = root / "inputs"
    configs = root / "configs"
    outputs = root / "outputs"
    inputs.mkdir(parents=True, exist_ok=True)
    configs.mkdir(parents=True, exist_ok=True)

    cohort_rows = []
    for person_id in range(1, 31):
        cohort_rows.append(
            {
                "person_id": person_id,
                "age": 25 + person_id,
                "sex": "female" if person_id % 2 else "male",
                "observation_start_date": "2020-01-01",
                "observation_end_date": "2024-12-31",
            }
        )

    condition_rows = []
    for person_id in range(1, 31):
        condition_rows.extend(
            [
                {
                    "person_id": person_id,
                    "condition_concept_id": 300,
                    "condition_concept_name": "Routine health maintenance",
                    "condition_start_date": "2020-06-01",
                },
                {
                    "person_id": person_id,
                    "condition_concept_id": 301,
                    "condition_concept_name": "Primary care visit",
                    "condition_start_date": "2021-01-01",
                },
                {
                    "person_id": person_id,
                    "condition_concept_id": 302,
                    "condition_concept_name": "Follow-up examination",
                    "condition_start_date": "2024-01-05",
                },
            ]
        )

    primary_cases = {1, 2, 9, 10, 11, 12, 13, 14}
    rhabdo_dates = {
        1: "2021-03-01",
        2: "2021-04-01",
        3: "2021-05-01",
        4: "2021-06-01",
        5: "2021-07-01",
        9: "2021-08-01",
        10: "2021-09-01",
        11: "2021-10-01",
        12: "2021-11-01",
        13: "2021-12-01",
        14: "2022-01-01",
    }
    for person_id, date in rhabdo_dates.items():
        condition_rows.append(
            {
                "person_id": person_id,
                "condition_concept_id": 100,
                "condition_concept_name": "Rhabdomyolysis",
                "condition_start_date": date,
            }
        )

    condition_rows.extend(
        [
            {"person_id": 2, "condition_concept_id": 410, "condition_concept_name": "Heat exhaustion", "condition_start_date": "2021-04-01"},
            {"person_id": 3, "condition_concept_id": 400, "condition_concept_name": "Crush injury", "condition_start_date": "2021-05-01"},
            {"person_id": 4, "condition_concept_id": 401, "condition_concept_name": "Sepsis", "condition_start_date": "2021-05-20"},
            {"person_id": 7, "condition_concept_id": 100, "condition_concept_name": "Rhabdomyolysis", "condition_start_date": "2020-12-01"},
            {"person_id": 15, "condition_concept_id": 420, "condition_concept_name": "Chronic kidney disease", "condition_start_date": "2020-09-01"},
            {"person_id": 16, "condition_concept_id": 421, "condition_concept_name": "Diabetes mellitus", "condition_start_date": "2020-10-01"},
            {"person_id": 17, "condition_concept_id": 422, "condition_concept_name": "Myopathy", "condition_start_date": "2020-11-01"},
            {"person_id": 18, "condition_concept_id": 423, "condition_concept_name": "Remote sepsis", "condition_start_date": "2020-08-01"},
            {"person_id": 19, "condition_concept_id": 424, "condition_concept_name": "Acute kidney injury", "condition_start_date": "2022-06-01"},
            {"person_id": 30, "condition_concept_id": 303, "condition_concept_name": "Second prior condition", "condition_start_date": "2020-08-01"},
        ]
    )

    measurement_rows = []
    for person_id in range(1, 31):
        measurement_rows.append(
            {
                "person_id": person_id,
                "measurement_concept_id": 901,
                "measurement_concept_name": "Creatinine",
                "value_as_number": 0.8 + person_id / 100.0,
                "measurement_date": "2020-12-15",
            }
        )
    for person_id in primary_cases | {3, 4}:
        rhabdo_date = pd.Timestamp(rhabdo_dates[person_id])
        measurement_rows.append(
            {
                "person_id": person_id,
                "measurement_concept_id": 900,
                "measurement_concept_name": "Creatine kinase",
                "value_as_number": 6500 + person_id,
                "measurement_date": (rhabdo_date + pd.Timedelta(days=1)).strftime("%Y-%m-%d"),
            }
        )
    measurement_rows.extend(
        [
            {"person_id": 5, "measurement_concept_id": 900, "measurement_concept_name": "Creatine kinase", "value_as_number": 1200, "measurement_date": "2021-07-02"},
            {"person_id": 6, "measurement_concept_id": 900, "measurement_concept_name": "Creatine kinase", "value_as_number": 7100, "measurement_date": "2021-08-02"},
            {"person_id": 8, "measurement_concept_id": 900, "measurement_concept_name": "Creatine kinase", "value_as_number": 9000, "measurement_date": "2020-12-15"},
            {"person_id": 20, "measurement_concept_id": 900, "measurement_concept_name": "Creatine kinase", "value_as_number": 100, "measurement_date": "2021-08-02"},
            {"person_id": 30, "measurement_concept_id": 900, "measurement_concept_name": "Creatine kinase", "value_as_number": 100, "measurement_date": "2020-12-31"},
        ]
    )

    ancestry_rows = []
    for person_id in range(1, 31):
        row = {"research_id": person_id, "ancestry_pred": "EUR" if person_id <= 15 else "AFR"}
        for pc in range(1, 11):
            row[f"PC{pc}"] = person_id * pc / 100.0
        ancestry_rows.append(row)

    paths = {
        "root": str(root),
        "output_dir": str(outputs),
        "cohort_table": _write_table(inputs / "cohort.tsv", cohort_rows),
        "condition_table": _write_table(inputs / "condition.tsv", condition_rows),
        "measurement_table": _write_table(inputs / "measurement.tsv", measurement_rows),
        "ancestry_table": _write_table(inputs / "ancestry.tsv", ancestry_rows),
    }
    workbench = {"storage_root": str(outputs), "workspace_cdr": "test.CDR"}
    phenotype = {
        "tables": {
            "cohort_table": paths["cohort_table"],
            "condition_table": paths["condition_table"],
            "measurement_table": paths["measurement_table"],
            "ancestry_table": paths["ancestry_table"],
        },
        "person_id_column": "person_id",
        "age_column": "age",
        "sex_column": "sex",
        "observation_start_column": "observation_start_date",
        "observation_end_column": "observation_end_date",
        "condition_concept_column": "condition_concept_id",
        "condition_date_column": "condition_start_date",
        "measurement_concept_column": "measurement_concept_id",
        "measurement_value_column": "value_as_number",
        "measurement_date_column": "measurement_date",
        "ancestry_person_id_column": "research_id",
        "ancestry_label_column": "ancestry_pred",
        "pc_columns": [f"PC{i}" for i in range(1, 11)],
        "definite": {
            "condition_concept_ids": [100],
            "measurement_concept_ids": [900],
            "measurement_terms": ["creatine kinase"],
            "measurement_min": 5000,
            "require_condition": True,
            "require_measurement": True,
            "measurement_window_start_days": -7,
            "measurement_window_end_days": 45,
        },
        "broad": {"condition_concept_ids": [100], "require_condition": True, "require_measurement": False},
        "clinical_cofactors": [
            {"name": "crush_injury", "condition_terms": ["crush injury"]},
            {"name": "major_trauma", "condition_terms": ["major trauma"]},
            {"name": "sepsis", "condition_terms": ["sepsis"]},
            {"name": "renal_injury", "condition_terms": ["acute kidney injury"]},
            {"name": "kidney_disease", "condition_terms": ["chronic kidney disease"]},
            {"name": "heat_illness", "condition_terms": ["heat exhaustion", "heat stroke"]},
            {"name": "dehydration", "condition_terms": ["dehydration"]},
            {"name": "myopathy_muscle_disease", "condition_terms": ["myopathy"]},
            {"name": "diabetes", "condition_terms": ["diabetes mellitus"]},
        ],
    }
    cohort = {"control_ratio": 0, "exact_match_columns": [], "primary_case_tier": "eir_primary"}
    panel = {"genes_of_interest": ["RYR1"], "burden_target_genes": ["RYR1"], "a_priori_variants": []}
    analysis = {
        "analysis_name": "eir_clinical_v1_test",
        "output_dir": str(outputs),
        "matched_outcome_column": "eir_primary_case",
        "run_stage1": False,
        "run_stage2": False,
        "run_stage3": False,
        "run_stage4": False,
    }
    paths["workbench"] = _write_yaml(configs / "workbench.yaml", workbench)
    paths["phenotype"] = _write_yaml(configs / "phenotype.yaml", phenotype)
    paths["cohort"] = _write_yaml(configs / "cohort.yaml", cohort)
    paths["panel"] = _write_yaml(configs / "panel.yaml", panel)
    paths["analysis"] = _write_yaml(configs / "analysis.yaml", analysis)
    return paths


def _config():
    paths = build_eir_demo_project_tree()
    return paths, load_project_config(
        workbench_path=paths["workbench"],
        phenotype_path=paths["phenotype"],
        cohort_path=paths["cohort"],
        panel_path=paths["panel"],
        analysis_path=paths["analysis"],
    )


class EIRPhenotypeTests(unittest.TestCase):
    def test_baseline_outcome_exclusions_and_controls(self) -> None:
        _, config = _config()
        cohort = build_eir_cohort(config).set_index("person_id")

        self.assertEqual(int(cohort.loc["1", "eir_primary_case"]), 1)
        self.assertEqual(int(cohort.loc["2", "eir_exertion_heat_dehydration_coded_case"]), 1)
        self.assertEqual(int(cohort.loc["3", "excluded_periindex_trauma"]), 1)
        self.assertEqual(int(cohort.loc["3", "eir_primary_case"]), 0)
        self.assertEqual(int(cohort.loc["4", "excluded_pre_or_same_day_sepsis"]), 1)
        self.assertEqual(int(cohort.loc["4", "eir_primary_case"]), 0)
        self.assertEqual(int(cohort.loc["5", "eir_diagnosis_only_case"]), 1)
        self.assertEqual(int(cohort.loc["5", "eir_primary_case"]), 0)
        self.assertEqual(int(cohort.loc["6", "ck_only_no_rhabdo"]), 1)
        self.assertEqual(int(cohort.loc["6", "eligible_control"]), 0)
        self.assertEqual(int(cohort.loc["7", "incident_denominator"]), 0)
        self.assertEqual(int(cohort.loc["8", "incident_denominator"]), 0)
        self.assertEqual(int(cohort.loc["15", "eligible_control"]), 1)
        self.assertEqual(int(cohort.loc["20", "eir_ck_tested_control"]), 1)
        self.assertEqual(str(pd.Timestamp(cohort.loc["30", "baseline_date"]).date()), "2020-12-31")

    def test_characterization_and_model_outputs_are_written(self) -> None:
        _, config = _config()
        effective, paths, cohort = build_eir_cohort_artifacts(config)
        _, _, characterized = characterize_eir_artifacts(effective)
        outputs = run_eir_clinical_model(effective, paths)

        self.assertFalse(cohort.empty)
        self.assertFalse(characterized["model_input"].empty)
        self.assertEqual(set(outputs["metrics"]["evaluation_set"]), {"train", "test"})
        self.assertTrue(Path(eir_consort_counts_path(paths)).exists())
        self.assertTrue(Path(eir_model_input_path(paths)).exists())
        self.assertTrue(Path(eir_table1_path(paths)).exists())
        self.assertTrue(Path(eir_missingness_path(paths)).exists())
        self.assertTrue(Path(clinical_model_metrics_path(paths)).exists())
        self.assertTrue(Path(clinical_model_coefficients_path(paths)).exists())
        self.assertTrue(Path(clinical_model_predictions_path(paths)).exists())
        self.assertTrue(Path(clinical_model_report_path(paths)).exists())

        features = "\n".join(outputs["coefficients"]["feature"].astype(str).tolist())
        self.assertNotIn("observation_days", features)
        self.assertNotIn("condition_record", features)
        self.assertNotIn("periindex", features)
        self.assertNotIn("postindex", features)

    def test_cli_help_includes_eir_commands(self) -> None:
        parser = _build_parser()
        parsed = parser.parse_args(["build-eir-cohort", "--dry-run", "--max-tib", "0.5", "--write-sql", "eir.sql"])
        self.assertEqual(parsed.command, "build-eir-cohort")
        self.assertTrue(parsed.dry_run)
        self.assertEqual(parsed.max_tib, 0.5)
        self.assertEqual(parsed.write_sql, "eir.sql")
        parsed = parser.parse_args(["characterize-eir-cohort"])
        self.assertEqual(parsed.command, "characterize-eir-cohort")
        parsed = parser.parse_args(["run-eir-clinical-model", "--run-sparse"])
        self.assertEqual(parsed.command, "run-eir-clinical-model")
        self.assertTrue(parsed.run_sparse)
        parsed = parser.parse_args(["incident-rhabdo-feasibility", "--dry-run", "--max-tib", "0.5", "--microarray-fam", "arrays.fam", "--from-cohort-tsv", "built.tsv"])
        self.assertEqual(parsed.command, "incident-rhabdo-feasibility")
        self.assertTrue(parsed.dry_run)
        self.assertEqual(parsed.max_tib, 0.5)
        self.assertEqual(parsed.microarray_fam, "arrays.fam")
        self.assertEqual(parsed.from_cohort_tsv, "built.tsv")

    def test_eir_dry_run_noops_for_local_tables(self) -> None:
        _, config = _config()
        _, _, estimate = estimate_eir_cohort_artifacts(config, max_tib=0.5)

        self.assertEqual(estimate["mode"], "local")
        self.assertEqual(estimate["total_bytes_processed"], 0)
        self.assertEqual(estimate["estimated_query_cost_usd"], 0.0)
        self.assertEqual(estimate["maximum_bytes_billed"], int(0.5 * 1024**4))

    def test_bigquery_artifact_build_streams_to_tsv(self) -> None:
        _, config = _config()
        bigquery_config = replace(
            config,
            phenotype=replace(
                config.phenotype,
                tables=replace(config.phenotype.tables, cohort_table=None),
            ),
        )

        def fake_stream(*args, **kwargs):
            output_path = Path(kwargs["stream_tsv_path"])
            output_path.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(
                [
                    {"person_id": "1", "eir_primary_case": 1, "eligible_control": 0},
                    {"person_id": "2", "eir_primary_case": 0, "eligible_control": 1},
                    {"person_id": "3", "eir_primary_case": 0, "eligible_control": 1},
                ]
            ).to_csv(output_path, sep="\t", index=False)
            return {"row_count": 3, "total_bytes_billed": 123, "path": str(output_path)}

        with patch("aou_workbench.eir._build_eir_cohort_bigquery_aggregated", side_effect=fake_stream) as mocked:
            _, paths, summary = build_eir_cohort_artifacts(bigquery_config, max_tib=1.0, write_sql_path="eir.sql")

        self.assertEqual(len(summary), 3)
        self.assertEqual(int(summary["eir_primary_case"].sum()), 1)
        self.assertEqual(int(summary["eligible_control"].sum()), 2)
        self.assertTrue(Path(paths.built_cohort_tsv).exists())
        self.assertEqual(mocked.call_args.kwargs["stream_tsv_path"], paths.built_cohort_tsv)

    def test_incident_feasibility_local_outputs_and_keeps_sepsis_as_flag(self) -> None:
        paths_dict, config = _config()
        fam_path = Path(paths_dict["root"]) / "arrays.fam"
        fam_path.write_text(
            "\n".join(
                [
                    "0 1 0 0 0 NA",
                    "0 2 0 0 0 NA",
                    "0 4 0 0 0 NA",
                    "0 15 0 0 0 NA",
                    "0 20 0 0 0 NA",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

        _, paths, outputs = run_incident_feasibility(config, microarray_fam=str(fam_path))
        counts = outputs["feasibility_counts"].set_index("metric")["n"].astype(int)
        case_funnel = outputs["case_funnel"].set_index("metric")["n"].astype(int)
        microarray = outputs["microarray_overlap_counts"].set_index("metric")["n"].astype(int)

        self.assertEqual(int(counts["Primary non-trauma incident cases"]), 10)
        self.assertEqual(int(counts["CK-confirmed non-trauma cases"]), 9)
        self.assertEqual(int(case_funnel["Peri-index sepsis-flagged non-trauma cases"]), 1)
        self.assertEqual(int(microarray["Primary non-trauma incident cases with microarray"]), 3)
        self.assertEqual(int(microarray["Peri-index sepsis-flagged cases with microarray"]), 1)
        self.assertTrue(Path(incident_feasibility_counts_path(paths)).exists())
        self.assertTrue(Path(incident_case_funnel_path(paths)).exists())
        self.assertTrue(Path(incident_microarray_overlap_counts_path(paths)).exists())
        self.assertTrue(Path(incident_feasibility_report_path(paths)).exists())

    def test_incident_feasibility_sql_contains_primary_rules(self) -> None:
        _, config = _config()
        sql = render_incident_feasibility_sql(config)

        self.assertIn("INTERVAL 365 DAY", sql)
        self.assertIn("omop_condition_record_dates >= 2", sql)
        self.assertIn("prior_rhabdo_date IS NULL", sql)
        self.assertIn("prior_high_ck = 0", sql)
        self.assertIn("excluded_periindex_trauma = 0", sql)
        self.assertIn("periindex_sepsis_flag", sql)
        self.assertIn("ck_confirmed_nontrauma_case", sql)
        self.assertIn("eligible_control", sql)
        self.assertIn("summary AS", sql)
        self.assertIn("UNNEST", sql)

    def test_incident_feasibility_dry_run_uses_aggregate_sql(self) -> None:
        paths_dict, config = _config()
        bigquery_config = replace(
            config,
            phenotype=replace(
                config.phenotype,
                tables=replace(config.phenotype.tables, cohort_table=None),
            ),
        )
        sql_path = str(Path(paths_dict["root"]) / "incident.sql")

        with patch("aou_workbench.incident_feasibility.dry_run_bigquery_query", return_value={"total_bytes_processed": 1024, "total_tib_processed": 1024 / float(1024**4), "maximum_bytes_billed": None, "would_exceed_maximum_bytes_billed": False}) as mocked:
            _, _, estimate = estimate_incident_feasibility_artifacts(bigquery_config, max_tib=1.0, write_sql_path=sql_path)

        self.assertEqual(estimate["mode"], "bigquery")
        self.assertEqual(estimate["estimated_query_cost_usd"], 1024 / float(1024**4) * 6.25)
        self.assertIn("case_funnel", mocked.call_args.args[0])

    def test_incident_feasibility_can_reuse_existing_cohort_tsv(self) -> None:
        paths_dict, config = _config()
        cohort = build_eir_cohort(config)
        cohort_path = Path(paths_dict["root"]) / "built_cohort.tsv"
        cohort.to_csv(cohort_path, sep="\t", index=False)

        _, paths, outputs = run_incident_feasibility(config, from_cohort_tsv=str(cohort_path))

        counts = outputs["feasibility_counts"].set_index("metric")["n"].astype(int)
        self.assertEqual(int(counts["Primary non-trauma incident cases"]), 10)
        self.assertTrue(Path(incident_feasibility_counts_path(paths)).exists())


if __name__ == "__main__":
    unittest.main()
