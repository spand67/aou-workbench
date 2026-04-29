from __future__ import annotations

from pathlib import Path
import tempfile

import pandas as pd
import yaml


def _write_table(path: Path, rows: list[dict]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return str(path)


def build_demo_project_tree() -> dict[str, str]:
    root = Path(tempfile.mkdtemp(prefix="aou-workbench-test-"))
    inputs = root / "inputs"
    configs = root / "configs"
    outputs = root / "outputs"
    configs.mkdir(parents=True, exist_ok=True)
    (configs / "rhabdo").mkdir(parents=True, exist_ok=True)

    cohort_rows = []
    for person_id in range(1, 16):
        cohort_rows.append(
            {
                "person_id": person_id,
                "age": 35 + person_id,
                "sex": "female" if person_id in {1, 3, 4, 5, 6, 7, 8, 13, 14, 15} else "male",
                "observation_start_date": "2020-01-01",
                "observation_end_date": "2023-12-31",
            }
        )
    condition_rows = [
        {"person_id": 1, "condition_concept_id": 100, "condition_start_date": "2022-01-10"},
        {"person_id": 2, "condition_concept_id": 100, "condition_start_date": "2022-02-12"},
        {"person_id": 3, "condition_concept_id": 100, "condition_start_date": "2022-03-15"},
        {"person_id": 4, "condition_concept_id": 100, "condition_start_date": "2022-04-01"},
    ]
    measurement_rows = [
        {"person_id": 1, "measurement_concept_id": 900, "value_as_number": 6800, "measurement_date": "2022-01-11"},
        {"person_id": 2, "measurement_concept_id": 900, "value_as_number": 7200, "measurement_date": "2022-02-13"},
        {"person_id": 3, "measurement_concept_id": 900, "value_as_number": 5600, "measurement_date": "2022-03-16"},
        {"person_id": 4, "measurement_concept_id": 900, "value_as_number": 1800, "measurement_date": "2022-04-01"},
    ]
    ancestry_rows = []
    for person_id in range(1, 16):
        ancestry = "EUR" if person_id in {1, 3, 4, 5, 6, 7, 8, 13, 14, 15} else "AFR"
        row = {"research_id": person_id, "ancestry_pred": ancestry}
        for pc in range(1, 11):
            row[f"PC{pc}"] = round((person_id * pc) / 100.0, 4)
        ancestry_rows.append(row)
    clinical_rows = []
    for person_id in range(1, 16):
        clinical_rows.append(
            {
                "person_id": person_id,
                "statin_exposure": int(person_id in {1, 5, 9}),
                "crush_injury": int(person_id in {2}),
                "sepsis": int(person_id in {3, 10}),
                "renal_injury": int(person_id in {1, 3, 11}),
            }
        )

    stage1_rows = [
        {"person_id": 1, "variant_id": "11-5227002-T-A", "gene": "HBB", "dosage": 1},
        {"person_id": 4, "variant_id": "11-5227002-T-A", "gene": "HBB", "dosage": 1},
        {"person_id": 9, "variant_id": "11-5227002-T-A", "gene": "HBB", "dosage": 1},
        {"person_id": 2, "variant_id": "19-38451842-C-T", "gene": "RYR1", "dosage": 1},
        {"person_id": 3, "variant_id": "19-38499993-G-A", "gene": "RYR1", "dosage": 1},
        {"person_id": 5, "variant_id": "19-38499993-G-A", "gene": "RYR1", "dosage": 1},
        {"person_id": 1, "variant_id": "1-53202427-C-T", "gene": "CPT2", "dosage": 1},
        {"person_id": 3, "variant_id": "1-53202427-C-T", "gene": "CPT2", "dosage": 1},
        {"person_id": 10, "variant_id": "1-201060815-C-T", "gene": "CACNA1S", "dosage": 1},
        {"person_id": 2, "variant_id": "11-64526623-C-T", "gene": "PYGM", "dosage": 1},
        {"person_id": 12, "variant_id": "11-64526623-C-T", "gene": "PYGM", "dosage": 1},
    ]

    stage2_rows = [
        {"person_id": 1, "variant_id": "19-38451842-C-T", "gene": "RYR1", "dosage": 1, "clinvar_significance": "Pathogenic", "consequence": "missense_variant", "revel": 0.95, "max_af": 0.0001},
        {"person_id": 2, "variant_id": "1-53202427-C-T", "gene": "CPT2", "dosage": 1, "clinvar_significance": "Likely pathogenic", "consequence": "missense_variant", "revel": 0.86, "max_af": 0.0002},
        {"person_id": 3, "variant_id": "11-64526623-C-T", "gene": "PYGM", "dosage": 1, "clinvar_significance": "Uncertain significance", "consequence": "frameshift_variant", "revel": 0.10, "max_af": 0.0001},
        {"person_id": 9, "variant_id": "1-53202427-C-T", "gene": "CPT2", "dosage": 1, "clinvar_significance": "Pathogenic", "consequence": "missense_variant", "revel": 0.91, "max_af": 0.0002},
        {"person_id": 6, "variant_id": "19-38577955-G-A", "gene": "RYR1", "dosage": 1, "clinvar_significance": "Benign", "consequence": "missense_variant", "revel": 0.20, "max_af": 0.0002},
        {"person_id": 1, "variant_id": "19-38580094-C-T", "gene": "RYR1", "dosage": 1, "clinvar_significance": "Likely pathogenic", "consequence": "missense_variant", "revel": 0.88, "max_af": 0.0003},
    ]
    stage3_rows = stage2_rows + [
        {"person_id": 2, "variant_id": "19-38457545-C-T", "gene": "RYR1", "dosage": 1, "clinvar_significance": "Pathogenic", "consequence": "missense_variant", "revel": 0.93, "max_af": 0.0002},
        {"person_id": 3, "variant_id": "1-114693436-C-T", "gene": "AMPD1", "dosage": 1, "clinvar_significance": "Likely pathogenic", "consequence": "stop_gained", "revel": 0.22, "max_af": 0.0004},
        {"person_id": 11, "variant_id": "1-114693436-C-T", "gene": "AMPD1", "dosage": 1, "clinvar_significance": "Benign", "consequence": "stop_gained", "revel": 0.18, "max_af": 0.0004},
    ]

    genotype_rows = []
    variants = [
        ("1-100-A-G", "1", 100, "RYR1", [1, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]),
        ("1-200-C-T", "1", 200, "CPT2", [0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0]),
        ("19-300-G-A", "19", 300, "PYGM", [2, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]),
        ("7-400-T-C", "7", 400, "LPIN1", [0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0]),
        ("11-500-A-C", "11", 500, "CKM", [1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0]),
    ]
    for variant_id, chromosome, position, gene, dosages in variants:
        maf = sum(dosages) / (2 * len(dosages))
        for index, dosage in enumerate(dosages, start=1):
            genotype_rows.append(
                {
                    "person_id": index,
                    "variant_id": variant_id,
                    "chromosome": chromosome,
                    "position": position,
                    "gene": gene,
                    "dosage": dosage,
                    "maf": round(maf, 4),
                }
            )
    annotation_rows = [
        {"variant_id": "1-100-A-G", "nearest_gene": "RYR1", "annotation": "intronic"},
        {"variant_id": "1-200-C-T", "nearest_gene": "CPT2", "annotation": "missense"},
        {"variant_id": "19-300-G-A", "nearest_gene": "PYGM", "annotation": "missense"},
        {"variant_id": "7-400-T-C", "nearest_gene": "LPIN1", "annotation": "regulatory"},
        {"variant_id": "11-500-A-C", "nearest_gene": "CKM", "annotation": "synonymous"},
    ]

    paths = {
        "cohort_table": _write_table(inputs / "cohort.tsv", cohort_rows),
        "condition_table": _write_table(inputs / "condition.tsv", condition_rows),
        "measurement_table": _write_table(inputs / "measurement.tsv", measurement_rows),
        "ancestry_table": _write_table(inputs / "ancestry.tsv", ancestry_rows),
        "clinical_table": _write_table(inputs / "clinical.tsv", clinical_rows),
        "stage1_table": _write_table(inputs / "stage1.tsv", stage1_rows),
        "stage2_table": _write_table(inputs / "stage2.tsv", stage2_rows),
        "stage3_table": _write_table(inputs / "stage3.tsv", stage3_rows),
        "stage4_table": _write_table(inputs / "stage4.tsv", genotype_rows),
        "stage4_annotation_table": _write_table(inputs / "stage4_annotations.tsv", annotation_rows),
        "root": str(root),
        "output_dir": str(outputs),
    }
    max_unrelated_path = inputs / "max_unrelated.tsv"
    _write_table(max_unrelated_path, [{"person_id": str(person_id)} for person_id in range(1, 13)])
    paths["max_unrelated"] = str(max_unrelated_path)

    workbench = {
        "storage_root": str(outputs),
        "workspace_bucket": "gs://test-workspace-bucket",
        "workspace_cdr": "test-project.C2024Q3R9",
        "requester_pays_project": "test-project",
    }
    phenotype = {
        "tables": {
            "cohort_table": paths["cohort_table"],
            "condition_table": paths["condition_table"],
            "measurement_table": paths["measurement_table"],
            "ancestry_table": paths["ancestry_table"],
            "clinical_table": paths["clinical_table"],
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
        "clinical_person_id_column": "person_id",
        "clinical_cofactor_columns": ["statin_exposure", "crush_injury", "sepsis", "renal_injury"],
        "definite": {
            "condition_concept_ids": [100],
            "measurement_terms": ["creatine kinase"],
            "measurement_concept_ids": [900],
            "measurement_min": 5000,
            "require_condition": True,
            "require_measurement": True,
        },
        "probable": {
            "condition_concept_ids": [100],
            "require_condition": True,
            "require_measurement": False,
        },
        "min_observation_days": 180,
    }
    cohort = {
        "control_ratio": 2,
        "minimum_controls": 1,
        "age_tolerance_years": 8,
        "index_window_days": 365,
        "exact_match_columns": ["is_female", "ancestry_pred"],
        "primary_case_tier": "definite",
        "covariate_columns": ["age_at_index", "is_female", "pc1", "pc2", "pc3", "statin_exposure"],
    }
    panel = {
        "genes_of_interest": ["RYR1", "CPT2", "PYGM", "AMPD1"],
        "burden_target_genes": ["RYR1", "CPT2", "PYGM", "AMPD1"],
        "a_priori_variants": [
            {"label": "HBB sickle trait", "gene": "HBB", "contig": "chr11", "position": 5227002, "ref": "T", "alt": "A", "rsid": "rs334"},
            {"label": "RYR1 Arg401Cys", "gene": "RYR1", "contig": "chr19", "position": 38451842, "ref": "C", "alt": "T", "rsid": "rs193922764"},
            {"label": "RYR1 Gly2434Arg", "gene": "RYR1", "contig": "chr19", "position": 38499993, "ref": "G", "alt": "A", "rsid": "rs121918593"},
            {"label": "CPT2 Ser113Leu", "gene": "CPT2", "contig": "chr1", "position": 53202427, "ref": "C", "alt": "T", "rsid": "rs74315294"},
            {"label": "CACNA1S Arg1086His", "gene": "CACNA1S", "contig": "chr1", "position": 201060815, "ref": "C", "alt": "T", "rsid": "rs1800559"},
            {"label": "PYGM Arg50Ter", "gene": "PYGM", "contig": "chr11", "position": 64526623, "ref": "C", "alt": "T", "rsid": "rs116987552"},
        ],
    }
    analysis = {
        "analysis_name": "unit_test_rhabdo",
        "output_dir": str(outputs),
        "matched_outcome_column": "analysis_case",
        "stage1": {
            "variant_table": paths["stage1_table"],
            "covariates": ["age_at_index", "is_female", "pc1", "pc2", "pc3", "statin_exposure"],
        },
        "stage2": {
            "variant_table": paths["stage2_table"],
            "covariates": ["age_at_index", "is_female", "pc1", "pc2", "pc3", "statin_exposure"],
        },
        "stage3": {
            "variant_table": paths["stage3_table"],
            "mode": "targeted",
            "covariates": ["age_at_index", "is_female", "pc1", "pc2", "pc3", "statin_exposure"],
        },
        "stage4": {
            "genotype_table": paths["stage4_table"],
            "annotation_table": paths["stage4_annotation_table"],
            "covariates": ["age_at_index", "is_female", "pc1", "pc2", "pc3", "statin_exposure"],
            "min_minor_allele_count": 2,
        },
    }

    config_paths = {
        "workbench": str(configs / "workbench.yaml"),
        "phenotype": str(configs / "rhabdo" / "phenotype.yaml"),
        "cohort": str(configs / "rhabdo" / "cohort.yaml"),
        "panel": str(configs / "rhabdo" / "panel.yaml"),
        "analysis": str(configs / "rhabdo" / "analysis.yaml"),
    }
    for path, payload in (
        (config_paths["workbench"], workbench),
        (config_paths["phenotype"], phenotype),
        (config_paths["cohort"], cohort),
        (config_paths["panel"], panel),
        (config_paths["analysis"], analysis),
    ):
        Path(path).write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")

    stage1_manifest_path = Path(paths["stage1_table"]).with_suffix(".samples.tsv")
    _write_table(stage1_manifest_path, [{"person_id": str(person_id)} for person_id in range(1, 16)])

    paths.update(config_paths)
    return paths
