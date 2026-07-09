from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from aou_workbench.cohort import build_rhabdo_cohort
from aou_workbench.config import load_project_config
from aou_workbench.paths import build_output_paths
from aou_workbench.stage2_plp_panel import run_stage2_plp_panel
from aou_workbench.stage2_prepare import (
    _collapse_vat_annotations,
    _stage2_candidate_cache_path,
    stage2_sample_manifest_path,
)
from tests.support import build_demo_project_tree


class Stage2PlpPanelTests(unittest.TestCase):
    def test_stage2_candidate_cache_path_uses_workspace_bucket_and_mask_settings(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )

        strict = _stage2_candidate_cache_path(
            config,
            variant_table_path=paths["stage2_table"],
            genes=["CPT2", "RYR1"],
            max_af=0.001,
            revel_min=0.8,
            plof_terms=config.analysis.stage2.plof_terms,
            clinvar_plp_terms=config.analysis.stage2.clinvar_plp_terms,
        )
        loose = _stage2_candidate_cache_path(
            config,
            variant_table_path=paths["stage2_table"],
            genes=["CPT2", "RYR1"],
            max_af=0.01,
            revel_min=0.8,
            plof_terms=config.analysis.stage2.plof_terms,
            clinvar_plp_terms=config.analysis.stage2.clinvar_plp_terms,
        )

        self.assertTrue(strict.startswith("gs://test-workspace-bucket/aou-workbench/unit-test-rhabdo/stage2/"))
        self.assertTrue(strict.endswith(".ht"))
        self.assertNotEqual(strict, loose)

    def test_collapse_vat_annotations_combines_transcripts_and_preserves_max_scores(self) -> None:
        vat = pd.DataFrame(
            [
                {
                    "vid": "19-38451842-C-T",
                    "contig": "chr19",
                    "position": "38451842",
                    "ref_allele": "C",
                    "alt_allele": "T",
                    "gene_symbol": "RYR1",
                    "consequence": "missense_variant",
                    "revel": "0.82",
                    "gnomad_max_af": "0.0002",
                    "gvs_max_af": "0.0003",
                    "gvs_all_af": "0.0001",
                    "clinvar_classification": "Pathogenic",
                    "aa_change": "p.Arg401Cys",
                    "dbsnp_rsid": "rs193922764",
                },
                {
                    "vid": "19-38451842-C-T",
                    "contig": "chr19",
                    "position": "38451842",
                    "ref_allele": "C",
                    "alt_allele": "T",
                    "gene_symbol": "RYR1",
                    "consequence": "splice_region_variant",
                    "revel": "0.91",
                    "gnomad_max_af": "0.0001",
                    "gvs_max_af": "0.0004",
                    "gvs_all_af": "",
                    "clinvar_classification": "Likely_pathogenic",
                    "aa_change": "p.Arg401Cys",
                    "dbsnp_rsid": "rs193922764",
                },
                {
                    "vid": "",
                    "contig": "chr2",
                    "position": "27765084",
                    "ref_allele": "G",
                    "alt_allele": "A",
                    "gene_symbol": "CPT2",
                    "consequence": "stop_gained",
                    "revel": "",
                    "gnomad_max_af": "0.00005",
                    "gvs_max_af": "",
                    "gvs_all_af": "",
                    "clinvar_classification": "",
                    "aa_change": "",
                    "dbsnp_rsid": "",
                },
            ]
        )

        collapsed = _collapse_vat_annotations(vat)

        self.assertEqual(collapsed.shape[0], 2)
        ryr1 = collapsed.loc[collapsed["gene"] == "RYR1"].iloc[0]
        self.assertEqual(ryr1["variant_id"], "19-38451842-C-T")
        self.assertEqual(ryr1["contig"], "chr19")
        self.assertEqual(ryr1["position"], 38451842)
        self.assertEqual(ryr1["revel"], 0.91)
        self.assertEqual(ryr1["max_af"], 0.0004)
        self.assertEqual(ryr1["clinvar_significance"], "Likely_pathogenic;Pathogenic")
        self.assertEqual(ryr1["consequence"], "missense_variant;splice_region_variant")

        cpt2 = collapsed.loc[collapsed["gene"] == "CPT2"].iloc[0]
        self.assertEqual(cpt2["variant_id"], "2-27765084-G-A")

    def test_stage2_runs_case_vs_rest_comparisons_using_manifest_filtered_samples(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        cohort_df = build_rhabdo_cohort(config)
        pd.DataFrame({"person_id": cohort_df["person_id"].astype(str)}).to_csv(
            stage2_sample_manifest_path(paths["stage2_table"]),
            sep="\t",
            index=False,
        )

        variant_df, gene_df, person_df = run_stage2_plp_panel(config, cohort_df, build_output_paths(config))

        self.assertFalse(variant_df.empty)
        self.assertFalse(gene_df.empty)
        self.assertFalse(person_df.empty)
        self.assertIn("omop_rhabdo_vs_non_rhabdo", set(gene_df["comparison"]))
        self.assertIn("omop_rhabdo_plus_ck_vs_non_rhabdo", set(gene_df["comparison"]))
        self.assertTrue(set(gene_df["gene"]).issubset(set(config.panel.genes_of_interest)))

    def test_stage2_requires_sample_manifest(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        cohort_df = build_rhabdo_cohort(config)
        manifest_path = Path(stage2_sample_manifest_path(paths["stage2_table"]))
        if manifest_path.exists():
            manifest_path.unlink()

        with self.assertRaises(RuntimeError):
            run_stage2_plp_panel(config, cohort_df, build_output_paths(config))


if __name__ == "__main__":
    unittest.main()
