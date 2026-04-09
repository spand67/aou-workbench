from __future__ import annotations

import unittest

import pandas as pd

from aou_workbench.config import load_project_config
from aou_workbench.stage1_prepare import (
    _analysis_person_ids,
    _collapse_stage1_rows,
    _matched_person_ids,
    _panel_targets_frame,
    _target_interval_strings,
)
from tests.support import build_demo_project_tree


class Stage1PrepareTests(unittest.TestCase):
    def setUp(self) -> None:
        paths = build_demo_project_tree()
        self.config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )

    def test_panel_targets_frame_preserves_variant_metadata(self) -> None:
        frame = _panel_targets_frame(self.config.panel.a_priori_variants)

        self.assertEqual(frame.shape[0], 6)
        self.assertEqual(
            frame["variant_id"].tolist(),
            [
                "11-5227002-T-A",
                "19-38451842-C-T",
                "19-38499993-G-A",
                "1-53202427-C-T",
                "1-201060815-C-T",
                "11-64526623-C-T",
            ],
        )
        self.assertTrue((frame["exact_test_model"] == "carrier_vs_noncarrier").all())

    def test_analysis_person_ids_sorts_and_deduplicates(self) -> None:
        matched_df = pd.DataFrame({"person_id": [5, "2", 5, 1, "10"]})

        self.assertEqual(_analysis_person_ids(matched_df), ["1", "10", "2", "5"])
        self.assertEqual(_matched_person_ids(matched_df), ["1", "10", "2", "5"])

    def test_target_interval_strings_pad_and_sort_targets(self) -> None:
        frame = _panel_targets_frame(self.config.panel.a_priori_variants)

        self.assertEqual(
            _target_interval_strings(frame, pad_bp=5),
            [
                "chr11:5226997-5227007",
                "chr11:64526618-64526628",
                "chr19:38451837-38451847",
                "chr19:38499988-38499998",
                "chr1:201060810-201060820",
                "chr1:53202422-53202432",
            ],
        )

    def test_collapse_stage1_rows_keeps_max_dosage_and_merges_callsets(self) -> None:
        frame = pd.DataFrame(
            [
                {
                    "person_id": "1",
                    "variant_id": "11-5227002-T-A",
                    "gene": "HBB",
                    "label": "HBB sickle trait",
                    "rsid": "rs334",
                    "source": "literature",
                    "evidence_tier": "candidate_common",
                    "exact_test_model": "carrier_vs_noncarrier",
                    "dosage": 1.0,
                    "callset": "wgs_vds",
                },
                {
                    "person_id": "1",
                    "variant_id": "11-5227002-T-A",
                    "gene": "HBB",
                    "label": "HBB sickle trait",
                    "rsid": "rs334",
                    "source": "literature",
                    "evidence_tier": "candidate_common",
                    "exact_test_model": "carrier_vs_noncarrier",
                    "dosage": 2.0,
                    "callset": "legacy_exact",
                },
                {
                    "person_id": "2",
                    "variant_id": "19-38451842-C-T",
                    "gene": "RYR1",
                    "label": "RYR1 Arg401Cys",
                    "rsid": "rs193922764",
                    "source": "a_priori",
                    "evidence_tier": "core_pathogenic",
                    "exact_test_model": "carrier_vs_noncarrier",
                    "dosage": 1.0,
                    "callset": "wgs_vds",
                },
            ]
        )

        collapsed = _collapse_stage1_rows(frame)

        self.assertEqual(collapsed.shape[0], 2)
        first = collapsed.loc[collapsed["person_id"] == "1"].iloc[0]
        self.assertEqual(first["dosage"], 2.0)
        self.assertEqual(first["callset"], "legacy_exact,wgs_vds")
        self.assertEqual(first["variant_id"], "11-5227002-T-A")


if __name__ == "__main__":
    unittest.main()
