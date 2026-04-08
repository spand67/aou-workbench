from __future__ import annotations

import unittest

import pandas as pd

from aou_workbench.config import load_project_config
from aou_workbench.stage1_prepare import _callset_root, _collapse_stage1_rows, _panel_targets_frame, _split_mt_path
from tests.support import build_demo_project_tree


class Stage1PrepareTests(unittest.TestCase):
    def test_split_mt_path_normalizes_known_layouts(self) -> None:
        self.assertEqual(
            _split_mt_path("gs://bucket/v8/wgs/short_read/snpindel/exome/multiMT/hail.mt"),
            "gs://bucket/v8/wgs/short_read/snpindel/exome/splitMT/hail.mt",
        )
        self.assertEqual(
            _split_mt_path("gs://bucket/v8/wgs/short_read/snpindel/acaf_threshold/hail.mt"),
            "gs://bucket/v8/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt",
        )
        self.assertEqual(
            _callset_root("gs://bucket/v8/wgs/short_read/snpindel/clinvar/splitMT/hail.mt"),
            "gs://bucket/v8/wgs/short_read/snpindel/clinvar",
        )

    def test_panel_targets_frame_uses_project_variant_ids(self) -> None:
        paths = build_demo_project_tree()
        config = load_project_config(
            workbench_path=paths["workbench"],
            phenotype_path=paths["phenotype"],
            cohort_path=paths["cohort"],
            panel_path=paths["panel"],
            analysis_path=paths["analysis"],
        )
        frame = _panel_targets_frame(config.panel.a_priori_variants)
        self.assertIn("11-5227002-T-A", set(frame["variant_id"]))
        self.assertIn("19-38451842-C-T", set(frame["variant_id"]))

    def test_collapse_stage1_rows_deduplicates_across_callsets(self) -> None:
        frame = pd.DataFrame(
            [
                {
                    "person_id": "101",
                    "variant_id": "11-5227002-T-A",
                    "gene": "HBB",
                    "label": "HBB sickle trait",
                    "rsid": "rs334",
                    "source": "Deuster2013",
                    "evidence_tier": "candidate_common",
                    "exact_test_model": "carrier_vs_noncarrier",
                    "dosage": 1.0,
                    "callset": "acaf",
                },
                {
                    "person_id": "101",
                    "variant_id": "11-5227002-T-A",
                    "gene": "HBB",
                    "label": "HBB sickle trait",
                    "rsid": "rs334",
                    "source": "Deuster2013",
                    "evidence_tier": "candidate_common",
                    "exact_test_model": "carrier_vs_noncarrier",
                    "dosage": 2.0,
                    "callset": "clinvar",
                },
            ]
        )
        collapsed = _collapse_stage1_rows(frame)
        self.assertEqual(len(collapsed), 1)
        self.assertEqual(collapsed.loc[0, "dosage"], 2.0)
        self.assertEqual(collapsed.loc[0, "callset"], "acaf,clinvar")


if __name__ == "__main__":
    unittest.main()
