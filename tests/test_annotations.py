from __future__ import annotations

import unittest

import pandas as pd

from aou_workbench.annotations import annotate_variant_masks
from aou_workbench.config import DEFAULT_CLINVAR_PLP_TERMS, DEFAULT_PLOF_TERMS


class AnnotationMaskTests(unittest.TestCase):
    def test_splice_region_variant_is_not_counted_as_plof_without_other_evidence(self) -> None:
        variants = pd.DataFrame(
            [
                {
                    "variant_id": "splice-region-only",
                    "clinvar_significance": "",
                    "consequence": "splice_region_variant",
                    "revel": 0.2,
                    "max_af": 0.0001,
                },
                {
                    "variant_id": "splice-donor",
                    "clinvar_significance": "",
                    "consequence": "splice_donor_variant",
                    "revel": 0.2,
                    "max_af": 0.0001,
                },
                {
                    "variant_id": "splice-acceptor",
                    "clinvar_significance": "",
                    "consequence": "splice_acceptor_variant",
                    "revel": 0.2,
                    "max_af": 0.0001,
                },
                {
                    "variant_id": "clinvar-plp-splice-region",
                    "clinvar_significance": "Pathogenic",
                    "consequence": "splice_region_variant",
                    "revel": 0.2,
                    "max_af": 0.0001,
                },
            ]
        )

        annotated = annotate_variant_masks(
            variants,
            clinvar_column="clinvar_significance",
            consequence_column="consequence",
            revel_column="revel",
            af_column="max_af",
            max_af=0.01,
            revel_min=0.8,
            plof_terms=DEFAULT_PLOF_TERMS,
            clinvar_plp_terms=DEFAULT_CLINVAR_PLP_TERMS,
        ).set_index("variant_id")

        self.assertFalse(bool(annotated.loc["splice-region-only", "is_plof"]))
        self.assertFalse(bool(annotated.loc["splice-region-only", "mask_primary"]))
        self.assertTrue(bool(annotated.loc["splice-donor", "is_plof"]))
        self.assertTrue(bool(annotated.loc["splice-acceptor", "is_plof"]))
        self.assertTrue(bool(annotated.loc["clinvar-plp-splice-region", "mask_primary"]))


if __name__ == "__main__":
    unittest.main()
