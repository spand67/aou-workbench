from __future__ import annotations

import unittest

import pandas as pd

from aou_workbench.stage1_prior_variants import _variant_exposure


class Stage1PriorVariantTests(unittest.TestCase):
    def test_variant_exposure_uses_homozygous_alt_model_when_requested(self) -> None:
        subset = pd.DataFrame(
            [
                {"person_id": "1", "dosage": 1.0},
                {"person_id": "2", "dosage": 2.0},
                {"person_id": "3", "dosage": 0.0},
            ]
        )

        exposure, threshold = _variant_exposure(subset, exact_test_model="hom_alt_vs_rest")

        self.assertEqual(threshold, 1.0)
        self.assertEqual(exposure.to_dict(), {"1": 0.0, "2": 1.0, "3": 0.0})

    def test_variant_exposure_preserves_dosage_for_carrier_model(self) -> None:
        subset = pd.DataFrame(
            [
                {"person_id": "1", "dosage": 1.0},
                {"person_id": "1", "dosage": 2.0},
                {"person_id": "2", "dosage": 1.0},
            ]
        )

        exposure, threshold = _variant_exposure(subset, exact_test_model="carrier_vs_noncarrier")

        self.assertEqual(threshold, 1.0)
        self.assertEqual(exposure.to_dict(), {"1": 2.0, "2": 1.0})


if __name__ == "__main__":
    unittest.main()
