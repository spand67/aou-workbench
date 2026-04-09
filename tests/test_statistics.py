from __future__ import annotations

import unittest
from unittest import mock

import numpy as np
import pandas as pd

from aou_workbench.statistics import run_binary_logistic_regression, summarize_binary_exposure


class StatisticsTests(unittest.TestCase):
    def test_binary_exposure_matches_string_exposures_to_integer_person_ids(self) -> None:
        sample_df = pd.DataFrame(
            [
                {"person_id": 1, "analysis_case": 1},
                {"person_id": 2, "analysis_case": 0},
                {"person_id": 3, "analysis_case": 1},
                {"person_id": 4, "analysis_case": 0},
            ]
        )
        exposure = pd.Series({"1": 1.0, "2": 1.0})

        result = summarize_binary_exposure(
            sample_df,
            exposure,
            outcome_column="analysis_case",
            carrier_threshold=1.0,
        )

        self.assertEqual(result["case_carriers"], 1)
        self.assertEqual(result["control_carriers"], 1)

    def test_logistic_regression_matches_string_exposures_to_integer_person_ids(self) -> None:
        sample_df = pd.DataFrame(
            [
                {"person_id": 1, "analysis_case": 1, "age_at_index": 45.0},
                {"person_id": 2, "analysis_case": 0, "age_at_index": 46.0},
                {"person_id": 3, "analysis_case": 1, "age_at_index": 47.0},
                {"person_id": 4, "analysis_case": 0, "age_at_index": 48.0},
            ]
        )
        exposure = pd.Series({"1": 1.0, "2": 0.0, "3": 1.0, "4": 0.0})

        result = run_binary_logistic_regression(
            sample_df,
            exposure,
            outcome_column="analysis_case",
            covariates=("age_at_index",),
        )

        self.assertEqual(result["n_samples"], 4)

    def test_logistic_regression_returns_nan_when_optimizer_fails(self) -> None:
        sample_df = pd.DataFrame(
            [
                {"person_id": "1", "analysis_case": 1, "age_at_index": 45.0},
                {"person_id": "2", "analysis_case": 0, "age_at_index": 46.0},
                {"person_id": "3", "analysis_case": 1, "age_at_index": 47.0},
                {"person_id": "4", "analysis_case": 0, "age_at_index": 48.0},
            ]
        )
        exposure = pd.Series({"1": 1.0, "2": 0.0, "3": 1.0, "4": 0.0})

        failed_opt = mock.Mock(
            success=False,
            message="did not converge",
            x=np.array([0.0, 1.0]),
            hess_inv=np.eye(2),
        )

        with mock.patch("aou_workbench.statistics.minimize", return_value=failed_opt):
            result = run_binary_logistic_regression(
                sample_df,
                exposure,
                outcome_column="analysis_case",
                covariates=("age_at_index",),
            )

        self.assertTrue(np.isnan(result["beta"]))
        self.assertTrue(np.isnan(result["se"]))
        self.assertTrue(np.isnan(result["odds_ratio"]))
        self.assertTrue(np.isnan(result["regression_p"]))


if __name__ == "__main__":
    unittest.main()
