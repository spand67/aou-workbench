"""Shared statistics helpers used by stage analyses."""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2, fisher_exact, norm


def bh_fdr(pvals: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([])
    order = np.argsort(np.where(np.isnan(p), np.inf, p))
    ranked = p[order]
    q = np.full(n, np.nan, dtype=float)
    prev = 1.0
    for index in range(n - 1, -1, -1):
        rank = index + 1
        current = ranked[index]
        if np.isnan(current):
            value = np.nan
        else:
            value = min(prev, current * n / rank)
            prev = value
        q[index] = value
    output = np.full(n, np.nan, dtype=float)
    output[order] = q
    return output


def resolve_covariates(df: pd.DataFrame, requested: Sequence[str]) -> list[str]:
    return [column for column in requested if column in df.columns]


def _design_matrix(
    sample_df: pd.DataFrame,
    exposure_by_person: pd.Series,
    *,
    outcome_column: str,
    covariates: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    subset = sample_df[["person_id", outcome_column, *covariates]].copy()
    subset = subset.dropna(subset=[outcome_column])
    subset["exposure"] = subset["person_id"].map(exposure_by_person).fillna(0.0)
    subset = subset.dropna()
    y = subset[outcome_column].astype(float).to_numpy()
    columns = [np.ones(len(subset), dtype=float), subset["exposure"].astype(float).to_numpy()]
    for covariate in covariates:
        columns.append(subset[covariate].astype(float).to_numpy())
    x = np.column_stack(columns)
    return x, y, subset


def _logistic_nll(beta: np.ndarray, x: np.ndarray, y: np.ndarray) -> float:
    eta = x @ beta
    p = 1.0 / (1.0 + np.exp(-eta))
    eps = 1e-9
    return -np.sum(y * np.log(p + eps) + (1.0 - y) * np.log(1.0 - p + eps))


def summarize_binary_exposure(
    sample_df: pd.DataFrame,
    exposure_by_person: pd.Series,
    *,
    outcome_column: str,
    carrier_threshold: float = 1.0,
) -> dict[str, float]:
    analysis = sample_df[["person_id", outcome_column]].copy()
    analysis = analysis.dropna(subset=[outcome_column])
    analysis["carrier"] = analysis["person_id"].map(exposure_by_person).fillna(0.0) >= carrier_threshold
    cases = analysis[analysis[outcome_column] == 1]
    controls = analysis[analysis[outcome_column] == 0]
    case_carriers = int(cases["carrier"].sum())
    control_carriers = int(controls["carrier"].sum())
    case_noncarriers = int(len(cases) - case_carriers)
    control_noncarriers = int(len(controls) - control_carriers)
    odds_ratio, p_value = fisher_exact(
        [[case_carriers, case_noncarriers], [control_carriers, control_noncarriers]],
        alternative="two-sided",
    )
    return {
        "n_cases": int(len(cases)),
        "n_controls": int(len(controls)),
        "case_carriers": case_carriers,
        "control_carriers": control_carriers,
        "case_noncarriers": case_noncarriers,
        "control_noncarriers": control_noncarriers,
        "case_carrier_fraction": case_carriers / len(cases) if len(cases) else np.nan,
        "control_carrier_fraction": control_carriers / len(controls) if len(controls) else np.nan,
        "fisher_odds_ratio": odds_ratio,
        "fisher_p": p_value,
    }


def run_binary_logistic_regression(
    sample_df: pd.DataFrame,
    exposure_by_person: pd.Series,
    *,
    outcome_column: str,
    covariates: Sequence[str],
) -> dict[str, float]:
    covariates = resolve_covariates(sample_df, covariates)
    x, y, used = _design_matrix(
        sample_df,
        exposure_by_person,
        outcome_column=outcome_column,
        covariates=covariates,
    )
    if len(used) <= x.shape[1] or len(np.unique(y)) < 2:
        return {
            "beta": np.nan,
            "se": np.nan,
            "odds_ratio": np.nan,
            "regression_p": np.nan,
            "n_samples": int(len(used)),
            "covariates": ";".join(covariates),
        }
    try:
        opt = minimize(_logistic_nll, np.zeros(x.shape[1]), args=(x, y), method="BFGS")
        if not opt.success:
            raise RuntimeError(opt.message)
        beta = opt.x
        cov = np.asarray(opt.hess_inv.todense() if hasattr(opt.hess_inv, "todense") else opt.hess_inv)
        se = math.sqrt(max(cov[1, 1], 0.0))
        z = beta[1] / se if se > 0 else np.nan
        p_value = 2 * (1 - norm.cdf(abs(z))) if np.isfinite(z) else np.nan
        return {
            "beta": float(beta[1]),
            "se": float(se),
            "odds_ratio": float(math.exp(beta[1])) if np.isfinite(beta[1]) else np.nan,
            "regression_p": float(p_value) if np.isfinite(p_value) else np.nan,
            "n_samples": int(len(used)),
            "covariates": ";".join(covariates),
        }
    except Exception:
        return {
            "beta": np.nan,
            "se": np.nan,
            "odds_ratio": np.nan,
            "regression_p": np.nan,
            "n_samples": int(len(used)),
            "covariates": ";".join(covariates),
        }


def genomic_control_lambda(p_values: Sequence[float]) -> float:
    valid = [float(value) for value in p_values if pd.notna(value) and 0 < float(value) <= 1]
    if not valid:
        return float("nan")
    chisq = chi2.isf(valid, df=1)
    median = float(np.median(chisq))
    return median / 0.4549364 if median >= 0 else float("nan")


__all__ = [
    "bh_fdr",
    "genomic_control_lambda",
    "resolve_covariates",
    "run_binary_logistic_regression",
    "summarize_binary_exposure",
]
