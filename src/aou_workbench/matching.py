"""Matched-control cohort construction."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import ProjectConfig


def _control_anchor_date(controls: pd.DataFrame) -> pd.Series:
    midpoint = controls["obs_start_date"] + (
        (controls["obs_end_date"] - controls["obs_start_date"]) / 2
    )
    return controls["baseline_index_date"].combine_first(midpoint)


def match_case_controls(cohort_df: pd.DataFrame, config: ProjectConfig) -> pd.DataFrame:
    cases = cohort_df[cohort_df["case_tier"] == config.cohort.primary_case_tier].copy()
    controls = cohort_df[cohort_df["eligible_control"]].copy()
    if cases.empty or controls.empty:
        return pd.DataFrame(columns=list(cohort_df.columns) + ["analysis_case", "match_group_id"])

    controls["control_anchor_date"] = _control_anchor_date(controls)
    used_controls: set[str] = set()
    matched_rows: list[pd.Series] = []

    for case in cases.sort_values(["index_date", "person_id"]).itertuples(index=False):
        pool = controls.copy()
        if not config.cohort.allow_reuse_controls:
            pool = pool[~pool["person_id"].isin(used_controls)]
        for column in config.cohort.exact_match_columns:
            if column in pool.columns:
                pool = pool[pool[column].fillna("__MISSING__") == getattr(case, column)]
        case_index_date = getattr(case, "index_date")
        if pd.notna(case_index_date):
            pool = pool[
                (pool["obs_start_date"].isna() | (pool["obs_start_date"] <= case_index_date))
                & (pool["obs_end_date"].isna() | (pool["obs_end_date"] >= case_index_date))
            ]
            pool["index_distance_days"] = (
                pool["control_anchor_date"] - case_index_date
            ).abs().dt.days.fillna(config.cohort.index_window_days + 1)
            pool = pool[pool["index_distance_days"] <= config.cohort.index_window_days]
        else:
            pool["index_distance_days"] = config.cohort.index_window_days
        pool["age_distance_years"] = (pool["age_at_index"] - getattr(case, "age_at_index")).abs()
        pool = pool[pool["age_distance_years"] <= config.cohort.age_tolerance_years]
        pool = pool.sort_values(["age_distance_years", "index_distance_days", "person_id"])
        chosen = pool.head(config.cohort.control_ratio).copy()
        if len(chosen) < config.cohort.minimum_controls:
            if config.cohort.require_complete_matches:
                continue
            if chosen.empty:
                continue

        match_group_id = f"match-{case.person_id}"
        case_row = pd.Series(case._asdict()).copy()
        case_row["analysis_case"] = 1
        case_row["matched_case_person_id"] = case.person_id
        case_row["match_group_id"] = match_group_id
        case_row["match_role"] = "case"
        case_row["match_rank"] = 0
        case_row["matched_control_count"] = int(len(chosen))
        case_row["match_complete"] = bool(len(chosen) >= config.cohort.control_ratio)
        matched_rows.append(case_row)

        for rank, control in enumerate(chosen.itertuples(index=False), start=1):
            control_row = pd.Series(control._asdict()).copy()
            control_row["analysis_case"] = 0
            control_row["matched_case_person_id"] = case.person_id
            control_row["match_group_id"] = match_group_id
            control_row["match_role"] = "control"
            control_row["match_rank"] = rank
            control_row["matched_control_count"] = int(len(chosen))
            control_row["match_complete"] = bool(len(chosen) >= config.cohort.control_ratio)
            control_row["control_source_index_date"] = control_row.get("baseline_index_date")
            control_row["index_date"] = case_index_date
            control_row["case_tier"] = "control"
            matched_rows.append(control_row)
            if not config.cohort.allow_reuse_controls:
                used_controls.add(str(control.person_id))

    if not matched_rows:
        return pd.DataFrame(columns=list(cohort_df.columns) + ["analysis_case", "match_group_id"])
    matched = pd.DataFrame(matched_rows).sort_values(["match_group_id", "match_rank"]).reset_index(drop=True)
    return matched


def matching_qc_summary(matched_df: pd.DataFrame) -> dict[str, int]:
    if matched_df.empty:
        return {"matched_cases": 0, "matched_controls": 0, "match_groups": 0}
    return {
        "matched_cases": int((matched_df["analysis_case"] == 1).sum()),
        "matched_controls": int((matched_df["analysis_case"] == 0).sum()),
        "match_groups": int(matched_df["match_group_id"].nunique()),
    }


__all__ = ["match_case_controls", "matching_qc_summary"]
