"""Cohort construction for tiered rhabdomyolysis analyses."""

from __future__ import annotations

import ast
from typing import Any

import numpy as np
import pandas as pd

from .config import CaseTierRule, ProjectConfig
from .io_utils import parse_date, read_table


def _normalize_sex(series: pd.Series) -> pd.Series:
    lowered = series.astype(str).str.lower().str.strip()
    return lowered.map(
        {
            "female": 1.0,
            "f": 1.0,
            "2": 1.0,
            "male": 0.0,
            "m": 0.0,
            "1": 0.0,
        }
    )


def _prepare_baseline(config: ProjectConfig) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.cohort_table).copy()
    base = pd.DataFrame({"person_id": raw[config.phenotype.person_id_column].astype(str)})
    base["obs_start_date"] = parse_date(raw[config.phenotype.observation_start_column])
    base["obs_end_date"] = parse_date(raw[config.phenotype.observation_end_column])
    if config.phenotype.cohort_index_date_column and config.phenotype.cohort_index_date_column in raw.columns:
        base["baseline_index_date"] = parse_date(raw[config.phenotype.cohort_index_date_column])
    else:
        base["baseline_index_date"] = pd.NaT
    if config.phenotype.birth_date_column and config.phenotype.birth_date_column in raw.columns:
        base["birth_date"] = parse_date(raw[config.phenotype.birth_date_column])
    else:
        base["birth_date"] = pd.NaT
    base["age_raw"] = (
        pd.to_numeric(raw[config.phenotype.age_column], errors="coerce")
        if config.phenotype.age_column in raw.columns
        else np.nan
    )
    base["is_female"] = (
        _normalize_sex(raw[config.phenotype.sex_column])
        if config.phenotype.sex_column in raw.columns
        else np.nan
    )
    return base


def _parse_pca_features(value: Any, count: int) -> list[float]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return [np.nan] * count
    if isinstance(value, list):
        parsed = value
    else:
        text = str(value).strip()
        if not text:
            return [np.nan] * count
        if text.startswith("[") and text.endswith("]"):
            try:
                parsed = ast.literal_eval(text)
            except Exception:
                parsed = [item.strip() for item in text[1:-1].split(",")]
        else:
            parsed = [item.strip() for item in text.split(",")]
    out: list[float] = []
    for item in parsed[:count]:
        try:
            out.append(float(item))
        except (TypeError, ValueError):
            out.append(np.nan)
    while len(out) < count:
        out.append(np.nan)
    return out


def _prepare_ancestry(config: ProjectConfig) -> pd.DataFrame:
    if not config.phenotype.tables.ancestry_table:
        return pd.DataFrame(columns=["person_id", "ancestry_pred", *[f"pc{i}" for i in range(1, 11)]])
    raw = read_table(config.phenotype.tables.ancestry_table).copy()
    ancestry = pd.DataFrame({"person_id": raw[config.phenotype.ancestry_person_id_column].astype(str)})
    ancestry["ancestry_pred"] = raw[config.phenotype.ancestry_label_column].astype(str)
    pc_targets = [f"pc{i}" for i in range(1, len(config.phenotype.pc_columns) + 1)]
    missing_pc_columns = [column for column in config.phenotype.pc_columns if column not in raw.columns]
    if not missing_pc_columns:
        for source, target in zip(config.phenotype.pc_columns, pc_targets):
            ancestry[target] = pd.to_numeric(raw[source], errors="coerce")
        return ancestry
    pca_features = raw.get(config.phenotype.ancestry_pca_features_column, pd.Series([None] * len(raw)))
    parsed = pca_features.apply(lambda value: _parse_pca_features(value, len(pc_targets)))
    for index, target in enumerate(pc_targets):
        ancestry[target] = parsed.apply(lambda values: values[index] if index < len(values) else np.nan)
    return ancestry


def _prepare_clinical(config: ProjectConfig) -> pd.DataFrame:
    columns = config.phenotype.clinical_cofactor_columns
    if not config.phenotype.tables.clinical_table or not columns:
        return pd.DataFrame(columns=["person_id", *columns])
    raw = read_table(config.phenotype.tables.clinical_table).copy()
    output = pd.DataFrame({"person_id": raw[config.phenotype.clinical_person_id_column].astype(str)})
    for column in columns:
        output[column] = raw[column] if column in raw.columns else np.nan
    return output


def _aggregate_condition_hits(config: ProjectConfig, concept_ids: tuple[int, ...]) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.condition_table).copy()
    if not concept_ids:
        return pd.DataFrame(columns=["person_id", "condition_date"])
    filtered = raw[raw[config.phenotype.condition_concept_column].isin(concept_ids)].copy()
    if filtered.empty:
        return pd.DataFrame(columns=["person_id", "condition_date"])
    filtered["person_id"] = filtered[config.phenotype.person_id_column].astype(str)
    filtered["condition_date"] = parse_date(filtered[config.phenotype.condition_date_column])
    return filtered.groupby("person_id", as_index=False)["condition_date"].min()


def _aggregate_measurement_hits(config: ProjectConfig, tier: CaseTierRule) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.measurement_table).copy()
    if not tier.measurement_concept_ids:
        return pd.DataFrame(columns=["person_id", "measurement_date", "measurement_value"])
    filtered = raw[raw[config.phenotype.measurement_concept_column].isin(tier.measurement_concept_ids)].copy()
    if tier.measurement_min is not None:
        values = pd.to_numeric(filtered[config.phenotype.measurement_value_column], errors="coerce")
        filtered = filtered[values >= tier.measurement_min].copy()
    if filtered.empty:
        return pd.DataFrame(columns=["person_id", "measurement_date", "measurement_value"])
    filtered["person_id"] = filtered[config.phenotype.person_id_column].astype(str)
    filtered["measurement_date"] = parse_date(filtered[config.phenotype.measurement_date_column])
    filtered["measurement_value"] = pd.to_numeric(filtered[config.phenotype.measurement_value_column], errors="coerce")
    return (
        filtered.sort_values(["person_id", "measurement_date"])
        .groupby("person_id", as_index=False)
        .agg(measurement_date=("measurement_date", "min"), measurement_value=("measurement_value", "max"))
    )


def _apply_case_rule(base: pd.DataFrame, config: ProjectConfig, tier: CaseTierRule, prefix: str) -> pd.DataFrame:
    output = base.merge(_aggregate_condition_hits(config, tier.condition_concept_ids), on="person_id", how="left")
    output = output.merge(_aggregate_measurement_hits(config, tier), on="person_id", how="left")
    output = output.rename(
        columns={
            "condition_date": f"{prefix}_condition_date",
            "measurement_date": f"{prefix}_measurement_date",
            "measurement_value": f"{prefix}_measurement_value",
        }
    )
    output[f"{prefix}_has_condition"] = output[f"{prefix}_condition_date"].notna()
    output[f"{prefix}_has_measurement"] = output[f"{prefix}_measurement_date"].notna()
    rule_hit = pd.Series(True, index=output.index)
    if tier.require_condition:
        rule_hit &= output[f"{prefix}_has_condition"]
    if tier.require_measurement:
        rule_hit &= output[f"{prefix}_has_measurement"]
    if tier.require_condition and tier.require_measurement:
        day_delta = (
            output[f"{prefix}_measurement_date"] - output[f"{prefix}_condition_date"]
        ).abs().dt.days
        rule_hit &= day_delta.fillna(np.inf) <= tier.lookback_days
    elif not tier.require_condition and not tier.require_measurement:
        rule_hit &= output[f"{prefix}_has_condition"] | output[f"{prefix}_has_measurement"]
    output[f"{prefix}_rule_hit"] = rule_hit
    output[f"{prefix}_index_date"] = output[
        [f"{prefix}_condition_date", f"{prefix}_measurement_date"]
    ].min(axis=1)
    return output


def build_rhabdo_cohort(config: ProjectConfig) -> pd.DataFrame:
    base = _prepare_baseline(config)
    base = base.merge(_prepare_ancestry(config), on="person_id", how="left")
    base = base.merge(_prepare_clinical(config), on="person_id", how="left")

    definite = _apply_case_rule(base, config, config.phenotype.definite, "definite")
    probable = _apply_case_rule(definite, config, config.phenotype.probable, "probable")
    probable["case_tier"] = np.where(
        probable["definite_rule_hit"],
        "definite",
        np.where(probable["probable_rule_hit"], "probable", "control"),
    )
    probable["rhabdo_case"] = (probable["case_tier"] != "control").astype(int)
    probable["rhabdo_primary_case"] = (probable["case_tier"] == "definite").astype(int)
    probable["index_date"] = probable["definite_index_date"].combine_first(probable["probable_index_date"])
    probable["index_date"] = probable["index_date"].combine_first(probable["baseline_index_date"])
    observation_days = (probable["obs_end_date"] - probable["obs_start_date"]).dt.days
    probable["observation_days"] = observation_days
    probable["age_at_index"] = probable["age_raw"]
    has_birth = probable["birth_date"].notna() & probable["index_date"].notna()
    probable.loc[has_birth, "age_at_index"] = (
        (probable.loc[has_birth, "index_date"] - probable.loc[has_birth, "birth_date"]).dt.days / 365.25
    )

    control_exclusions = _aggregate_condition_hits(config, config.phenotype.control_exclusion_concept_ids)
    probable["control_excluded"] = probable["person_id"].isin(control_exclusions["person_id"])
    probable["eligible_control"] = (
        (probable["case_tier"] == "control")
        & (~probable["control_excluded"])
        & (probable["observation_days"].fillna(-1) >= config.phenotype.min_observation_days)
    )
    probable = probable.sort_values(["case_tier", "index_date", "person_id"]).reset_index(drop=True)
    return probable


def cohort_qc_summary(cohort_df: pd.DataFrame) -> dict[str, Any]:
    case_counts = cohort_df["case_tier"].value_counts(dropna=False).to_dict()
    return {
        "n_people": int(len(cohort_df)),
        "case_counts": {str(key): int(value) for key, value in case_counts.items()},
        "eligible_controls": int(cohort_df.get("eligible_control", pd.Series(dtype=int)).sum()),
    }


__all__ = ["build_rhabdo_cohort", "cohort_qc_summary"]
