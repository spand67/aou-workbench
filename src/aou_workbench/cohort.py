"""Cohort construction for tiered rhabdomyolysis analyses."""

from __future__ import annotations

import ast
from typing import Any

import numpy as np
import pandas as pd

from .config import CaseTierRule, ProjectConfig
from .io_utils import parse_date, query_bigquery_dataframe, read_table
from .phenotype_sql import render_baseline_sql, render_case_tier_sql, render_clinical_cofactors_sql


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


def _prepare_baseline_local(config: ProjectConfig) -> pd.DataFrame:
    if not config.phenotype.tables.cohort_table:
        raise ValueError("A local cohort_table is required for non-BigQuery cohort builds.")
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
    base["year_of_birth"] = pd.to_numeric(raw.get("year_of_birth"), errors="coerce")
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


def _prepare_baseline_bigquery(config: ProjectConfig) -> pd.DataFrame:
    frame = query_bigquery_dataframe(render_baseline_sql(config))
    if frame.empty:
        return frame
    frame["obs_start_date"] = parse_date(frame["obs_start_date"])
    frame["obs_end_date"] = parse_date(frame["obs_end_date"])
    if "baseline_index_date" in frame.columns:
        frame["baseline_index_date"] = parse_date(frame["baseline_index_date"])
    if "birth_date" not in frame.columns:
        frame["birth_date"] = pd.NaT
    frame["year_of_birth"] = pd.to_numeric(frame.get("year_of_birth"), errors="coerce")
    frame["age_raw"] = pd.to_numeric(frame.get("age_raw"), errors="coerce")
    frame["is_female"] = pd.to_numeric(frame.get("is_female"), errors="coerce")
    return frame


def _prepare_denominator_local(config: ProjectConfig) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.condition_table).copy()
    if raw.empty:
        return pd.DataFrame(
            columns=[
                "person_id",
                "omop_condition_record_dates",
                "omop_condition_record_source",
                "eligible_ehr_denominator",
            ]
        )
    frame = pd.DataFrame(
        {
            "person_id": raw[config.phenotype.person_id_column].astype(str),
            "condition_date": parse_date(raw[config.phenotype.condition_date_column]),
        }
    )
    source_column = next(
        (
            column
            for column in ("condition_source_vocabulary_id", "source_vocabulary_id", "vocabulary_id")
            if column in raw.columns
        ),
        None,
    )
    if source_column:
        source_values = raw[source_column].astype(str).str.upper()
        use_icd = source_values.str.startswith("ICD").any()
    else:
        source_values = pd.Series("", index=raw.index)
        use_icd = False
    if use_icd:
        frame = frame[source_values.str.startswith("ICD")].copy()
        source = "icd_source_condition_records"
    else:
        source = "all_condition_records"
    counts = (
        frame.dropna(subset=["condition_date"])
        .groupby("person_id", as_index=False)["condition_date"]
        .nunique()
        .rename(columns={"condition_date": "omop_condition_record_dates"})
    )
    counts["omop_condition_record_source"] = source
    counts["eligible_ehr_denominator"] = counts["omop_condition_record_dates"] >= 2
    return counts


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


def _prepare_ancestry_table(config: ProjectConfig) -> pd.DataFrame:
    if not config.phenotype.tables.ancestry_table:
        return pd.DataFrame(columns=["person_id", "ancestry_pred", *[f"pc{i}" for i in range(1, 11)]])
    raw = read_table(config.phenotype.tables.ancestry_table).copy()
    person_column = (
        config.phenotype.ancestry_person_id_column
        if config.phenotype.ancestry_person_id_column in raw.columns
        else config.phenotype.person_id_column
    )
    ancestry = pd.DataFrame({"person_id": raw[person_column].astype(str)})
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


def _aggregate_condition_hits_local(
    config: ProjectConfig,
    concept_ids: tuple[int, ...],
    concept_terms: tuple[str, ...] = (),
) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.condition_table).copy()
    mask = pd.Series(False, index=raw.index)
    if concept_ids:
        mask |= raw[config.phenotype.condition_concept_column].isin(concept_ids)
    if concept_terms and "condition_concept_name" in raw.columns:
        lowered = raw["condition_concept_name"].astype(str).str.lower()
        term_mask = pd.Series(False, index=raw.index)
        for term in concept_terms:
            term_mask |= lowered.str.contains(term.lower(), regex=False)
        mask |= term_mask
    filtered = raw[mask].copy()
    if filtered.empty:
        return pd.DataFrame(columns=["person_id", "condition_date"])
    filtered["person_id"] = filtered[config.phenotype.person_id_column].astype(str)
    filtered["condition_date"] = parse_date(filtered[config.phenotype.condition_date_column])
    return filtered.groupby("person_id", as_index=False)["condition_date"].min()


def _aggregate_measurement_hits_local(config: ProjectConfig, tier: CaseTierRule) -> pd.DataFrame:
    filtered = _aggregate_measurement_events_local(config, tier)
    if filtered.empty:
        return pd.DataFrame(columns=["person_id", "measurement_date", "measurement_value"])
    return (
        filtered.sort_values(["person_id", "measurement_date"])
        .groupby("person_id", as_index=False)
        .agg(measurement_date=("measurement_date", "min"), measurement_value=("measurement_value", "max"))
    )


def _aggregate_measurement_events_local(config: ProjectConfig, tier: CaseTierRule) -> pd.DataFrame:
    raw = read_table(config.phenotype.tables.measurement_table).copy()
    mask = pd.Series(False, index=raw.index)
    if tier.measurement_concept_ids:
        mask |= raw[config.phenotype.measurement_concept_column].isin(tier.measurement_concept_ids)
    if tier.measurement_terms and "measurement_concept_name" in raw.columns:
        lowered = raw["measurement_concept_name"].astype(str).str.lower()
        term_mask = pd.Series(False, index=raw.index)
        for term in tier.measurement_terms:
            term_mask |= lowered.str.contains(term.lower(), regex=False)
        mask |= term_mask
    filtered = raw[mask].copy()
    if tier.measurement_min is not None:
        values = pd.to_numeric(filtered[config.phenotype.measurement_value_column], errors="coerce")
        filtered = filtered[values >= tier.measurement_min].copy()
    if filtered.empty:
        return pd.DataFrame(columns=["person_id", "measurement_date", "measurement_value"])
    filtered["person_id"] = filtered[config.phenotype.person_id_column].astype(str)
    filtered["measurement_date"] = parse_date(filtered[config.phenotype.measurement_date_column])
    filtered["measurement_value"] = pd.to_numeric(filtered[config.phenotype.measurement_value_column], errors="coerce")
    return filtered[["person_id", "measurement_date", "measurement_value"]].copy()


def _case_tier_hits_local(config: ProjectConfig, tier: CaseTierRule) -> pd.DataFrame:
    conditions = _aggregate_condition_hits_local(
        config,
        tier.condition_concept_ids,
        tier.condition_terms,
    )
    measurements = _aggregate_measurement_events_local(config, tier)
    empty_measurement = pd.DataFrame(
        {
            "measurement_date": pd.Series(pd.NaT, dtype="datetime64[ns]"),
            "measurement_value": pd.Series(dtype=float),
        }
    )
    if tier.require_condition and not tier.require_measurement:
        return pd.concat([conditions.reset_index(drop=True), empty_measurement.iloc[: len(conditions)].reset_index(drop=True)], axis=1)
    if tier.require_measurement and not tier.require_condition:
        return _aggregate_measurement_hits_local(config, tier).assign(condition_date=pd.NaT)[
            ["person_id", "condition_date", "measurement_date", "measurement_value"]
        ]
    if not tier.require_condition and not tier.require_measurement:
        condition_hits = pd.concat(
            [conditions.reset_index(drop=True), empty_measurement.iloc[: len(conditions)].reset_index(drop=True)],
            axis=1,
        )
        measurement_hits = _aggregate_measurement_hits_local(config, tier).assign(condition_date=pd.NaT)[
            ["person_id", "condition_date", "measurement_date", "measurement_value"]
        ]
        return pd.concat([condition_hits, measurement_hits], ignore_index=True).drop_duplicates("person_id")
    if conditions.empty or measurements.empty:
        return pd.DataFrame(columns=["person_id", "condition_date", "measurement_date", "measurement_value"])
    merged = conditions.merge(measurements, on="person_id", how="inner")
    delta = (merged["measurement_date"] - merged["condition_date"]).dt.days
    merged = merged[
        delta.between(
            tier.measurement_window_start_days,
            tier.measurement_window_end_days,
            inclusive="both",
        )
    ].copy()
    if merged.empty:
        return pd.DataFrame(columns=["person_id", "condition_date", "measurement_date", "measurement_value"])
    return (
        merged.sort_values(["person_id", "measurement_date"])
        .groupby("person_id", as_index=False)
        .agg(
            condition_date=("condition_date", "first"),
            measurement_date=("measurement_date", "first"),
            measurement_value=("measurement_value", "max"),
        )
    )


def _high_ck_rule(config: ProjectConfig) -> CaseTierRule:
    return CaseTierRule(
        name="high_ck",
        measurement_concept_ids=config.phenotype.definite.measurement_concept_ids,
        measurement_terms=config.phenotype.definite.measurement_terms,
        measurement_min=config.phenotype.definite.measurement_min,
        require_condition=False,
        require_measurement=True,
    )


def _prepare_clinical_local(config: ProjectConfig) -> pd.DataFrame:
    columns = config.phenotype.clinical_cofactor_columns
    if config.phenotype.tables.clinical_table and columns:
        raw = read_table(config.phenotype.tables.clinical_table).copy()
        output = pd.DataFrame({"person_id": raw[config.phenotype.clinical_person_id_column].astype(str)})
        for column in columns:
            output[column] = raw[column] if column in raw.columns else np.nan
        return output
    if not config.phenotype.clinical_cofactors:
        return pd.DataFrame(columns=["person_id"])
    output: pd.DataFrame | None = None
    for rule in config.phenotype.clinical_cofactors:
        hits = _aggregate_condition_hits_local(config, rule.condition_concept_ids, rule.condition_terms)
        column = pd.DataFrame({"person_id": hits["person_id"], rule.name: 1})
        output = column if output is None else output.merge(column, on="person_id", how="outer")
    if output is None:
        return pd.DataFrame(columns=["person_id"])
    for rule in config.phenotype.clinical_cofactors:
        output[rule.name] = output[rule.name].fillna(0).astype(int)
    return output


def _prepare_clinical_bigquery(config: ProjectConfig) -> pd.DataFrame:
    if not config.phenotype.clinical_cofactors:
        return pd.DataFrame(columns=["person_id"])
    frame = query_bigquery_dataframe(render_clinical_cofactors_sql(config))
    if frame.empty:
        return frame
    for rule in config.phenotype.clinical_cofactors:
        if rule.name in frame.columns:
            frame[rule.name] = pd.to_numeric(frame[rule.name], errors="coerce").fillna(0).astype(int)
    return frame


def _merge_case_hits(base: pd.DataFrame, tier_hits: pd.DataFrame, tier: CaseTierRule, prefix: str) -> pd.DataFrame:
    output = base.merge(tier_hits, on="person_id", how="left")
    if "condition_date" in output.columns:
        output["condition_date"] = parse_date(output["condition_date"])
    if "measurement_date" in output.columns:
        output["measurement_date"] = parse_date(output["measurement_date"])
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
        ).dt.days
        rule_hit &= day_delta.between(
            tier.measurement_window_start_days,
            tier.measurement_window_end_days,
            inclusive="both",
        ).fillna(False)
    elif not tier.require_condition and not tier.require_measurement:
        rule_hit &= output[f"{prefix}_has_condition"] | output[f"{prefix}_has_measurement"]
    output[f"{prefix}_rule_hit"] = rule_hit
    output[f"{prefix}_index_date"] = output[f"{prefix}_condition_date"].combine_first(output[f"{prefix}_measurement_date"])
    return output


def _merge_high_ck_hits(base: pd.DataFrame, high_ck_hits: pd.DataFrame) -> pd.DataFrame:
    output = base.merge(high_ck_hits, on="person_id", how="left")
    if "condition_date" in output.columns:
        output = output.drop(columns=["condition_date"])
    output = output.rename(
        columns={
            "measurement_date": "high_ck_measurement_date",
            "measurement_value": "high_ck_measurement_value",
        }
    )
    output["high_ck_has_measurement"] = output["high_ck_measurement_date"].notna()
    return output


def _finalize_cohort(
    cohort_df: pd.DataFrame,
    config: ProjectConfig,
    *,
    control_exclusion_ids: set[str] | None = None,
) -> pd.DataFrame:
    output = cohort_df.copy()
    output["omop_condition_record_dates"] = pd.to_numeric(
        output.get("omop_condition_record_dates", 0),
        errors="coerce",
    ).fillna(0).astype(int)
    output["omop_condition_record_source"] = output.get(
        "omop_condition_record_source",
        pd.Series("all_condition_records", index=output.index),
    ).fillna("all_condition_records")
    output["eligible_ehr_denominator"] = output["omop_condition_record_dates"] >= 2
    output["broad_rule_hit"] = output["broad_rule_hit"] & output["eligible_ehr_denominator"]
    output["definite_rule_hit"] = output["definite_rule_hit"] & output["broad_rule_hit"]
    output["high_ck_has_measurement"] = output.get("high_ck_has_measurement", False).fillna(False).astype(bool)
    output["high_ck_without_rhabdo"] = (
        output["eligible_ehr_denominator"]
        & output["high_ck_has_measurement"]
        & (~output["broad_rule_hit"])
    )
    output["case_tier"] = np.select(
        [
            ~output["eligible_ehr_denominator"],
            output["definite_rule_hit"],
            output["broad_rule_hit"],
            output["high_ck_without_rhabdo"],
        ],
        ["excluded_denominator", "definite", "broad", "indeterminate_ck_only"],
        default="control",
    )
    output["broad_rhabdo_case"] = output["broad_rule_hit"].astype(int)
    output["definite_rhabdo_case"] = output["definite_rule_hit"].astype(int)
    output["rhabdo_case"] = output["broad_rhabdo_case"]
    output["rhabdo_primary_case"] = output["definite_rhabdo_case"]
    output["index_date"] = output["broad_index_date"].combine_first(output["definite_index_date"])
    if "baseline_index_date" in output.columns:
        output["index_date"] = output["index_date"].combine_first(output["baseline_index_date"])
    observation_days = (output["obs_end_date"] - output["obs_start_date"]).dt.days
    output["observation_days"] = observation_days
    output["age_at_index"] = output["age_raw"]
    has_birth = output.get("birth_date", pd.Series(pd.NaT, index=output.index)).notna() & output["index_date"].notna()
    output.loc[has_birth, "age_at_index"] = (
        (output.loc[has_birth, "index_date"] - output.loc[has_birth, "birth_date"]).dt.days / 365.25
    )
    has_yob = output.get("year_of_birth", pd.Series(np.nan, index=output.index)).notna() & output["index_date"].notna()
    output.loc[has_yob, "age_at_index"] = (
        output.loc[has_yob, "index_date"].dt.year - output.loc[has_yob, "year_of_birth"]
    )
    control_exclusion_ids = control_exclusion_ids or set()
    output["control_excluded"] = output["person_id"].isin(control_exclusion_ids)
    output["eligible_control"] = (
        (output["case_tier"] == "control")
        & (~output["control_excluded"])
        & (output["observation_days"].fillna(-1) >= config.phenotype.min_observation_days)
    )
    for rule in config.phenotype.clinical_cofactors:
        if rule.name not in output.columns:
            output[rule.name] = 0
    output = output.sort_values(["case_tier", "index_date", "person_id"]).reset_index(drop=True)
    return output


def _build_local_cohort(config: ProjectConfig) -> pd.DataFrame:
    base = _prepare_baseline_local(config)
    base = base.merge(_prepare_denominator_local(config), on="person_id", how="left")
    base = base.merge(_prepare_ancestry_table(config), on="person_id", how="left")
    base = base.merge(_prepare_clinical_local(config), on="person_id", how="left")
    broad_hits = _case_tier_hits_local(config, config.phenotype.broad)
    broad = _merge_case_hits(base, broad_hits, config.phenotype.broad, "broad")
    definite_hits = _case_tier_hits_local(config, config.phenotype.definite)
    definite = _merge_case_hits(base, definite_hits, config.phenotype.definite, "definite")
    cohort = broad.merge(
        definite[
            [
                "person_id",
                "definite_condition_date",
                "definite_measurement_date",
                "definite_measurement_value",
                "definite_has_condition",
                "definite_has_measurement",
                "definite_rule_hit",
                "definite_index_date",
            ]
        ],
        on="person_id",
        how="left",
    )
    cohort = _merge_high_ck_hits(cohort, _aggregate_measurement_hits_local(config, _high_ck_rule(config)))
    control_exclusions = _aggregate_condition_hits_local(
        config,
        config.phenotype.control_exclusion_concept_ids,
    )
    return _finalize_cohort(cohort, config, control_exclusion_ids=set(control_exclusions["person_id"]))


def _build_bigquery_cohort(config: ProjectConfig) -> pd.DataFrame:
    base = _prepare_baseline_bigquery(config)
    if base.empty:
        return base
    ancestry = _prepare_ancestry_table(config)
    if not ancestry.empty:
        base = base.drop(columns=[column for column in ancestry.columns if column != "person_id" and column in base.columns])
        base = base.merge(ancestry, on="person_id", how="left")
    clinical = _prepare_clinical_bigquery(config)
    if not clinical.empty:
        base = base.merge(clinical, on="person_id", how="left")
    broad = _merge_case_hits(
        base,
        query_bigquery_dataframe(render_case_tier_sql(config, config.phenotype.broad)),
        config.phenotype.broad,
        "broad",
    )
    definite = _merge_case_hits(
        base,
        query_bigquery_dataframe(render_case_tier_sql(config, config.phenotype.definite)),
        config.phenotype.definite,
        "definite",
    )
    cohort = broad.merge(
        definite[
            [
                "person_id",
                "definite_condition_date",
                "definite_measurement_date",
                "definite_measurement_value",
                "definite_has_condition",
                "definite_has_measurement",
                "definite_rule_hit",
                "definite_index_date",
            ]
        ],
        on="person_id",
        how="left",
    )
    cohort = _merge_high_ck_hits(cohort, query_bigquery_dataframe(render_case_tier_sql(config, _high_ck_rule(config))))
    control_exclusion_ids: set[str] = set()
    if config.phenotype.control_exclusion_concept_ids:
        exclusion_rule = CaseTierRule(
            name="control_exclusion",
            condition_concept_ids=config.phenotype.control_exclusion_concept_ids,
            require_condition=True,
            require_measurement=False,
        )
        exclusions = query_bigquery_dataframe(render_case_tier_sql(config, exclusion_rule))
        control_exclusion_ids = set(exclusions.get("person_id", pd.Series(dtype=str)).astype(str))
    return _finalize_cohort(cohort, config, control_exclusion_ids=control_exclusion_ids)


def build_rhabdo_cohort(config: ProjectConfig) -> pd.DataFrame:
    if config.phenotype.tables.cohort_table:
        return _build_local_cohort(config)
    return _build_bigquery_cohort(config)


def cohort_qc_summary(cohort_df: pd.DataFrame) -> dict[str, Any]:
    case_counts = cohort_df["case_tier"].value_counts(dropna=False).to_dict()
    return {
        "n_people": int(len(cohort_df)),
        "eligible_ehr_denominator": int(cohort_df.get("eligible_ehr_denominator", pd.Series(dtype=bool)).sum()),
        "case_counts": {str(key): int(value) for key, value in case_counts.items()},
        "broad_rhabdo_cases": int(cohort_df.get("broad_rhabdo_case", pd.Series(dtype=int)).sum()),
        "definite_rhabdo_cases": int(cohort_df.get("definite_rhabdo_case", pd.Series(dtype=int)).sum()),
        "indeterminate_ck_only": int((cohort_df.get("case_tier", pd.Series(dtype=str)) == "indeterminate_ck_only").sum()),
        "eligible_controls": int(cohort_df.get("eligible_control", pd.Series(dtype=int)).sum()),
    }


__all__ = ["build_rhabdo_cohort", "cohort_qc_summary"]
