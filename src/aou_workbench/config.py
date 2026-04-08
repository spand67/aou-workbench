"""Configuration models and YAML loading for the AoU workbench."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Iterable, Mapping, Sequence

from .io_utils import load_yaml, slugify, stable_hash

DEFAULT_GENOMICS_BUCKET = "gs://vwb-aou-datasets-controlled"
DEFAULT_WORKSPACE_CDR = "fc-aou-cdr-prod.C2024Q3R9"
DEFAULT_WGS_VDS_PATH = f"{DEFAULT_GENOMICS_BUCKET}/v8/wgs/short_read/snpindel/vds/hail.vds"
DEFAULT_VAT_PATH = (
    f"{DEFAULT_GENOMICS_BUCKET}/v8/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz"
)
DEFAULT_EXOME_MT_PATH = (
    f"{DEFAULT_GENOMICS_BUCKET}/v8/wgs/short_read/snpindel/exome/hail.mt"
)
DEFAULT_CLINVAR_MT_PATH = (
    f"{DEFAULT_GENOMICS_BUCKET}/v8/wgs/short_read/snpindel/clinvar/hail.mt"
)
DEFAULT_ACAF_MT_PATH = (
    f"{DEFAULT_GENOMICS_BUCKET}/v8/wgs/short_read/snpindel/acaf_threshold/hail.mt"
)
DEFAULT_COVARIATES = (
    "age_at_index",
    "is_female",
    "pc1",
    "pc2",
    "pc3",
    "pc4",
    "pc5",
)
DEFAULT_PLOF_TERMS = (
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "stop_lost",
    "transcript_ablation",
)
DEFAULT_CLINVAR_PLP_TERMS = (
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
)
DEFAULT_RHABDO_GENES = (
    "RYR1",
    "CPT2",
    "PYGM",
    "LPIN1",
    "CACNA1S",
    "AMPD1",
    "CKM",
)


def _as_tuple(values: Sequence[Any] | None) -> tuple[Any, ...]:
    return tuple(values or ())


def _as_string_tuple(values: Sequence[Any] | None) -> tuple[str, ...]:
    return tuple(str(value) for value in (values or ()))


def _as_int_tuple(values: Sequence[Any] | None) -> tuple[int, ...]:
    return tuple(int(value) for value in (values or ()))


@dataclass(frozen=True)
class WorkbenchConfig:
    workspace_bucket: str | None = None
    workspace_cdr: str | None = None
    requester_pays_project: str | None = None
    requester_pays_buckets: tuple[str, ...] = ()
    genomics_bucket: str = DEFAULT_GENOMICS_BUCKET
    storage_root: str = "./outputs/aou-workbench"
    wgs_vds_path: str = DEFAULT_WGS_VDS_PATH
    vat_path: str = DEFAULT_VAT_PATH
    exome_mt_path: str = DEFAULT_EXOME_MT_PATH
    clinvar_mt_path: str = DEFAULT_CLINVAR_MT_PATH
    acaf_mt_path: str = DEFAULT_ACAF_MT_PATH
    flagged_samples_path: str | None = None
    max_unrelated_path: str | None = None

    def __post_init__(self) -> None:
        requester_pays = tuple(
            bucket[5:] if bucket.startswith("gs://") else bucket
            for bucket in self.requester_pays_buckets
        )
        object.__setattr__(self, "requester_pays_buckets", requester_pays)
        if self.genomics_bucket and not self.requester_pays_buckets:
            object.__setattr__(
                self,
                "requester_pays_buckets",
                (self.genomics_bucket[5:],) if self.genomics_bucket.startswith("gs://") else (self.genomics_bucket,),
            )


@dataclass(frozen=True)
class CaseTierRule:
    name: str
    condition_concept_ids: tuple[int, ...] = ()
    condition_terms: tuple[str, ...] = ()
    measurement_concept_ids: tuple[int, ...] = ()
    measurement_terms: tuple[str, ...] = ()
    measurement_min: float | None = None
    require_condition: bool = True
    require_measurement: bool = False
    lookback_days: int = 30


@dataclass(frozen=True)
class ClinicalCofactorRule:
    name: str
    condition_concept_ids: tuple[int, ...] = ()
    condition_terms: tuple[str, ...] = ()


@dataclass(frozen=True)
class PhenotypeTables:
    cohort_table: str | None = None
    person_table: str = "person"
    observation_table: str = "observation_period"
    condition_table: str = "condition_occurrence"
    measurement_table: str = "measurement"
    concept_table: str = "concept"
    ancestry_table: str | None = None
    clinical_table: str | None = None


@dataclass(frozen=True)
class PhenotypeConfig:
    tables: PhenotypeTables
    person_id_column: str = "person_id"
    cohort_index_date_column: str | None = None
    birth_date_column: str | None = None
    age_column: str = "age"
    sex_column: str = "sex"
    observation_start_column: str = "observation_start_date"
    observation_end_column: str = "observation_end_date"
    condition_concept_column: str = "condition_concept_id"
    condition_date_column: str = "condition_start_date"
    measurement_concept_column: str = "measurement_concept_id"
    measurement_value_column: str = "value_as_number"
    measurement_date_column: str = "measurement_date"
    ancestry_person_id_column: str = "research_id"
    ancestry_label_column: str = "ancestry_pred"
    ancestry_pca_features_column: str = "pca_features"
    pc_columns: tuple[str, ...] = tuple(f"PC{i}" for i in range(1, 11))
    clinical_person_id_column: str = "person_id"
    clinical_cofactor_columns: tuple[str, ...] = ()
    clinical_cofactors: tuple[ClinicalCofactorRule, ...] = ()
    definite: CaseTierRule = field(
        default_factory=lambda: CaseTierRule(name="definite")
    )
    probable: CaseTierRule = field(
        default_factory=lambda: CaseTierRule(
            name="probable",
            require_condition=False,
            require_measurement=True,
        )
    )
    control_exclusion_concept_ids: tuple[int, ...] = ()
    min_observation_days: int = 180


@dataclass(frozen=True)
class CohortConfig:
    control_ratio: int = 4
    minimum_controls: int = 1
    age_tolerance_years: float = 5.0
    index_window_days: int = 365
    exact_match_columns: tuple[str, ...] = ("is_female", "ancestry_pred")
    primary_case_tier: str = "definite"
    sensitivity_case_tiers: tuple[str, ...] = ("probable",)
    covariate_columns: tuple[str, ...] = DEFAULT_COVARIATES
    require_complete_matches: bool = False
    allow_reuse_controls: bool = False


@dataclass(frozen=True)
class TargetVariant:
    label: str
    gene: str
    contig: str
    position: int
    ref: str
    alt: str
    rsid: str | None = None
    source: str | None = None
    evidence_tier: str | None = None
    exact_test_model: str = "carrier_vs_noncarrier"

    @property
    def variant_id(self) -> str:
        return f"{self.contig.replace('chr', '')}-{self.position}-{self.ref}-{self.alt}"


@dataclass(frozen=True)
class PanelConfig:
    genes_of_interest: tuple[str, ...] = DEFAULT_RHABDO_GENES
    burden_target_genes: tuple[str, ...] = DEFAULT_RHABDO_GENES
    a_priori_variants: tuple[TargetVariant, ...] = ()


@dataclass(frozen=True)
class Stage1Config:
    variant_table: str
    person_id_column: str = "person_id"
    variant_id_column: str = "variant_id"
    dosage_column: str = "dosage"
    gene_column: str = "gene"
    carrier_min_dosage: float = 1.0
    covariates: tuple[str, ...] = DEFAULT_COVARIATES


@dataclass(frozen=True)
class Stage2Config:
    variant_table: str
    person_id_column: str = "person_id"
    variant_id_column: str = "variant_id"
    gene_column: str = "gene"
    dosage_column: str = "dosage"
    clinvar_column: str = "clinvar_significance"
    consequence_column: str = "consequence"
    revel_column: str = "revel"
    af_column: str = "max_af"
    max_af: float = 0.01
    revel_min: float = 0.8
    plof_terms: tuple[str, ...] = DEFAULT_PLOF_TERMS
    clinvar_plp_terms: tuple[str, ...] = DEFAULT_CLINVAR_PLP_TERMS
    covariates: tuple[str, ...] = DEFAULT_COVARIATES


@dataclass(frozen=True)
class Stage3Config:
    variant_table: str
    mode: str = "targeted"
    run_exome_wide: bool = False
    person_id_column: str = "person_id"
    variant_id_column: str = "variant_id"
    gene_column: str = "gene"
    dosage_column: str = "dosage"
    clinvar_column: str = "clinvar_significance"
    consequence_column: str = "consequence"
    revel_column: str = "revel"
    af_column: str = "max_af"
    max_af: float = 0.01
    revel_min: float = 0.8
    masks: tuple[str, ...] = ("pLoF", "ClinVar_PLP", "pLoF_or_REVEL_0_8")
    plof_terms: tuple[str, ...] = DEFAULT_PLOF_TERMS
    clinvar_plp_terms: tuple[str, ...] = DEFAULT_CLINVAR_PLP_TERMS
    covariates: tuple[str, ...] = DEFAULT_COVARIATES


@dataclass(frozen=True)
class Stage4Config:
    genotype_table: str
    annotation_table: str | None = None
    person_id_column: str = "person_id"
    variant_id_column: str = "variant_id"
    dosage_column: str = "dosage"
    chromosome_column: str = "chromosome"
    position_column: str = "position"
    gene_column: str = "gene"
    af_column: str = "maf"
    covariates: tuple[str, ...] = DEFAULT_COVARIATES
    min_maf: float = 0.01
    min_minor_allele_count: float = 5.0
    pvalue_threshold: float = 5e-8
    qvalue_threshold: float = 0.1
    lead_hit_window_bp: int = 250_000


@dataclass(frozen=True)
class ReportingConfig:
    include_all_by_all_comparison: bool = False
    all_by_all_table: str | None = None
    top_n: int = 10


@dataclass(frozen=True)
class AnalysisConfig:
    analysis_name: str
    output_dir: str
    matched_outcome_column: str = "analysis_case"
    run_stage1: bool = True
    run_stage2: bool = True
    run_stage3: bool = True
    run_stage4: bool = True
    stage1: Stage1Config | None = None
    stage2: Stage2Config | None = None
    stage3: Stage3Config | None = None
    stage4: Stage4Config | None = None
    reporting: ReportingConfig = field(default_factory=ReportingConfig)

    @property
    def analysis_slug(self) -> str:
        return slugify(self.analysis_name)


@dataclass(frozen=True)
class ProjectConfig:
    workbench: WorkbenchConfig
    phenotype: PhenotypeConfig
    cohort: CohortConfig
    panel: PanelConfig
    analysis: AnalysisConfig

    @property
    def config_hash(self) -> str:
        return stable_hash(self.to_dict())

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _load_workbench(payload: Mapping[str, Any]) -> WorkbenchConfig:
    return WorkbenchConfig(
        workspace_bucket=payload.get("workspace_bucket"),
        workspace_cdr=payload.get("workspace_cdr"),
        requester_pays_project=payload.get("requester_pays_project"),
        requester_pays_buckets=_as_string_tuple(payload.get("requester_pays_buckets")),
        genomics_bucket=payload.get("genomics_bucket", DEFAULT_GENOMICS_BUCKET),
        storage_root=payload.get("storage_root", "./outputs/aou-workbench"),
        wgs_vds_path=payload.get("wgs_vds_path", DEFAULT_WGS_VDS_PATH),
        vat_path=payload.get("vat_path", DEFAULT_VAT_PATH),
        exome_mt_path=payload.get("exome_mt_path", DEFAULT_EXOME_MT_PATH),
        clinvar_mt_path=payload.get("clinvar_mt_path", DEFAULT_CLINVAR_MT_PATH),
        acaf_mt_path=payload.get("acaf_mt_path", DEFAULT_ACAF_MT_PATH),
        flagged_samples_path=payload.get("flagged_samples_path"),
        max_unrelated_path=payload.get("max_unrelated_path"),
    )


def _load_case_tier(name: str, payload: Mapping[str, Any] | None) -> CaseTierRule:
    payload = payload or {}
    return CaseTierRule(
        name=name,
        condition_concept_ids=_as_int_tuple(payload.get("condition_concept_ids")),
        condition_terms=_as_string_tuple(payload.get("condition_terms")),
        measurement_concept_ids=_as_int_tuple(payload.get("measurement_concept_ids")),
        measurement_terms=_as_string_tuple(payload.get("measurement_terms")),
        measurement_min=payload.get("measurement_min"),
        require_condition=bool(payload.get("require_condition", True)),
        require_measurement=bool(payload.get("require_measurement", False)),
        lookback_days=int(payload.get("lookback_days", 30)),
    )


def _load_clinical_cofactor(payload: Mapping[str, Any]) -> ClinicalCofactorRule:
    return ClinicalCofactorRule(
        name=payload["name"],
        condition_concept_ids=_as_int_tuple(payload.get("condition_concept_ids")),
        condition_terms=_as_string_tuple(payload.get("condition_terms")),
    )


def _load_phenotype(payload: Mapping[str, Any]) -> PhenotypeConfig:
    tables_payload = payload.get("tables", {})
    return PhenotypeConfig(
        tables=PhenotypeTables(
            cohort_table=tables_payload.get("cohort_table"),
            person_table=tables_payload.get("person_table", "person"),
            observation_table=tables_payload.get("observation_table", "observation_period"),
            condition_table=tables_payload.get("condition_table", "condition_occurrence"),
            measurement_table=tables_payload.get("measurement_table", "measurement"),
            concept_table=tables_payload.get("concept_table", "concept"),
            ancestry_table=tables_payload.get("ancestry_table", "person_ext"),
            clinical_table=tables_payload.get("clinical_table"),
        ),
        person_id_column=payload.get("person_id_column", "person_id"),
        cohort_index_date_column=payload.get("cohort_index_date_column"),
        birth_date_column=payload.get("birth_date_column"),
        age_column=payload.get("age_column", "age"),
        sex_column=payload.get("sex_column", "sex"),
        observation_start_column=payload.get("observation_start_column", "observation_start_date"),
        observation_end_column=payload.get("observation_end_column", "observation_end_date"),
        condition_concept_column=payload.get("condition_concept_column", "condition_concept_id"),
        condition_date_column=payload.get("condition_date_column", "condition_start_date"),
        measurement_concept_column=payload.get("measurement_concept_column", "measurement_concept_id"),
        measurement_value_column=payload.get("measurement_value_column", "value_as_number"),
        measurement_date_column=payload.get("measurement_date_column", "measurement_date"),
        ancestry_person_id_column=payload.get("ancestry_person_id_column", "research_id"),
        ancestry_label_column=payload.get("ancestry_label_column", "ancestry_pred"),
        ancestry_pca_features_column=payload.get("ancestry_pca_features_column", "pca_features"),
        pc_columns=_as_string_tuple(payload.get("pc_columns")) or tuple(f"PC{i}" for i in range(1, 11)),
        clinical_person_id_column=payload.get("clinical_person_id_column", "person_id"),
        clinical_cofactor_columns=_as_string_tuple(payload.get("clinical_cofactor_columns")),
        clinical_cofactors=tuple(
            _load_clinical_cofactor(item) for item in payload.get("clinical_cofactors", [])
        ),
        definite=_load_case_tier("definite", payload.get("definite")),
        probable=_load_case_tier("probable", payload.get("probable")),
        control_exclusion_concept_ids=_as_int_tuple(payload.get("control_exclusion_concept_ids")),
        min_observation_days=int(payload.get("min_observation_days", 180)),
    )


def _load_cohort(payload: Mapping[str, Any]) -> CohortConfig:
    return CohortConfig(
        control_ratio=int(payload.get("control_ratio", 4)),
        minimum_controls=int(payload.get("minimum_controls", 1)),
        age_tolerance_years=float(payload.get("age_tolerance_years", 5.0)),
        index_window_days=int(payload.get("index_window_days", 365)),
        exact_match_columns=_as_string_tuple(payload.get("exact_match_columns")) or ("is_female", "ancestry_pred"),
        primary_case_tier=payload.get("primary_case_tier", "definite"),
        sensitivity_case_tiers=_as_string_tuple(payload.get("sensitivity_case_tiers")) or ("probable",),
        covariate_columns=_as_string_tuple(payload.get("covariate_columns")) or DEFAULT_COVARIATES,
        require_complete_matches=bool(payload.get("require_complete_matches", False)),
        allow_reuse_controls=bool(payload.get("allow_reuse_controls", False)),
    )


def _load_target_variant(payload: Mapping[str, Any]) -> TargetVariant:
    return TargetVariant(
        label=payload["label"],
        gene=payload["gene"],
        contig=payload["contig"],
        position=int(payload["position"]),
        ref=payload["ref"],
        alt=payload["alt"],
        rsid=payload.get("rsid"),
        source=payload.get("source"),
        evidence_tier=payload.get("evidence_tier"),
        exact_test_model=payload.get("exact_test_model", "carrier_vs_noncarrier"),
    )


def _load_panel(payload: Mapping[str, Any]) -> PanelConfig:
    return PanelConfig(
        genes_of_interest=_as_string_tuple(payload.get("genes_of_interest")) or DEFAULT_RHABDO_GENES,
        burden_target_genes=_as_string_tuple(payload.get("burden_target_genes")) or DEFAULT_RHABDO_GENES,
        a_priori_variants=tuple(
            _load_target_variant(item) for item in payload.get("a_priori_variants", [])
        ),
    )


def _load_stage1(payload: Mapping[str, Any]) -> Stage1Config:
    return Stage1Config(
        variant_table=payload["variant_table"],
        person_id_column=payload.get("person_id_column", "person_id"),
        variant_id_column=payload.get("variant_id_column", "variant_id"),
        dosage_column=payload.get("dosage_column", "dosage"),
        gene_column=payload.get("gene_column", "gene"),
        carrier_min_dosage=float(payload.get("carrier_min_dosage", 1.0)),
        covariates=_as_string_tuple(payload.get("covariates")) or DEFAULT_COVARIATES,
    )


def _load_stage2(payload: Mapping[str, Any]) -> Stage2Config:
    return Stage2Config(
        variant_table=payload["variant_table"],
        person_id_column=payload.get("person_id_column", "person_id"),
        variant_id_column=payload.get("variant_id_column", "variant_id"),
        gene_column=payload.get("gene_column", "gene"),
        dosage_column=payload.get("dosage_column", "dosage"),
        clinvar_column=payload.get("clinvar_column", "clinvar_significance"),
        consequence_column=payload.get("consequence_column", "consequence"),
        revel_column=payload.get("revel_column", "revel"),
        af_column=payload.get("af_column", "max_af"),
        max_af=float(payload.get("max_af", 0.01)),
        revel_min=float(payload.get("revel_min", 0.8)),
        plof_terms=_as_string_tuple(payload.get("plof_terms")) or DEFAULT_PLOF_TERMS,
        clinvar_plp_terms=_as_string_tuple(payload.get("clinvar_plp_terms")) or DEFAULT_CLINVAR_PLP_TERMS,
        covariates=_as_string_tuple(payload.get("covariates")) or DEFAULT_COVARIATES,
    )


def _load_stage3(payload: Mapping[str, Any]) -> Stage3Config:
    return Stage3Config(
        variant_table=payload["variant_table"],
        mode=payload.get("mode", "targeted"),
        run_exome_wide=bool(payload.get("run_exome_wide", False)),
        person_id_column=payload.get("person_id_column", "person_id"),
        variant_id_column=payload.get("variant_id_column", "variant_id"),
        gene_column=payload.get("gene_column", "gene"),
        dosage_column=payload.get("dosage_column", "dosage"),
        clinvar_column=payload.get("clinvar_column", "clinvar_significance"),
        consequence_column=payload.get("consequence_column", "consequence"),
        revel_column=payload.get("revel_column", "revel"),
        af_column=payload.get("af_column", "max_af"),
        max_af=float(payload.get("max_af", 0.01)),
        revel_min=float(payload.get("revel_min", 0.8)),
        masks=_as_string_tuple(payload.get("masks")) or ("pLoF", "ClinVar_PLP", "pLoF_or_REVEL_0_8"),
        plof_terms=_as_string_tuple(payload.get("plof_terms")) or DEFAULT_PLOF_TERMS,
        clinvar_plp_terms=_as_string_tuple(payload.get("clinvar_plp_terms")) or DEFAULT_CLINVAR_PLP_TERMS,
        covariates=_as_string_tuple(payload.get("covariates")) or DEFAULT_COVARIATES,
    )


def _load_stage4(payload: Mapping[str, Any]) -> Stage4Config:
    return Stage4Config(
        genotype_table=payload["genotype_table"],
        annotation_table=payload.get("annotation_table"),
        person_id_column=payload.get("person_id_column", "person_id"),
        variant_id_column=payload.get("variant_id_column", "variant_id"),
        dosage_column=payload.get("dosage_column", "dosage"),
        chromosome_column=payload.get("chromosome_column", "chromosome"),
        position_column=payload.get("position_column", "position"),
        gene_column=payload.get("gene_column", "gene"),
        af_column=payload.get("af_column", "maf"),
        covariates=_as_string_tuple(payload.get("covariates")) or DEFAULT_COVARIATES,
        min_maf=float(payload.get("min_maf", 0.01)),
        min_minor_allele_count=float(payload.get("min_minor_allele_count", 5.0)),
        pvalue_threshold=float(payload.get("pvalue_threshold", 5e-8)),
        qvalue_threshold=float(payload.get("qvalue_threshold", 0.1)),
        lead_hit_window_bp=int(payload.get("lead_hit_window_bp", 250_000)),
    )


def _load_reporting(payload: Mapping[str, Any] | None) -> ReportingConfig:
    payload = payload or {}
    return ReportingConfig(
        include_all_by_all_comparison=bool(payload.get("include_all_by_all_comparison", False)),
        all_by_all_table=payload.get("all_by_all_table"),
        top_n=int(payload.get("top_n", 10)),
    )


def _load_analysis(payload: Mapping[str, Any]) -> AnalysisConfig:
    return AnalysisConfig(
        analysis_name=payload["analysis_name"],
        output_dir=payload["output_dir"],
        matched_outcome_column=payload.get("matched_outcome_column", "analysis_case"),
        run_stage1=bool(payload.get("run_stage1", True)),
        run_stage2=bool(payload.get("run_stage2", True)),
        run_stage3=bool(payload.get("run_stage3", True)),
        run_stage4=bool(payload.get("run_stage4", True)),
        stage1=_load_stage1(payload["stage1"]) if payload.get("stage1") else None,
        stage2=_load_stage2(payload["stage2"]) if payload.get("stage2") else None,
        stage3=_load_stage3(payload["stage3"]) if payload.get("stage3") else None,
        stage4=_load_stage4(payload["stage4"]) if payload.get("stage4") else None,
        reporting=_load_reporting(payload.get("reporting")),
    )


def load_project_config(
    *,
    workbench_path: str,
    phenotype_path: str,
    cohort_path: str,
    panel_path: str,
    analysis_path: str,
) -> ProjectConfig:
    return ProjectConfig(
        workbench=_load_workbench(load_yaml(workbench_path)),
        phenotype=_load_phenotype(load_yaml(phenotype_path)),
        cohort=_load_cohort(load_yaml(cohort_path)),
        panel=_load_panel(load_yaml(panel_path)),
        analysis=_load_analysis(load_yaml(analysis_path)),
    )


def build_default_target_panel() -> tuple[TargetVariant, ...]:
    return (
        TargetVariant(
            label="HBB sickle trait",
            gene="HBB",
            rsid="rs334",
            contig="chr11",
            position=5227002,
            ref="T",
            alt="A",
            source="Deuster2013",
            evidence_tier="candidate_common",
        ),
        TargetVariant(
            label="RYR1 Arg401Cys",
            gene="RYR1",
            rsid="rs193922764",
            contig="chr19",
            position=38451842,
            ref="C",
            alt="T",
            source="a_priori_RYR1",
            evidence_tier="core_pathogenic",
        ),
        TargetVariant(
            label="RYR1 Gly2434Arg",
            gene="RYR1",
            rsid="rs121918593",
            contig="chr19",
            position=38499993,
            ref="G",
            alt="A",
            source="a_priori_RYR1",
            evidence_tier="core_pathogenic",
        ),
        TargetVariant(
            label="CPT2 Ser113Leu",
            gene="CPT2",
            rsid="rs74315294",
            contig="chr1",
            position=53202427,
            ref="C",
            alt="T",
            source="CPT2_myopathic_form",
            evidence_tier="core_pathogenic",
        ),
        TargetVariant(
            label="CACNA1S Arg1086His",
            gene="CACNA1S",
            rsid="rs1800559",
            contig="chr1",
            position=201060815,
            ref="C",
            alt="T",
            source="MH_ER_related",
            evidence_tier="core_pathogenic",
        ),
        TargetVariant(
            label="PYGM Arg50Ter",
            gene="PYGM",
            rsid="rs116987552",
            contig="chr11",
            position=64526623,
            ref="C",
            alt="T",
            source="McArdle_candidate",
            evidence_tier="core_pathogenic",
        ),
    )


__all__ = [
    "AnalysisConfig",
    "CaseTierRule",
    "CohortConfig",
    "DEFAULT_CLINVAR_PLP_TERMS",
    "DEFAULT_COVARIATES",
    "DEFAULT_PLOF_TERMS",
    "DEFAULT_WORKSPACE_CDR",
    "PanelConfig",
    "PhenotypeConfig",
    "PhenotypeTables",
    "ProjectConfig",
    "ReportingConfig",
    "Stage1Config",
    "Stage2Config",
    "Stage3Config",
    "Stage4Config",
    "TargetVariant",
    "WorkbenchConfig",
    "build_default_target_panel",
    "load_project_config",
]
