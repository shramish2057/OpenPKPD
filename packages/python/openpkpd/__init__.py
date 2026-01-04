"""
OpenPKPD - Professional PK/PD Simulation Platform

A Python interface to the OpenPKPD Julia core for pharmacokinetic and
pharmacodynamic simulations.

Features:
- Single PK simulations:
  - One-compartment: IV bolus, IV infusion, oral first-order
  - Two-compartment: IV bolus, oral first-order
  - Three-compartment: IV bolus
  - Advanced absorption: Transit compartment model
  - Nonlinear elimination: Michaelis-Menten kinetics
- Coupled PK-PD simulations:
  - Direct Emax model
  - Sigmoid Emax model (Hill equation)
  - Biophase equilibration (effect compartment)
  - Indirect response turnover model
- Population simulations with IIV, IOV, and covariates
- Parameter estimation (FOCE-I, SAEM, Laplacian)
- Visual Predictive Check (VPC, pcVPC)
- Non-Compartmental Analysis (NCA)
- CDISC/SDTM data import
- NONMEM/Monolix model import
- Residual error models
- Sensitivity analysis
- Artifact replay for reproducibility
- PK/PD metrics (Cmax, AUC, Tmax, Emin, time below threshold)

Quick Start:
    >>> import openpkpd
    >>> openpkpd.init_julia()  # Initialize once per session

    >>> # Run a simple PK simulation
    >>> result = openpkpd.simulate_pk_iv_bolus(
    ...     cl=1.0, v=10.0,
    ...     doses=[{"time": 0.0, "amount": 100.0}],
    ...     t0=0.0, t1=24.0,
    ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
    ... )

    >>> # Run IV infusion
    >>> result = openpkpd.simulate_pk_iv_bolus(
    ...     cl=1.0, v=10.0,
    ...     doses=[{"time": 0.0, "amount": 100.0, "duration": 1.0}],  # 1-hour infusion
    ...     t0=0.0, t1=24.0,
    ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
    ... )

    >>> # Run parameter estimation
    >>> from openpkpd.estimation import estimate, EstimationConfig
    >>> config = EstimationConfig(method="foce", theta_init=[1.0, 10.0])
    >>> result = estimate(observed_data, "OneCompIVBolus", config, grid)

    >>> # Run VPC
    >>> from openpkpd.vpc import compute_vpc, VPCConfig
    >>> vpc_result = compute_vpc(observed, pop_spec, grid)

For more information, see the docstrings for individual functions.
"""

# Core initialization and utilities
from ._core import (
    init_julia,
    version,
    SensitivityMetrics,
    SensitivityResult,
    PopulationSensitivityResult,
)

# Artifact operations
from .artifacts import (
    replay_artifact,
    write_single_artifact,
)

# PK/PD metrics
from .metrics import (
    cmax,
    auc_trapezoid,
    emin,
    time_below,
    auc_above_baseline,
    tmax,
    half_life,
)

# Simulations - One-compartment PK
from .simulations.pk_onecomp import (
    simulate_pk_iv_bolus,
    simulate_pk_oral_first_order,
)

# Simulations - Two-compartment PK
from .simulations.pk_twocomp import (
    simulate_pk_twocomp_iv_bolus,
    simulate_pk_twocomp_oral,
)

# Simulations - Three-compartment PK
from .simulations.pk_threecomp import (
    simulate_pk_threecomp_iv_bolus,
)

# Simulations - Advanced PK
from .simulations.pk_advanced import (
    simulate_pk_transit_absorption,
    simulate_pk_michaelis_menten,
)

# Simulations - Coupled PKPD
from .simulations.pkpd import (
    simulate_pkpd_direct_emax,
    simulate_pkpd_indirect_response,
    simulate_pkpd_sigmoid_emax,
    simulate_pkpd_biophase_equilibration,
)

# Simulations - Population
from .simulations.population import (
    simulate_population_iv_bolus,
    simulate_population_oral,
)

# Simulations - Sensitivity
from .simulations.sensitivity import (
    run_sensitivity,
)

# NCA - Non-Compartmental Analysis
from .nca import (
    run_nca,
    run_population_nca,
    summarize_population_nca,
    NCAConfig,
    NCAResult,
    nca_cmax,
    nca_tmax,
    auc_0_t,
    auc_0_inf,
    estimate_lambda_z,
    nca_half_life,
    bioequivalence_90ci,
    tost_analysis,
    be_conclusion,
    geometric_mean_ratio,
)

# Parameter Estimation (NLME)
from .estimation import (
    estimate,
    EstimationConfig,
    EstimationResult,
    IndividualEstimate,
)

# Visual Predictive Check
from .vpc import (
    compute_vpc,
    compute_pcvpc,
    VPCConfig,
    VPCResult,
    VPCBin,
    VPCPercentile,
)

# CDISC Data Import
from .data import (
    read_cdisc_pc,
    read_cdisc_ex,
    read_cdisc_dm,
    read_cdisc_dataset,
    cdisc_to_observed_data,
    validate_cdisc_dataset,
    SubjectData,
    ObservedData,
)

# Model Import (NONMEM/Monolix)
from .import_ import (
    import_nonmem,
    import_monolix,
    import_model,
    parse_nonmem_control,
    parse_monolix_project,
    ImportedModel,
    NONMEMControlFile,
    MonolixProject,
)


__all__ = [
    # Core
    "init_julia",
    "version",

    # Data classes
    "SensitivityMetrics",
    "SensitivityResult",
    "PopulationSensitivityResult",

    # Artifacts
    "replay_artifact",
    "write_single_artifact",

    # Single PK - One-compartment
    "simulate_pk_iv_bolus",
    "simulate_pk_oral_first_order",

    # Single PK - Two-compartment
    "simulate_pk_twocomp_iv_bolus",
    "simulate_pk_twocomp_oral",

    # Single PK - Three-compartment
    "simulate_pk_threecomp_iv_bolus",

    # Single PK - Advanced models
    "simulate_pk_transit_absorption",
    "simulate_pk_michaelis_menten",

    # Coupled PKPD - Direct PD models
    "simulate_pkpd_direct_emax",
    "simulate_pkpd_sigmoid_emax",
    "simulate_pkpd_biophase_equilibration",

    # Coupled PKPD - ODE-based PD models
    "simulate_pkpd_indirect_response",

    # Population
    "simulate_population_iv_bolus",
    "simulate_population_oral",

    # Sensitivity
    "run_sensitivity",

    # Metrics
    "cmax",
    "auc_trapezoid",
    "emin",
    "time_below",
    "auc_above_baseline",
    "tmax",
    "half_life",

    # NCA
    "run_nca",
    "run_population_nca",
    "summarize_population_nca",
    "NCAConfig",
    "NCAResult",
    "nca_cmax",
    "nca_tmax",
    "auc_0_t",
    "auc_0_inf",
    "estimate_lambda_z",
    "nca_half_life",
    "bioequivalence_90ci",
    "tost_analysis",
    "be_conclusion",
    "geometric_mean_ratio",

    # Estimation (NLME)
    "estimate",
    "EstimationConfig",
    "EstimationResult",
    "IndividualEstimate",

    # VPC
    "compute_vpc",
    "compute_pcvpc",
    "VPCConfig",
    "VPCResult",
    "VPCBin",
    "VPCPercentile",

    # CDISC Data
    "read_cdisc_pc",
    "read_cdisc_ex",
    "read_cdisc_dm",
    "read_cdisc_dataset",
    "cdisc_to_observed_data",
    "validate_cdisc_dataset",
    "SubjectData",
    "ObservedData",

    # Model Import
    "import_nonmem",
    "import_monolix",
    "import_model",
    "parse_nonmem_control",
    "parse_monolix_project",
    "ImportedModel",
    "NONMEMControlFile",
    "MonolixProject",
]
