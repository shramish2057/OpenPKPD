module OpenPKPDCore

# ------------------------------------------------------------------
# External dependencies
# ------------------------------------------------------------------
using SciMLBase
using DifferentialEquations

const OPENPKPD_VERSION = "0.1.0"
export OPENPKPD_VERSION

# ------------------------------------------------------------------
# Core specs and shared types
# ------------------------------------------------------------------
include("specs/time_covariates.jl")
include("specs/specs.jl")
include("specs/error_models.jl")
include("specs/sensitivity.jl")

# ------------------------------------------------------------------
# PK model definitions
# ------------------------------------------------------------------
include("models/onecomp_iv_bolus.jl")
include("models/onecomp_oral_first_order.jl")
include("models/twocomp_iv_bolus.jl")
include("models/twocomp_oral.jl")
include("models/threecomp_iv_bolus.jl")
include("models/transit_absorption.jl")
include("models/michaelis_menten.jl")
include("models/pk_interface.jl")

# ------------------------------------------------------------------
# PD model definitions
# ------------------------------------------------------------------
include("pd/direct_emax.jl")
include("pd/sigmoid_emax.jl")
include("pd/indirect_response_turnover.jl")
include("pd/biophase_equilibration.jl")

# ------------------------------------------------------------------
# Numerical semantics 
# These define versioned scientific meaning
# ------------------------------------------------------------------
include("engine/semantics.jl")
include("engine/solver_semantics.jl")
include("engine/semantics_fingerprint.jl")

# ------------------------------------------------------------------
# Perturbation + sensitivity core
# ------------------------------------------------------------------
include("engine/perturb.jl")
include("engine/sensitivity_metrics.jl")
include("engine/sensitivity.jl")
include("engine/sensitivity_population.jl")

# ------------------------------------------------------------------
# Core simulation engine
# ------------------------------------------------------------------
include("engine/events.jl")
include("engine/callbacks.jl")
include("engine/infusion.jl")
include("engine/residual_error.jl")
include("engine/solve.jl")

# ------------------------------------------------------------------
# PK–PD execution layers
# ------------------------------------------------------------------
include("engine/pkpd.jl")
include("engine/pkpd_coupled.jl")

# ------------------------------------------------------------------
# Population engine
# Defines PopulationSpec, PopulationResult, IIV, etc.
# ------------------------------------------------------------------
include("engine/iov.jl")
include("engine/segment_sim.jl")
include("engine/segment_sim_pkpd.jl")
include("engine/time_covariates.jl")
include("engine/covariates.jl")
include("engine/population.jl")

# ------------------------------------------------------------------
# Serialization
# ------------------------------------------------------------------
include("serialization/schema.jl")
include("serialization/serialize.jl")
include("serialization/deserialize.jl")
include("serialization/serialize_population.jl")
include("serialization/deserialize_population.jl")
include("serialization/serialize_sensitivity.jl")
include("serialization/deserialize_sensitivity.jl")
include("serialization/serialize_error.jl")
include("serialization/deserialize_error.jl")

include("analysis/exposure.jl")
include("analysis/response_metrics.jl")

# ------------------------------------------------------------------
# Model Import (NONMEM, Monolix)
# ------------------------------------------------------------------
include("import/import.jl")

# ------------------------------------------------------------------
# Data Module (CDISC/SDTM support)
# ------------------------------------------------------------------
include("data/data.jl")

# ------------------------------------------------------------------
# Visual Predictive Check (VPC)
# Requires data module for ObservedData
# ------------------------------------------------------------------
include("analysis/vpc_types.jl")
include("analysis/binning.jl")
include("analysis/vpc.jl")

# ------------------------------------------------------------------
# Parameter Estimation (NLME)
# FOCE-I, SAEM, Laplacian methods
# ------------------------------------------------------------------
include("estimation/estimation.jl")

# ------------------------------------------------------------------
# NCA (Non-Compartmental Analysis)
# FDA/EMA compliant NCA metrics
# ------------------------------------------------------------------
include("nca/nca.jl")

# ------------------------------------------------------------------
# Clinical Trial Simulation
# Parallel, crossover, dose-escalation, adaptive designs
# Virtual population, power analysis, bioequivalence
# ------------------------------------------------------------------
include("trial/trial.jl")

end
