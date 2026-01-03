# Specifications are pure data. No solver logic and no hidden defaults.
# -------------------------
# Shared validation helpers
# -------------------------

function _require_positive(name::String, x::Float64)
    if !(x > 0.0)
        error("Expected positive value for $(name), got $(x)")
    end
    return nothing
end

# -------------------------
# PK specifications
# -------------------------
export ModelKind,
    OneCompIVBolus, OneCompOralFirstOrder, OneCompIVBolusParams, OneCompOralFirstOrderParams
export DoseEvent, ModelSpec
export SolverSpec, SimGrid, SimResult

abstract type ModelKind end

struct OneCompIVBolus <: ModelKind end
struct OneCompOralFirstOrder <: ModelKind end

struct DoseEvent
    time::Float64
    amount::Float64
end

struct OneCompIVBolusParams
    CL::Float64
    V::Float64
end

struct OneCompOralFirstOrderParams
    Ka::Float64
    CL::Float64
    V::Float64
end

struct ModelSpec{K<:ModelKind,P}
    kind::K
    name::String
    params::P
    doses::Vector{DoseEvent}
end

struct SolverSpec
    alg::Symbol
    reltol::Float64
    abstol::Float64
    maxiters::Int
end

struct SimGrid
    t0::Float64
    t1::Float64
    saveat::Vector{Float64}
end

struct SimResult
    t::Vector{Float64}
    states::Dict{Symbol,Vector{Float64}}
    observations::Dict{Symbol,Vector{Float64}}
    metadata::Dict{String,Any}
end

# -------------------------
# PD specifications
# -------------------------

export PDModelKind, DirectEmax, DirectEmaxParams, PDSpec

abstract type PDModelKind end

"""
Direct Emax PD model.

Effect(C) = E0 + (Emax * C) / (EC50 + C)
"""
struct DirectEmax <: PDModelKind end

struct DirectEmaxParams
    E0::Float64
    Emax::Float64
    EC50::Float64
end

export IndirectResponseTurnover, IndirectResponseTurnoverParams

"""
Indirect response turnover PD model with inhibition of Kout.

States:
- R(t): response

Effect:
I(C) = (Imax * C) / (IC50 + C)

Dynamics:
dR/dt = Kin - Kout * (1 - I(C)) * R
"""
struct IndirectResponseTurnover <: PDModelKind end

struct IndirectResponseTurnoverParams
    Kin::Float64
    Kout::Float64
    R0::Float64
    Imax::Float64
    IC50::Float64
end

"""
PD specification container.

input_observation:
- which observation key from the upstream system is used as input, usually :conc

output_observation:
- name of the produced PD observable, default :effect is typical
"""
struct PDSpec{K<:PDModelKind,P}
    kind::K
    name::String
    params::P
    input_observation::Symbol
    output_observation::Symbol
end

# -------------------------
# Population specifications
# -------------------------

export RandomEffectKind, LogNormalIIV, IIVSpec, PopulationSpec, IndividualCovariates

abstract type RandomEffectKind end

"""
Log-normal inter-individual variability (IIV).

Parameter transform:
theta_i = theta_pop * exp(eta_i)

eta_i ~ Normal(0, omega^2)
"""
struct LogNormalIIV <: RandomEffectKind end

"""
IIV specification for a set of parameters.

omegas:
- Dict mapping parameter symbol to omega (standard deviation of eta)

seed:
- deterministic seed for RNG

n:
- number of individuals
"""
struct IIVSpec{K<:RandomEffectKind}
    kind::K
    omegas::Dict{Symbol,Float64}
    seed::UInt64
    n::Int
end

export IOVSpec, OccasionDefinition

"""
OccasionDefinition defines how dosing occasions are determined.

Supported v1 mode:
- :dose_times -> each unique dose time strictly greater than t0 starts a new occasion
- t0 is occasion 1
"""
struct OccasionDefinition
    mode::Symbol
end

"""
IOV specification.

pis:
- Dict mapping parameter symbol to pi (std dev of kappa)

seed:
- deterministic seed for IOV RNG stream (separate from IIV)

occasion_def:
- how occasions are determined
"""
struct IOVSpec{K<:RandomEffectKind}
    kind::K
    pis::Dict{Symbol,Float64}
    seed::UInt64
    occasion_def::OccasionDefinition
end

export CovariateEffectKind, LinearCovariate, PowerCovariate, ExpCovariate
export CovariateEffect, CovariateModel

abstract type CovariateEffectKind end

"""
Linear covariate model:
theta_i = theta_pop * (1 + beta * (cov - ref))
"""
struct LinearCovariate <: CovariateEffectKind end

"""
Power covariate model:
theta_i = theta_pop * (cov / ref) ^ beta
"""
struct PowerCovariate <: CovariateEffectKind end

"""
Exponential covariate model:
theta_i = theta_pop * exp(beta * (cov - ref))
"""
struct ExpCovariate <: CovariateEffectKind end

struct CovariateEffect{K<:CovariateEffectKind}
    kind::K
    param::Symbol
    covariate::Symbol
    beta::Float64
    ref::Float64
end

"""
A set of covariate effects applied to parameters.
"""
struct CovariateModel
    name::String
    effects::Vector{CovariateEffect}
end

"""
Optional covariates per individual.

values:
- static covariates as Dict{Symbol, Float64}

time_varying:
- optional TimeVaryingCovariates for time-dependent covariates
"""
struct IndividualCovariates
    values::Dict{Symbol,Float64}
    time_varying::Union{Nothing,TimeVaryingCovariates}
end

"""
Population simulation specification.

base_model_spec:
- the typical value model (population thetas) with doses, grid, solver handled outside
iiv:
- optional IIV spec (can be nothing)
covariates:
- optional vector aligned with n individuals (can be empty)
"""
struct PopulationSpec{MS}
    base_model_spec::MS
    iiv::Union{Nothing,IIVSpec}
    iov::Union{Nothing,IOVSpec}
    covariate_model::Union{Nothing,CovariateModel}
    covariates::Vector{IndividualCovariates}
end
