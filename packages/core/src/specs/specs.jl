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
export TwoCompIVBolus, TwoCompIVBolusParams, TwoCompOral, TwoCompOralParams
export ThreeCompIVBolus, ThreeCompIVBolusParams
export TransitAbsorption, TransitAbsorptionParams
export MichaelisMentenElimination, MichaelisMentenEliminationParams
export DoseEvent, ModelSpec
export SolverSpec, SimGrid, SimResult

abstract type ModelKind end

struct OneCompIVBolus <: ModelKind end
struct OneCompOralFirstOrder <: ModelKind end

"""
Dose event specification supporting both bolus and infusion administration.

For IV bolus: duration = 0.0 (default), dose is instantaneously added to compartment
For IV infusion: duration > 0.0, dose is delivered at constant rate over duration

Fields:
- time: Start time of dose administration
- amount: Total drug amount to be administered
- duration: Infusion duration in time units (0.0 = instantaneous bolus)

The infusion rate is computed as: rate = amount / duration

Examples:
- DoseEvent(0.0, 100.0) - 100 mg bolus at t=0
- DoseEvent(0.0, 100.0, 0.0) - same as above, explicit bolus
- DoseEvent(0.0, 100.0, 1.0) - 100 mg infused over 1 hour (rate = 100 mg/h)
- DoseEvent(0.0, 100.0, 0.5) - 100 mg infused over 30 min (rate = 200 mg/h)
"""
struct DoseEvent
    time::Float64
    amount::Float64
    duration::Float64

    function DoseEvent(time::Float64, amount::Float64, duration::Float64=0.0)
        duration >= 0.0 || error("DoseEvent duration must be >= 0, got $(duration)")
        amount >= 0.0 || error("DoseEvent amount must be >= 0, got $(amount)")
        new(time, amount, duration)
    end
end

# Note: Two-argument constructor DoseEvent(time, amount) works via the default
# duration=0.0 in the inner constructor above.

"""
Check if dose event is a bolus (instantaneous) administration.
"""
is_bolus(dose::DoseEvent) = dose.duration == 0.0

"""
Check if dose event is an infusion (zero-order input).
"""
is_infusion(dose::DoseEvent) = dose.duration > 0.0

"""
Compute the infusion rate for a dose event.
Returns the rate in amount/time units. For bolus doses, returns Inf.
"""
function infusion_rate(dose::DoseEvent)
    if is_bolus(dose)
        return Inf
    else
        return dose.amount / dose.duration
    end
end

"""
Get the end time of a dose event (when infusion stops).
For bolus doses, this equals the start time.
"""
dose_end_time(dose::DoseEvent) = dose.time + dose.duration

export is_bolus, is_infusion, infusion_rate, dose_end_time

struct OneCompIVBolusParams
    CL::Float64
    V::Float64
end

struct OneCompOralFirstOrderParams
    Ka::Float64
    CL::Float64
    V::Float64
end

# -------------------------
# Two-compartment PK models
# -------------------------

"""
Two-compartment IV bolus PK model.

States:
- A_central: Amount in central compartment
- A_peripheral: Amount in peripheral compartment

Parameters:
- CL: Clearance from central compartment (volume/time)
- V1: Volume of central compartment
- Q: Inter-compartmental clearance (volume/time)
- V2: Volume of peripheral compartment

Micro-constants:
- k10 = CL/V1 (elimination rate constant)
- k12 = Q/V1 (central to peripheral rate constant)
- k21 = Q/V2 (peripheral to central rate constant)

Dynamics:
dA_central/dt = -k10*A_central - k12*A_central + k21*A_peripheral
dA_peripheral/dt = k12*A_central - k21*A_peripheral
"""
struct TwoCompIVBolus <: ModelKind end

struct TwoCompIVBolusParams
    CL::Float64   # Clearance
    V1::Float64   # Central volume
    Q::Float64    # Inter-compartmental clearance
    V2::Float64   # Peripheral volume
end

"""
Two-compartment oral first-order absorption PK model.

States:
- A_gut: Amount in gut compartment
- A_central: Amount in central compartment
- A_peripheral: Amount in peripheral compartment

Parameters:
- Ka: Absorption rate constant (1/time)
- CL: Clearance from central compartment (volume/time)
- V1: Volume of central compartment
- Q: Inter-compartmental clearance (volume/time)
- V2: Volume of peripheral compartment

Dynamics:
dA_gut/dt = -Ka*A_gut
dA_central/dt = Ka*A_gut - (CL/V1)*A_central - (Q/V1)*A_central + (Q/V2)*A_peripheral
dA_peripheral/dt = (Q/V1)*A_central - (Q/V2)*A_peripheral
"""
struct TwoCompOral <: ModelKind end

struct TwoCompOralParams
    Ka::Float64   # Absorption rate constant
    CL::Float64   # Clearance
    V1::Float64   # Central volume
    Q::Float64    # Inter-compartmental clearance
    V2::Float64   # Peripheral volume
end

# -------------------------
# Three-compartment PK model
# -------------------------

"""
Three-compartment IV bolus PK model (mammillary).

States:
- A_central: Amount in central compartment
- A_periph1: Amount in first peripheral (shallow) compartment
- A_periph2: Amount in second peripheral (deep) compartment

Parameters:
- CL: Clearance from central compartment (volume/time)
- V1: Volume of central compartment
- Q2: Inter-compartmental clearance to shallow peripheral (volume/time)
- V2: Volume of shallow peripheral compartment
- Q3: Inter-compartmental clearance to deep peripheral (volume/time)
- V3: Volume of deep peripheral compartment

Dynamics:
dA_central/dt = -(CL/V1)*A_central - (Q2/V1)*A_central + (Q2/V2)*A_periph1
                - (Q3/V1)*A_central + (Q3/V3)*A_periph2
dA_periph1/dt = (Q2/V1)*A_central - (Q2/V2)*A_periph1
dA_periph2/dt = (Q3/V1)*A_central - (Q3/V3)*A_periph2
"""
struct ThreeCompIVBolus <: ModelKind end

struct ThreeCompIVBolusParams
    CL::Float64   # Clearance
    V1::Float64   # Central volume
    Q2::Float64   # Inter-compartmental clearance (shallow)
    V2::Float64   # Shallow peripheral volume
    Q3::Float64   # Inter-compartmental clearance (deep)
    V3::Float64   # Deep peripheral volume
end

# -------------------------
# Transit absorption model
# -------------------------

"""
Transit compartment absorption model.

This model implements a chain of transit compartments before the absorption
compartment, providing a delayed and more physiological absorption profile.
Based on Savic et al. (2007) transit compartment model.

States:
- Transit[1:N]: Amount in each transit compartment
- A_central: Amount in central compartment

Parameters:
- N: Number of transit compartments (integer >= 1)
- Ktr: Transit rate constant (1/time) - same for all transit compartments
- Ka: Absorption rate constant from last transit to central (1/time)
- CL: Clearance (volume/time)
- V: Volume of distribution

Dynamics:
For i = 1: dTransit[1]/dt = -Ktr * Transit[1]  (receives dose)
For i > 1: dTransit[i]/dt = Ktr * Transit[i-1] - Ktr * Transit[i]
dA_central/dt = Ka * Transit[N] - (CL/V) * A_central

Note: Mean transit time (MTT) ≈ (N+1) / Ktr
"""
struct TransitAbsorption <: ModelKind end

struct TransitAbsorptionParams
    N::Int        # Number of transit compartments
    Ktr::Float64  # Transit rate constant
    Ka::Float64   # Absorption rate constant
    CL::Float64   # Clearance
    V::Float64    # Volume of distribution
end

# -------------------------
# Michaelis-Menten elimination model
# -------------------------

"""
One-compartment PK model with Michaelis-Menten (saturable) elimination.

This model describes nonlinear pharmacokinetics where the elimination
pathway becomes saturated at higher concentrations.

States:
- A_central: Amount in central compartment

Parameters:
- Vmax: Maximum elimination rate (mass/time)
- Km: Michaelis constant - concentration at half Vmax (mass/volume)
- V: Volume of distribution

Dynamics:
C = A_central / V
dA_central/dt = -Vmax * C / (Km + C)
             = -Vmax * A_central / (Km * V + A_central)

At low concentrations (C << Km): Approximates first-order with CL ≈ Vmax/Km
At high concentrations (C >> Km): Approximates zero-order with rate ≈ Vmax
"""
struct MichaelisMentenElimination <: ModelKind end

struct MichaelisMentenEliminationParams
    Vmax::Float64  # Maximum elimination rate
    Km::Float64    # Michaelis constant
    V::Float64     # Volume of distribution
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
export SigmoidEmax, SigmoidEmaxParams
export BiophaseEquilibration, BiophaseEquilibrationParams

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

"""
Sigmoid Emax PD model (Hill equation).

This model extends the direct Emax model with a Hill coefficient (gamma)
that controls the steepness of the concentration-effect relationship.

Effect(C) = E0 + (Emax * C^gamma) / (EC50^gamma + C^gamma)

Parameters:
- E0: Baseline effect (no drug)
- Emax: Maximum effect above baseline
- EC50: Concentration at 50% of maximum effect
- gamma: Hill coefficient (steepness parameter)
  - gamma = 1: Standard Emax model (hyperbolic)
  - gamma > 1: Steeper (more switch-like) response
  - gamma < 1: More gradual response

Note: gamma is also known as the Hill coefficient or slope factor.
Typical range is 0.5 to 5.
"""
struct SigmoidEmax <: PDModelKind end

struct SigmoidEmaxParams
    E0::Float64      # Baseline effect
    Emax::Float64    # Maximum effect
    EC50::Float64    # Concentration at 50% Emax
    gamma::Float64   # Hill coefficient
end

"""
Biophase equilibration (effect compartment) PD model.

This model introduces a hypothetical effect compartment to account for
temporal delays between plasma concentration and observed effect.
The effect compartment has no volume (doesn't affect PK) and equilibrates
with the plasma concentration via first-order kinetics.

States:
- Ce: Effect site concentration (hypothetical)

Parameters:
- ke0: Effect site equilibration rate constant (1/time)
- E0: Baseline effect
- Emax: Maximum effect
- EC50: Effect site concentration at 50% Emax

Dynamics:
dCe/dt = ke0 * (Cp - Ce)

where Cp is plasma concentration (from PK model)

Effect:
E(Ce) = E0 + (Emax * Ce) / (EC50 + Ce)

Note: t1/2,ke0 = ln(2)/ke0 is the equilibration half-life.
When ke0 is large, effect follows plasma concentration closely (direct effect).
When ke0 is small, there is significant hysteresis between PK and PD.
"""
struct BiophaseEquilibration <: PDModelKind end

struct BiophaseEquilibrationParams
    ke0::Float64    # Effect site equilibration rate constant
    E0::Float64     # Baseline effect
    Emax::Float64   # Maximum effect
    EC50::Float64   # Effect site EC50
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
