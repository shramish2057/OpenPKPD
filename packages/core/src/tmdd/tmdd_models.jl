# TMDD Model Implementations
#
# ODE systems and simulation functions for all TMDD model variants.
# Industry-standard implementations based on published literature.

using DifferentialEquations
using SciMLBase

# ============================================================================
# Full TMDD Model (Mager-Jusko)
# ============================================================================

"""
ODE system for full TMDD model (3 states: L, R, P).

States:
- u[1] = L (free drug amount)
- u[2] = R (free target amount)
- u[3] = P (drug-target complex amount)
"""
function _ode_tmdd_full!(du, u, p, t)
    L = u[1]  # Free drug
    R = u[2]  # Free target
    P = u[3]  # Complex

    kel = p.kel
    V = p.V
    kon = p.kon
    koff = p.koff
    ksyn = p.ksyn
    kdeg = p.kdeg
    kint = p.kint

    # Concentration of free drug
    C = L / V

    # Binding rate (mass action)
    binding_rate = kon * C * R

    # Free drug: elimination + binding - unbinding
    du[1] = -kel * L - binding_rate + koff * P

    # Free target: synthesis - degradation - binding + unbinding
    du[2] = ksyn - kdeg * R - binding_rate + koff * P

    # Complex: binding - unbinding - internalization
    du[3] = binding_rate - koff * P - kint * P

    return nothing
end

# ============================================================================
# QSS (Quasi-Steady-State) Approximation
# ============================================================================

"""
ODE system for QSS TMDD model (2 states: L, Rtot).

States:
- u[1] = L (free drug amount)
- u[2] = Rtot (total target amount = R + RL)
"""
function _ode_tmdd_qss!(du, u, p, t)
    L = u[1]     # Free drug
    Rtot = u[2]  # Total target

    kel = p.kel
    V = p.V
    KSS = p.KSS
    ksyn = p.ksyn
    kdeg = p.kdeg
    kint = p.kint

    # QSS: Complex at steady state
    # RL = L * Rtot / (KSS*V + L)
    KSS_V = KSS * V
    RL = L * Rtot / (KSS_V + L)
    R = Rtot - RL

    # Ensure non-negative
    RL = max(0.0, RL)
    R = max(0.0, R)

    # Free drug: elimination + complex internalization
    du[1] = -kel * L - kint * RL

    # Total target: synthesis - degradation of free - internalization of complex
    du[2] = ksyn - kdeg * R - kint * RL

    return nothing
end

# ============================================================================
# QE (Quasi-Equilibrium) Approximation
# ============================================================================

"""
ODE system for QE TMDD model (2 states: L, Rtot).

States:
- u[1] = L (free drug amount)
- u[2] = Rtot (total target amount)
"""
function _ode_tmdd_qe!(du, u, p, t)
    L = u[1]     # Free drug
    Rtot = u[2]  # Total target

    kel = p.kel
    V = p.V
    KD = p.KD
    ksyn = p.ksyn
    kdeg = p.kdeg  # Used for both R and RL (QE assumption)

    # QE: Complex at equilibrium
    # RL = L * Rtot / (KD*V + L)
    KD_V = KD * V
    RL = L * Rtot / (KD_V + L)

    # Ensure non-negative
    RL = max(0.0, RL)

    # Free drug: elimination + complex elimination
    du[1] = -kel * L - kdeg * RL

    # Total target: synthesis - degradation (same rate for both forms in QE)
    du[2] = ksyn - kdeg * Rtot

    return nothing
end

# ============================================================================
# Michaelis-Menten Approximation
# ============================================================================

"""
ODE system for MM TMDD approximation (1 state: L).

States:
- u[1] = L (drug amount)
"""
function _ode_tmdd_mm!(du, u, p, t)
    L = u[1]

    kel = p.kel
    V = p.V
    Vmax = p.Vmax
    Km = p.Km

    # Concentration
    C = L / V

    # Parallel linear and saturable elimination
    du[1] = -kel * L - (Vmax * C) / (Km + C)

    return nothing
end

# ============================================================================
# Rapid Binding (Wagner) Approximation
# ============================================================================

"""
ODE system for rapid binding TMDD model (2 states: Ltot, Rtot).

States:
- u[1] = Ltot (total drug amount = L + P)
- u[2] = Rtot (total target amount = R + P)
"""
function _ode_tmdd_rapid_binding!(du, u, p, t)
    Ltot = u[1]  # Total drug
    Rtot = u[2]  # Total target

    kel = p.kel
    V = p.V
    KD = p.KD
    ksyn = p.ksyn
    kdeg = p.kdeg
    kint = p.kint

    # Calculate free drug from quadratic solution
    # L = 0.5 * ((Ltot - Rtot - KD*V) + sqrt((Ltot - Rtot - KD*V)^2 + 4*KD*V*Ltot))
    KD_V = KD * V
    a = Ltot - Rtot - KD_V
    discriminant = a^2 + 4 * KD_V * Ltot
    L = 0.5 * (a + sqrt(max(0.0, discriminant)))
    L = max(0.0, L)

    # Complex and free target from conservation
    P = Ltot - L
    R = Rtot - P
    P = max(0.0, P)
    R = max(0.0, R)

    # Total drug: free drug elimination + complex internalization
    du[1] = -kel * L - kint * P

    # Total target: synthesis - degradation - internalization
    du[2] = ksyn - kdeg * R - kint * P

    return nothing
end

# ============================================================================
# Irreversible Binding
# ============================================================================

"""
ODE system for irreversible binding TMDD (2 states: L, R).

States:
- u[1] = L (free drug amount)
- u[2] = R (free target amount)
"""
function _ode_tmdd_irreversible!(du, u, p, t)
    L = u[1]  # Free drug
    R = u[2]  # Free target

    kel = p.kel
    V = p.V
    kon = p.kon
    ksyn = p.ksyn
    kdeg = p.kdeg

    # Concentration
    C = L / V

    # Irreversible binding (koff = 0)
    binding_rate = kon * C * R

    # Free drug: elimination + irreversible binding
    du[1] = -kel * L - binding_rate

    # Free target: synthesis - degradation - irreversible binding
    du[2] = ksyn - kdeg * R - binding_rate

    return nothing
end

# ============================================================================
# Two-Compartment Full TMDD
# ============================================================================

"""
ODE system for 2-compartment full TMDD (4 states: L, Lp, R, P).

States:
- u[1] = L (free drug in central)
- u[2] = Lp (drug in peripheral)
- u[3] = R (free target)
- u[4] = P (drug-target complex)
"""
function _ode_tmdd_2cpt_full!(du, u, p, t)
    L = u[1]   # Free drug (central)
    Lp = u[2]  # Drug (peripheral)
    R = u[3]   # Free target
    P = u[4]   # Complex

    kel = p.kel
    V1 = p.V1
    V2 = p.V2
    Q = p.Q
    kon = p.kon
    koff = p.koff
    ksyn = p.ksyn
    kdeg = p.kdeg
    kint = p.kint

    # Micro-constants
    k12 = Q / V1
    k21 = Q / V2

    # Concentration
    C = L / V1

    # Binding rate
    binding_rate = kon * C * R

    # Central drug: elimination + binding/unbinding + distribution
    du[1] = -kel * L - binding_rate + koff * P - k12 * L + k21 * Lp

    # Peripheral drug: distribution only
    du[2] = k12 * L - k21 * Lp

    # Free target: synthesis - degradation - binding + unbinding
    du[3] = ksyn - kdeg * R - binding_rate + koff * P

    # Complex: binding - unbinding - internalization
    du[4] = binding_rate - koff * P - kint * P

    return nothing
end

# ============================================================================
# Two-Compartment QSS TMDD
# ============================================================================

"""
ODE system for 2-compartment QSS TMDD (3 states: L, Lp, Rtot).

States:
- u[1] = L (free drug in central)
- u[2] = Lp (drug in peripheral)
- u[3] = Rtot (total target)
"""
function _ode_tmdd_2cpt_qss!(du, u, p, t)
    L = u[1]     # Free drug (central)
    Lp = u[2]    # Drug (peripheral)
    Rtot = u[3]  # Total target

    kel = p.kel
    V1 = p.V1
    V2 = p.V2
    Q = p.Q
    KSS = p.KSS
    ksyn = p.ksyn
    kdeg = p.kdeg
    kint = p.kint

    # Micro-constants
    k12 = Q / V1
    k21 = Q / V2

    # QSS: Complex at steady state
    KSS_V1 = KSS * V1
    RL = L * Rtot / (KSS_V1 + L)
    R = Rtot - RL

    # Ensure non-negative
    RL = max(0.0, RL)
    R = max(0.0, R)

    # Central drug: elimination + internalization + distribution
    du[1] = -kel * L - kint * RL - k12 * L + k21 * Lp

    # Peripheral drug: distribution only
    du[2] = k12 * L - k21 * Lp

    # Total target: synthesis - degradation - internalization
    du[3] = ksyn - kdeg * R - kint * RL

    return nothing
end

# ============================================================================
# Two-Compartment Clearance Parameterization
# ============================================================================

"""
ODE system for 2-compartment TMDD with clearance parameterization.

States:
- u[1] = L (free drug in central)
- u[2] = Lp (drug in peripheral)
- u[3] = Rtot (total target)
"""
function _ode_tmdd_2cpt_cl!(du, u, p, t)
    L = u[1]     # Free drug (central)
    Lp = u[2]    # Drug (peripheral)
    Rtot = u[3]  # Total target

    CL = p.CL
    V1 = p.V1
    V2 = p.V2
    Q = p.Q
    Kss = p.Kss
    Vmax = p.Vmax
    R0 = p.R0
    kdeg = p.kdeg

    # Derived parameters
    kel = CL / V1
    kint = Vmax / R0
    ksyn = kdeg * R0

    # Micro-constants
    k12 = Q / V1
    k21 = Q / V2

    # QSS: Complex at steady state
    Kss_V1 = Kss * V1
    RL = L * Rtot / (Kss_V1 + L)
    R = Rtot - RL

    # Ensure non-negative
    RL = max(0.0, RL)
    R = max(0.0, R)

    # Central drug: elimination + internalization + distribution
    du[1] = -kel * L - kint * RL - k12 * L + k21 * Lp

    # Peripheral drug: distribution only
    du[2] = k12 * L - k21 * Lp

    # Total target: synthesis - degradation - internalization
    du[3] = ksyn - kdeg * R - kint * RL

    return nothing
end

# ============================================================================
# Soluble Target TMDD
# ============================================================================

"""
ODE system for soluble target TMDD (3 states: L, R, P).

States:
- u[1] = L (free drug)
- u[2] = R (free soluble target)
- u[3] = P (drug-target complex)
"""
function _ode_tmdd_soluble!(du, u, p, t)
    L = u[1]  # Free drug
    R = u[2]  # Free target
    P = u[3]  # Complex

    kel = p.kel
    V = p.V
    kon = p.kon
    koff = p.koff
    ksyn = p.ksyn
    kdeg = p.kdeg
    kel_complex = p.kel_complex

    # Concentration
    C = L / V

    # Binding rate
    binding_rate = kon * C * R

    # Free drug: elimination + binding - unbinding
    du[1] = -kel * L - binding_rate + koff * P

    # Free target: synthesis - degradation - binding + unbinding
    du[2] = ksyn - kdeg * R - binding_rate + koff * P

    # Complex: binding - unbinding - elimination (different rate)
    du[3] = binding_rate - koff * P - kel_complex * P

    return nothing
end

# ============================================================================
# Internalization Model
# ============================================================================

"""
ODE system for TMDD with internalization/recycling (5 states).

States:
- u[1] = L (free drug)
- u[2] = R (free surface receptor)
- u[3] = P (surface drug-receptor complex)
- u[4] = Re (endosomal free receptor)
- u[5] = Pe (endosomal complex)
"""
function _ode_tmdd_internalization!(du, u, p, t)
    L = u[1]   # Free drug
    R = u[2]   # Surface receptor
    P = u[3]   # Surface complex
    Re = u[4]  # Endosomal receptor
    Pe = u[5]  # Endosomal complex

    kel = p.kel
    V = p.V
    kon = p.kon
    koff = p.koff
    ksyn = p.ksyn
    kint = p.kint
    kendo_R = p.kendo_R
    krec = p.krec
    kdeg_e = p.kdeg_e

    # Concentration
    C = L / V

    # Binding rate
    binding_rate = kon * C * R

    # Free drug: elimination + binding - unbinding
    du[1] = -kel * L - binding_rate + koff * P

    # Surface receptor: synthesis + recycling - binding + unbinding - endocytosis
    du[2] = ksyn + krec * Re - binding_rate + koff * P - kendo_R * R

    # Surface complex: binding - unbinding - internalization
    du[3] = binding_rate - koff * P - kint * P

    # Endosomal receptor: endocytosis - degradation - recycling
    du[4] = kendo_R * R - kdeg_e * Re - krec * Re

    # Endosomal complex: internalization - degradation
    du[5] = kint * P - kdeg_e * Pe

    return nothing
end

# ============================================================================
# State and Observation Extractors
# ============================================================================

"""
Get the number of ODE states for a TMDD model type.
"""
n_states(::TMDDFull) = 3
n_states(::TMDDQSS) = 2
n_states(::TMDDQE) = 2
n_states(::TMDDMM) = 1
n_states(::TMDDRapidBinding) = 2
n_states(::TMDDIrreversible) = 2
n_states(::TMDD2CptFull) = 4
n_states(::TMDD2CptQSS) = 3
n_states(::TMDD2CptCL) = 3
n_states(::TMDDSolubleTarget) = 3
n_states(::TMDDInternalization) = 5

"""
Get initial conditions for a TMDD model.
"""
function get_initial_state(spec::TMDDSpec{TMDDFull,TMDDFullParams})
    p = spec.params
    return [0.0, p.R0, 0.0]  # L=0, R=R0, P=0
end

function get_initial_state(spec::TMDDSpec{TMDDQSS,TMDDQSSParams})
    p = spec.params
    return [0.0, p.Rtot0]  # L=0, Rtot=Rtot0
end

function get_initial_state(spec::TMDDSpec{TMDDQE,TMDDQEParams})
    p = spec.params
    return [0.0, p.Rtot0]  # L=0, Rtot=Rtot0
end

function get_initial_state(spec::TMDDSpec{TMDDMM,TMDDMMParams})
    return [0.0]  # L=0
end

function get_initial_state(spec::TMDDSpec{TMDDRapidBinding,TMDDRapidBindingParams})
    p = spec.params
    return [0.0, p.Rtot0]  # Ltot=0, Rtot=Rtot0
end

function get_initial_state(spec::TMDDSpec{TMDDIrreversible,TMDDIrreversibleParams})
    p = spec.params
    return [0.0, p.R0]  # L=0, R=R0
end

function get_initial_state(spec::TMDDSpec{TMDD2CptFull,TMDD2CptFullParams})
    p = spec.params
    return [0.0, 0.0, p.R0, 0.0]  # L=0, Lp=0, R=R0, P=0
end

function get_initial_state(spec::TMDDSpec{TMDD2CptQSS,TMDD2CptQSSParams})
    p = spec.params
    return [0.0, 0.0, p.Rtot0]  # L=0, Lp=0, Rtot=Rtot0
end

function get_initial_state(spec::TMDDSpec{TMDD2CptCL,TMDD2CptCLParams})
    p = spec.params
    return [0.0, 0.0, p.R0]  # L=0, Lp=0, Rtot=R0
end

function get_initial_state(spec::TMDDSpec{TMDDSolubleTarget,TMDDSolubleTargetParams})
    p = spec.params
    return [0.0, p.R0, 0.0]  # L=0, R=R0, P=0
end

function get_initial_state(spec::TMDDSpec{TMDDInternalization,TMDDInternalizationParams})
    p = spec.params
    Re0 = p.kendo_R * p.R0 / (p.kdeg_e + p.krec)  # Steady state endosomal
    return [0.0, p.R0, 0.0, Re0, 0.0]  # L=0, R=R0, P=0, Re=steady-state, Pe=0
end

"""
Get ODE function for a TMDD model type.
"""
get_ode_function(::TMDDFull) = _ode_tmdd_full!
get_ode_function(::TMDDQSS) = _ode_tmdd_qss!
get_ode_function(::TMDDQE) = _ode_tmdd_qe!
get_ode_function(::TMDDMM) = _ode_tmdd_mm!
get_ode_function(::TMDDRapidBinding) = _ode_tmdd_rapid_binding!
get_ode_function(::TMDDIrreversible) = _ode_tmdd_irreversible!
get_ode_function(::TMDD2CptFull) = _ode_tmdd_2cpt_full!
get_ode_function(::TMDD2CptQSS) = _ode_tmdd_2cpt_qss!
get_ode_function(::TMDD2CptCL) = _ode_tmdd_2cpt_cl!
get_ode_function(::TMDDSolubleTarget) = _ode_tmdd_soluble!
get_ode_function(::TMDDInternalization) = _ode_tmdd_internalization!

"""
Get the dosing compartment index for a TMDD model.
"""
dosing_compartment(::TMDDFull) = 1
dosing_compartment(::TMDDQSS) = 1
dosing_compartment(::TMDDQE) = 1
dosing_compartment(::TMDDMM) = 1
dosing_compartment(::TMDDRapidBinding) = 1
dosing_compartment(::TMDDIrreversible) = 1
dosing_compartment(::TMDD2CptFull) = 1
dosing_compartment(::TMDD2CptQSS) = 1
dosing_compartment(::TMDD2CptCL) = 1
dosing_compartment(::TMDDSolubleTarget) = 1
dosing_compartment(::TMDDInternalization) = 1

"""
Get state names for a TMDD model type.
"""
state_names(::TMDDFull) = [:L, :R, :P]
state_names(::TMDDQSS) = [:L, :Rtot]
state_names(::TMDDQE) = [:L, :Rtot]
state_names(::TMDDMM) = [:L]
state_names(::TMDDRapidBinding) = [:Ltot, :Rtot]
state_names(::TMDDIrreversible) = [:L, :R]
state_names(::TMDD2CptFull) = [:L, :Lp, :R, :P]
state_names(::TMDD2CptQSS) = [:L, :Lp, :Rtot]
state_names(::TMDD2CptCL) = [:L, :Lp, :Rtot]
state_names(::TMDDSolubleTarget) = [:L, :R, :P]
state_names(::TMDDInternalization) = [:L, :R, :P, :Re, :Pe]

"""
Get central volume from TMDD parameters.
"""
get_central_volume(p::TMDDFullParams) = p.V
get_central_volume(p::TMDDQSSParams) = p.V
get_central_volume(p::TMDDQEParams) = p.V
get_central_volume(p::TMDDMMParams) = p.V
get_central_volume(p::TMDDRapidBindingParams) = p.V
get_central_volume(p::TMDDIrreversibleParams) = p.V
get_central_volume(p::TMDD2CptFullParams) = p.V1
get_central_volume(p::TMDD2CptQSSParams) = p.V1
get_central_volume(p::TMDD2CptCLParams) = p.V1
get_central_volume(p::TMDDSolubleTargetParams) = p.V
get_central_volume(p::TMDDInternalizationParams) = p.V

export n_states, get_initial_state, get_ode_function, dosing_compartment
export state_names, get_central_volume
