# TMDD Analysis Utilities
#
# Pharmacokinetic analysis functions specific to TMDD models.
# Includes target occupancy, steady-state calculations, and model selection.

using Statistics

export calculate_target_occupancy, calculate_tmdd_steady_state
export calculate_tmdd_half_lives, identify_tmdd_regime
export tmdd_exposure_metrics, compare_tmdd_approximations

# ============================================================================
# Target Occupancy Analysis
# ============================================================================

"""
    calculate_target_occupancy(C, KD) -> Float64

Calculate target occupancy at a given drug concentration.

Target Occupancy = C / (KD + C)

# Arguments
- `C`: Drug concentration
- `KD`: Equilibrium dissociation constant

# Returns
- Fraction of target occupied (0 to 1)
"""
function calculate_target_occupancy(C::Float64, KD::Float64)::Float64
    if C <= 0.0 || KD <= 0.0
        return 0.0
    end
    return C / (KD + C)
end

"""
Calculate EC50 for target occupancy.
Returns KD (by definition, 50% occupancy occurs at C = KD).
"""
target_occupancy_EC50(KD::Float64) = KD

"""
Calculate EC90 for target occupancy (90% occupancy).
"""
target_occupancy_EC90(KD::Float64) = 9.0 * KD

"""
Calculate EC99 for target occupancy (99% occupancy).
"""
target_occupancy_EC99(KD::Float64) = 99.0 * KD

# ============================================================================
# Steady-State Calculations
# ============================================================================

"""
    calculate_tmdd_steady_state(params::TMDDFullParams) -> NamedTuple

Calculate steady-state conditions for full TMDD model.

Returns:
- R_ss: Free target at steady state (no drug)
- Rtot_ss: Total target at steady state (no drug)
- t_half_target: Target half-life
"""
function calculate_tmdd_steady_state(params::TMDDFullParams)
    # At steady state without drug
    R_ss = params.ksyn / params.kdeg
    Rtot_ss = R_ss
    t_half_target = log(2) / params.kdeg

    # KD for reference
    KD = params.koff / params.kon

    return (
        R_ss = R_ss,
        Rtot_ss = Rtot_ss,
        t_half_target = t_half_target,
        KD = KD
    )
end

function calculate_tmdd_steady_state(params::TMDD2CptQSSParams)
    R_ss = params.ksyn / params.kdeg
    Rtot_ss = R_ss
    t_half_target = log(2) / params.kdeg

    return (
        R_ss = R_ss,
        Rtot_ss = Rtot_ss,
        t_half_target = t_half_target,
        KSS = params.KSS
    )
end

function calculate_tmdd_steady_state(params::TMDD2CptCLParams)
    R_ss = params.R0
    t_half_target = log(2) / params.kdeg

    return (
        R_ss = R_ss,
        Rtot_ss = R_ss,
        t_half_target = t_half_target,
        Kss = params.Kss
    )
end

# ============================================================================
# Half-Life Calculations
# ============================================================================

"""
    calculate_tmdd_half_lives(params, C) -> NamedTuple

Calculate effective half-lives for TMDD system at different concentrations.

At high concentrations: t1/2 ≈ ln(2) / kel (linear elimination dominates)
At low concentrations: t1/2 much shorter (target-mediated elimination dominates)

Returns multiple half-life estimates for different concentration regimes.
"""
function calculate_tmdd_half_lives(params::TMDDFullParams, C_high::Float64)
    # Linear elimination half-life (high concentration)
    t_half_linear = log(2) / params.kel

    # Target half-life
    t_half_target = log(2) / params.kdeg

    # Complex half-life
    t_half_complex = log(2) / params.kint

    # Effective clearance at high vs low concentrations
    V = params.V
    CL_linear = params.kel * V

    # Target-mediated clearance (max)
    R0 = params.R0
    CL_target_max = params.kint * R0 / params.kon * V  # Approximate

    return (
        t_half_linear = t_half_linear,
        t_half_target = t_half_target,
        t_half_complex = t_half_complex,
        CL_linear = CL_linear
    )
end

function calculate_tmdd_half_lives(params::TMDD2CptQSSParams)
    # Two-compartment: distribution and terminal phases
    V1 = params.V1
    V2 = params.V2
    Q = params.Q
    kel = params.kel

    k12 = Q / V1
    k21 = Q / V2

    # Eigenvalues of the system (simplified, ignoring target)
    a = kel + k12
    b = k21
    discriminant = sqrt((a - b)^2 + 4 * k12 * k21)

    alpha = 0.5 * (a + b + discriminant)
    beta = 0.5 * (a + b - discriminant)

    t_half_alpha = log(2) / alpha  # Distribution phase
    t_half_beta = log(2) / beta    # Terminal phase

    return (
        t_half_alpha = t_half_alpha,
        t_half_beta = t_half_beta,
        t_half_target = log(2) / params.kdeg
    )
end

# ============================================================================
# TMDD Regime Identification
# ============================================================================

"""
    identify_tmdd_regime(C, KSS, R_total) -> Symbol

Identify the TMDD kinetic regime based on drug concentration.

Returns:
- :linear - High concentration, linear kinetics dominate
- :mixed - Intermediate, both linear and TMDD
- :tmdd - Low concentration, TMDD kinetics dominate
"""
function identify_tmdd_regime(C::Float64, KSS::Float64, R_total::Float64)
    # Target occupancy
    occupancy = C / (KSS + C)

    # Ratio of free drug to target
    ratio = C / R_total

    if ratio > 10.0 || occupancy > 0.95
        return :linear  # Linear kinetics dominate
    elseif ratio < 0.1 || occupancy < 0.5
        return :tmdd    # TMDD kinetics dominate
    else
        return :mixed   # Mixed kinetics
    end
end

# ============================================================================
# Exposure Metrics
# ============================================================================

"""
    tmdd_exposure_metrics(result::TMDDSimResult) -> NamedTuple

Calculate PK exposure metrics from TMDD simulation.

Returns:
- AUC: Area under concentration curve
- Cmax: Maximum concentration
- Tmax: Time of maximum concentration
- Ctrough: Concentration at end of simulation
- AUC_occupancy: Area under target occupancy curve
"""
function tmdd_exposure_metrics(result::TMDDSimResult)
    t = result.t
    conc = result.observations[:conc]

    # AUC using trapezoidal rule
    AUC = 0.0
    for i in 2:length(t)
        dt = t[i] - t[i-1]
        AUC += 0.5 * (conc[i] + conc[i-1]) * dt
    end

    # Cmax and Tmax
    Cmax_idx = argmax(conc)
    Cmax = conc[Cmax_idx]
    Tmax = t[Cmax_idx]

    # Ctrough (last concentration)
    Ctrough = conc[end]

    # Target occupancy AUC (if available)
    AUC_occupancy = 0.0
    if haskey(result.observations, :target_occupancy)
        occ = result.observations[:target_occupancy]
        for i in 2:length(t)
            dt = t[i] - t[i-1]
            AUC_occupancy += 0.5 * (occ[i] + occ[i-1]) * dt
        end
    end

    # Mean target occupancy
    mean_occupancy = haskey(result.observations, :target_occupancy) ?
        mean(result.observations[:target_occupancy]) : NaN

    return (
        AUC = AUC,
        Cmax = Cmax,
        Tmax = Tmax,
        Ctrough = Ctrough,
        AUC_occupancy = AUC_occupancy,
        mean_occupancy = mean_occupancy
    )
end

# ============================================================================
# Model Comparison/Selection
# ============================================================================

"""
    compare_tmdd_approximations(full_result, approx_result) -> NamedTuple

Compare full TMDD model to an approximation (QSS, QE, MM).

Returns metrics for assessing approximation validity.
"""
function compare_tmdd_approximations(
    full_result::TMDDSimResult,
    approx_result::TMDDSimResult
)
    # Ensure same time points
    t_full = full_result.t
    t_approx = approx_result.t

    if t_full != t_approx
        @warn "Time grids differ between models"
    end

    # Compare concentrations
    conc_full = full_result.observations[:conc]
    conc_approx = approx_result.observations[:conc]

    # Relative error
    rel_errors = [
        abs(conc_full[i] - conc_approx[i]) / max(conc_full[i], 1e-10)
        for i in 1:min(length(conc_full), length(conc_approx))
    ]

    max_rel_error = maximum(rel_errors)
    mean_rel_error = mean(rel_errors)

    # Compare AUC
    metrics_full = tmdd_exposure_metrics(full_result)
    metrics_approx = tmdd_exposure_metrics(approx_result)

    AUC_rel_error = abs(metrics_full.AUC - metrics_approx.AUC) / metrics_full.AUC

    # Recommendation
    if max_rel_error < 0.05
        recommendation = :excellent
    elseif max_rel_error < 0.15
        recommendation = :acceptable
    else
        recommendation = :use_full_model
    end

    return (
        max_rel_error = max_rel_error,
        mean_rel_error = mean_rel_error,
        AUC_rel_error = AUC_rel_error,
        recommendation = recommendation
    )
end

# ============================================================================
# Dose Selection Guidance
# ============================================================================

"""
    target_trough_dose(params, target_occupancy, dosing_interval) -> Float64

Estimate dose needed to achieve target trough occupancy.

Simplified calculation based on steady-state assumptions.
"""
function target_trough_dose(
    params::TMDD2CptQSSParams,
    target_occupancy::Float64,
    dosing_interval::Float64
)
    # Target concentration for desired occupancy
    # occupancy = C / (KSS + C)
    # C = KSS * occupancy / (1 - occupancy)
    C_target = params.KSS * target_occupancy / (1.0 - target_occupancy)

    # Very rough estimate assuming linear decline
    # C_trough = C_max * exp(-kel * tau)
    # We want C_trough = C_target
    # C_max = C_target * exp(kel * tau)
    C_max_needed = C_target * exp(params.kel * dosing_interval)

    # Dose ≈ Cmax * V (very simplified)
    dose = C_max_needed * params.V1

    return dose
end

# ============================================================================
# Time Above Target Metrics
# ============================================================================

"""
    time_above_occupancy(result::TMDDSimResult, threshold) -> Float64

Calculate time (hours) above a target occupancy threshold.
"""
function time_above_occupancy(result::TMDDSimResult, threshold::Float64)
    if !haskey(result.observations, :target_occupancy)
        return NaN
    end

    t = result.t
    occ = result.observations[:target_occupancy]

    time_above = 0.0
    for i in 2:length(t)
        if occ[i] >= threshold && occ[i-1] >= threshold
            time_above += t[i] - t[i-1]
        elseif occ[i] >= threshold || occ[i-1] >= threshold
            # Interpolate crossing point
            time_above += 0.5 * (t[i] - t[i-1])
        end
    end

    return time_above
end

"""
    time_above_concentration(result::TMDDSimResult, threshold) -> Float64

Calculate time (hours) above a concentration threshold.
"""
function time_above_concentration(result::TMDDSimResult, threshold::Float64)
    t = result.t
    conc = result.observations[:conc]

    time_above = 0.0
    for i in 2:length(t)
        if conc[i] >= threshold && conc[i-1] >= threshold
            time_above += t[i] - t[i-1]
        elseif conc[i] >= threshold || conc[i-1] >= threshold
            time_above += 0.5 * (t[i] - t[i-1])
        end
    end

    return time_above
end

export target_occupancy_EC50, target_occupancy_EC90, target_occupancy_EC99
export target_trough_dose, time_above_occupancy, time_above_concentration
