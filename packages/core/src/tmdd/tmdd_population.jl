# TMDD Population Simulation
#
# Population-level simulation for TMDD models with IIV support.

using StableRNGs
using Statistics

export TMDDPopulationSpec, TMDDPopulationResult
export simulate_tmdd_population

# ============================================================================
# Population Specification
# ============================================================================

"""
Population specification for TMDD models.

Supports inter-individual variability (IIV) on TMDD parameters.
"""
struct TMDDPopulationSpec{K<:TMDDModelKind,P}
    base_spec::TMDDSpec{K,P}
    iiv::Union{Nothing,IIVSpec}
    n_subjects::Int
end

"""
Create a population spec from a base TMDD spec.
"""
function TMDDPopulationSpec(
    spec::TMDDSpec{K,P},
    iiv::IIVSpec
) where {K<:TMDDModelKind,P}
    return TMDDPopulationSpec(spec, iiv, iiv.n)
end

# ============================================================================
# Population Result
# ============================================================================

"""
Result container for TMDD population simulation.
"""
struct TMDDPopulationResult
    individual_results::Vector{TMDDSimResult}
    summary::Dict{Symbol,Dict{String,Vector{Float64}}}
    parameters::Dict{Symbol,Vector{Float64}}
    metadata::Dict{String,Any}
end

# ============================================================================
# Population Simulation
# ============================================================================

"""
    simulate_tmdd_population(pop_spec, grid, solver) -> TMDDPopulationResult

Simulate TMDD model for a population of subjects with IIV.

# Arguments
- `pop_spec::TMDDPopulationSpec`: Population specification
- `grid::SimGrid`: Simulation time grid
- `solver::SolverSpec`: ODE solver configuration

# Returns
- `TMDDPopulationResult`: Individual and summary results
"""
function simulate_tmdd_population(
    pop_spec::TMDDPopulationSpec{K,P},
    grid::SimGrid,
    solver::SolverSpec
) where {K<:TMDDModelKind,P}

    base_spec = pop_spec.base_spec
    n_subjects = pop_spec.n_subjects

    # Generate individual parameters
    individual_params = _generate_tmdd_individual_params(pop_spec)

    # Simulate each individual
    individual_results = Vector{TMDDSimResult}(undef, n_subjects)

    for i in 1:n_subjects
        # Create individual spec with perturbed parameters
        ind_spec = TMDDSpec(
            base_spec.kind,
            "$(base_spec.name)_ID$(i)",
            individual_params[i],
            base_spec.doses
        )

        # Simulate
        individual_results[i] = solve_tmdd(ind_spec, grid, solver)
    end

    # Compute summary statistics
    summary = _compute_population_summary(individual_results)

    # Collect parameter values
    param_values = _collect_param_values(individual_params)

    # Metadata
    metadata = Dict{String,Any}(
        "model_type" => string(K),
        "n_subjects" => n_subjects,
        "has_iiv" => pop_spec.iiv !== nothing
    )

    return TMDDPopulationResult(individual_results, summary, param_values, metadata)
end

"""
Generate individual TMDD parameters with IIV.
"""
function _generate_tmdd_individual_params(
    pop_spec::TMDDPopulationSpec{TMDDFull,TMDDFullParams}
)
    base_params = pop_spec.base_spec.params
    n = pop_spec.n_subjects
    iiv = pop_spec.iiv

    if iiv === nothing
        # No IIV - return copies of base params
        return [base_params for _ in 1:n]
    end

    rng = StableRNG(iiv.seed)
    omegas = iiv.omegas

    params = Vector{TMDDFullParams}(undef, n)

    for i in 1:n
        # Generate eta values
        kel = base_params.kel * exp(get(omegas, :kel, 0.0) * randn(rng))
        V = base_params.V * exp(get(omegas, :V, 0.0) * randn(rng))
        kon = base_params.kon * exp(get(omegas, :kon, 0.0) * randn(rng))
        koff = base_params.koff * exp(get(omegas, :koff, 0.0) * randn(rng))
        ksyn = base_params.ksyn * exp(get(omegas, :ksyn, 0.0) * randn(rng))
        kdeg = base_params.kdeg * exp(get(omegas, :kdeg, 0.0) * randn(rng))
        kint = base_params.kint * exp(get(omegas, :kint, 0.0) * randn(rng))
        R0 = base_params.R0 * exp(get(omegas, :R0, 0.0) * randn(rng))

        params[i] = TMDDFullParams(kel, V, kon, koff, ksyn, kdeg, kint, R0)
    end

    return params
end

function _generate_tmdd_individual_params(
    pop_spec::TMDDPopulationSpec{TMDD2CptQSS,TMDD2CptQSSParams}
)
    base_params = pop_spec.base_spec.params
    n = pop_spec.n_subjects
    iiv = pop_spec.iiv

    if iiv === nothing
        return [base_params for _ in 1:n]
    end

    rng = StableRNG(iiv.seed)
    omegas = iiv.omegas

    params = Vector{TMDD2CptQSSParams}(undef, n)

    for i in 1:n
        kel = base_params.kel * exp(get(omegas, :kel, 0.0) * randn(rng))
        V1 = base_params.V1 * exp(get(omegas, :V1, 0.0) * randn(rng))
        V2 = base_params.V2 * exp(get(omegas, :V2, 0.0) * randn(rng))
        Q = base_params.Q * exp(get(omegas, :Q, 0.0) * randn(rng))
        KSS = base_params.KSS * exp(get(omegas, :KSS, 0.0) * randn(rng))
        ksyn = base_params.ksyn * exp(get(omegas, :ksyn, 0.0) * randn(rng))
        kdeg = base_params.kdeg * exp(get(omegas, :kdeg, 0.0) * randn(rng))
        kint = base_params.kint * exp(get(omegas, :kint, 0.0) * randn(rng))
        Rtot0 = base_params.Rtot0 * exp(get(omegas, :Rtot0, 0.0) * randn(rng))

        params[i] = TMDD2CptQSSParams(kel, V1, V2, Q, KSS, ksyn, kdeg, kint, Rtot0)
    end

    return params
end

function _generate_tmdd_individual_params(
    pop_spec::TMDDPopulationSpec{TMDD2CptCL,TMDD2CptCLParams}
)
    base_params = pop_spec.base_spec.params
    n = pop_spec.n_subjects
    iiv = pop_spec.iiv

    if iiv === nothing
        return [base_params for _ in 1:n]
    end

    rng = StableRNG(iiv.seed)
    omegas = iiv.omegas

    params = Vector{TMDD2CptCLParams}(undef, n)

    for i in 1:n
        CL = base_params.CL * exp(get(omegas, :CL, 0.0) * randn(rng))
        V1 = base_params.V1 * exp(get(omegas, :V1, 0.0) * randn(rng))
        V2 = base_params.V2 * exp(get(omegas, :V2, 0.0) * randn(rng))
        Q = base_params.Q * exp(get(omegas, :Q, 0.0) * randn(rng))
        Kss = base_params.Kss * exp(get(omegas, :Kss, 0.0) * randn(rng))
        Vmax = base_params.Vmax * exp(get(omegas, :Vmax, 0.0) * randn(rng))
        R0 = base_params.R0 * exp(get(omegas, :R0, 0.0) * randn(rng))
        kdeg = base_params.kdeg * exp(get(omegas, :kdeg, 0.0) * randn(rng))

        params[i] = TMDD2CptCLParams(CL, V1, V2, Q, Kss, Vmax, R0, kdeg)
    end

    return params
end

# Generic fallback for other model types
function _generate_tmdd_individual_params(pop_spec::TMDDPopulationSpec)
    base_params = pop_spec.base_spec.params
    n = pop_spec.n_subjects
    return [base_params for _ in 1:n]
end

"""
Compute summary statistics across population.
"""
function _compute_population_summary(results::Vector{TMDDSimResult})
    if isempty(results)
        return Dict{Symbol,Dict{String,Vector{Float64}}}()
    end

    t = results[1].t
    n_times = length(t)
    n_subjects = length(results)

    summary = Dict{Symbol,Dict{String,Vector{Float64}}}()

    # Get all observation keys from first result
    obs_keys = keys(results[1].observations)

    for key in obs_keys
        # Collect values across subjects
        values_matrix = zeros(n_subjects, n_times)
        for (i, res) in enumerate(results)
            if haskey(res.observations, key)
                values_matrix[i, :] = res.observations[key]
            end
        end

        # Compute statistics
        mean_values = vec(mean(values_matrix, dims=1))
        std_values = vec(std(values_matrix, dims=1))
        median_values = vec(median(values_matrix, dims=1))

        # Percentiles
        p5 = [quantile(values_matrix[:, j], 0.05) for j in 1:n_times]
        p25 = [quantile(values_matrix[:, j], 0.25) for j in 1:n_times]
        p75 = [quantile(values_matrix[:, j], 0.75) for j in 1:n_times]
        p95 = [quantile(values_matrix[:, j], 0.95) for j in 1:n_times]

        summary[key] = Dict{String,Vector{Float64}}(
            "mean" => mean_values,
            "std" => std_values,
            "median" => median_values,
            "p5" => p5,
            "p25" => p25,
            "p75" => p75,
            "p95" => p95
        )
    end

    return summary
end

"""
Collect parameter values from individual params for analysis.
"""
function _collect_param_values(params::Vector{TMDDFullParams})
    n = length(params)
    return Dict{Symbol,Vector{Float64}}(
        :kel => [p.kel for p in params],
        :V => [p.V for p in params],
        :kon => [p.kon for p in params],
        :koff => [p.koff for p in params],
        :ksyn => [p.ksyn for p in params],
        :kdeg => [p.kdeg for p in params],
        :kint => [p.kint for p in params],
        :R0 => [p.R0 for p in params]
    )
end

function _collect_param_values(params::Vector{TMDD2CptQSSParams})
    n = length(params)
    return Dict{Symbol,Vector{Float64}}(
        :kel => [p.kel for p in params],
        :V1 => [p.V1 for p in params],
        :V2 => [p.V2 for p in params],
        :Q => [p.Q for p in params],
        :KSS => [p.KSS for p in params],
        :ksyn => [p.ksyn for p in params],
        :kdeg => [p.kdeg for p in params],
        :kint => [p.kint for p in params],
        :Rtot0 => [p.Rtot0 for p in params]
    )
end

function _collect_param_values(params::Vector{TMDD2CptCLParams})
    n = length(params)
    return Dict{Symbol,Vector{Float64}}(
        :CL => [p.CL for p in params],
        :V1 => [p.V1 for p in params],
        :V2 => [p.V2 for p in params],
        :Q => [p.Q for p in params],
        :Kss => [p.Kss for p in params],
        :Vmax => [p.Vmax for p in params],
        :R0 => [p.R0 for p in params],
        :kdeg => [p.kdeg for p in params]
    )
end

# Generic fallback
function _collect_param_values(params::Vector)
    return Dict{Symbol,Vector{Float64}}()
end
