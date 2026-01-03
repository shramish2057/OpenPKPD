using StableRNGs
using Distributions

export PopulationResult, PopulationSummary, simulate_population

"""
Population summary statistics for a single observation series.

Each vector is aligned with time grid t.
"""
struct PopulationSummary
    observation::Symbol
    probs::Vector{Float64}
    mean::Vector{Float64}
    median::Vector{Float64}
    quantiles::Dict{Float64,Vector{Float64}}
end

"""
Population simulation output.

individuals:
- Vector of SimResult, one per individual

params:
- Vector of Dicts describing realized individual parameters

metadata:
- includes seed, n, iiv kind, and omega map
"""
struct PopulationResult
    individuals::Vector{SimResult}
    params::Vector{Dict{Symbol,Float64}}
    summaries::Dict{Symbol,PopulationSummary}
    metadata::Dict{String,Any}
end

function validate(iiv::IIVSpec{LogNormalIIV})
    if iiv.n < 1
        error("IIVSpec n must be >= 1")
    end
    for (k, ω) in iiv.omegas
        _require_positive("omega for $(k)", ω)
    end
    return nothing
end

function _realize_params_log_normal(base_params, omegas::Dict{Symbol,Float64}, rng)
    # base_params is a typed struct, we produce a new typed struct with per-individual values
    T = typeof(base_params)
    fn = fieldnames(T)

    vals = Dict{Symbol,Float64}()

    for f in fn
        θ = Float64(getfield(base_params, f))
        if haskey(omegas, f)
            ω = omegas[f]
            η = rand(rng, Normal(0.0, ω))
            vals[f] = θ * exp(η)
        else
            vals[f] = θ
        end
    end

    # Reconstruct typed params in declared field order
    new_params = T((vals[f] for f in fn)...)

    return new_params, vals
end

function _mean(xs::Vector{Float64})
    s = 0.0
    for x in xs
        s += x
    end
    return s / length(xs)
end

function _median_sorted(xs_sorted::Vector{Float64})
    n = length(xs_sorted)
    if isodd(n)
        return xs_sorted[(n + 1) ÷ 2]
    end
    return 0.5 * (xs_sorted[n ÷ 2] + xs_sorted[n ÷ 2 + 1])
end

function _quantile_sorted(xs_sorted::Vector{Float64}, p::Float64)
    # Linear interpolation between order statistics, deterministic.
    n = length(xs_sorted)
    if n == 1
        return xs_sorted[1]
    end
    if p <= 0.0
        return xs_sorted[1]
    end
    if p >= 1.0
        return xs_sorted[end]
    end

    h = 1.0 + (n - 1) * p
    lo = floor(Int, h)
    hi = ceil(Int, h)

    if lo == hi
        return xs_sorted[lo]
    end

    w = h - lo
    return (1.0 - w) * xs_sorted[lo] + w * xs_sorted[hi]
end

function compute_population_summary(
    individuals::Vector{SimResult}, observation::Symbol; probs::Vector{Float64}=[0.05, 0.95]
)
    if isempty(individuals)
        error("Cannot compute summary for empty individuals")
    end

    t = individuals[1].t
    n_t = length(t)

    for (i, ind) in enumerate(individuals)
        if ind.t != t
            error(
                "All individuals must share identical time grid for summaries. Mismatch at index $(i)",
            )
        end
        if !haskey(ind.observations, observation)
            error("Observation $(observation) missing for individual $(i)")
        end
    end

    mean_v = Vector{Float64}(undef, n_t)
    median_v = Vector{Float64}(undef, n_t)
    qmap = Dict{Float64,Vector{Float64}}(p => Vector{Float64}(undef, n_t) for p in probs)

    scratch = Vector{Float64}(undef, length(individuals))

    for ti in 1:n_t
        for (j, ind) in enumerate(individuals)
            scratch[j] = ind.observations[observation][ti]
        end

        sort!(scratch)

        mean_v[ti] = _mean(scratch)
        median_v[ti] = _median_sorted(scratch)

        for p in probs
            qmap[p][ti] = _quantile_sorted(scratch, p)
        end
    end

    return PopulationSummary(observation, probs, mean_v, median_v, qmap)
end

"""
simulate_population runs deterministic population simulations for supported PK models.

Current support:
- OneCompIVBolusParams: CL, V
- OneCompOralFirstOrderParams: Ka, CL, V

IIV is log-normal, per-parameter independent (diagonal omega) in v1.

Covariates are accepted but not used yet, kept for forward compatibility.
"""
function simulate_population(pop::PopulationSpec, grid::SimGrid, solver::SolverSpec)
    base = pop.base_model_spec

    if pop.iiv === nothing
        n = 1
        seed = UInt64(0)
        iiv_kind = "none"
        omegas = Dict{Symbol,Float64}()
        rng = StableRNG(0)
    else
        validate(pop.iiv)
        n = pop.iiv.n
        seed = pop.iiv.seed
        iiv_kind = string(typeof(pop.iiv.kind))
        omegas = pop.iiv.omegas
        rng = StableRNG(seed)
    end

    # Covariates vector must be empty or length n
    if !isempty(pop.covariates) && length(pop.covariates) != n
        error("covariates must be empty or length n")
    end

    individuals = Vector{SimResult}(undef, n)
    realized_params = Vector{Dict{Symbol,Float64}}(undef, n)

    for i in 1:n
        if pop.iiv === nothing
            spec_i = base
            individuals[i] = simulate(spec_i, grid, solver)

            T = typeof(base.params)
            fn = fieldnames(T)
            d = Dict{Symbol,Float64}()
            for f in fn
                d[f] = Float64(getfield(base.params, f))
            end
            realized_params[i] = d
        else
            new_params, d = _realize_params_log_normal(base.params, omegas, rng)
            realized_params[i] = d

            spec_i = ModelSpec(base.kind, base.name * "_i$(i)", new_params, base.doses)
            individuals[i] = simulate(spec_i, grid, solver)
        end
    end

    # ✅ Now summaries can be computed safely
    summaries = Dict{Symbol,PopulationSummary}()

    if haskey(individuals[1].observations, :conc)
        summaries[:conc] = compute_population_summary(
            individuals, :conc; probs=[0.05, 0.95]
        )
    end

    metadata = Dict{String,Any}(
        "n" => n,
        "seed" => seed,
        "iiv_kind" => iiv_kind,
        "omegas" => Dict(String(k) => v for (k, v) in omegas),
        "engine_version" => "0.1.0",
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return PopulationResult(individuals, realized_params, summaries, metadata)
end
