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

function simulate_population(
    pop::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec;
    pd_spec::Union{Nothing,PDSpec}=nothing,
)
    base = pop.base_model_spec

    # Determine population size and RNG setup for IIV
    if pop.iiv === nothing
        if pop.covariate_model !== nothing
            n = length(pop.covariates)
            if n < 1
                error("covariates must be non-empty when covariate_model is provided")
            end
        else
            n = 1
        end

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

    if pop.iov !== nothing
        validate(pop.iov)
    end

    # Covariates length checks
    if pop.covariate_model !== nothing
        if length(pop.covariates) != n
            error("covariates must be length n when covariate_model is provided")
        end
    else
        if !isempty(pop.covariates) && length(pop.covariates) != n
            error("covariates must be empty or length n")
        end
    end

    individuals = Vector{SimResult}(undef, n)
    realized_params = Vector{Dict{Symbol,Float64}}(undef, n)

    # Precompute occasions once (same for all individuals) if IOV is enabled
    occ_starts = nothing
    if pop.iov !== nothing
        occ_starts = derive_occasions(base.doses, grid, pop.iov.occasion_def)
    end

    for i in 1:n
        # 1) Start from base params
        params_i = base.params

        # 2) Apply covariates (optional) and create a snapshot dict
        snap = Dict{Symbol,Float64}()

        if pop.covariate_model !== nothing
            covs = pop.covariates[i]
            params_i, snap = apply_covariates(base.params, pop.covariate_model, covs)
        else
            T = typeof(base.params)
            for f in fieldnames(T)
                snap[f] = Float64(getfield(base.params, f))
            end
        end

        # 3) Apply IIV (optional)
        if pop.iiv !== nothing
            params_i, d_iiv = _realize_params_log_normal(params_i, omegas, rng)
            realized = copy(snap)
            for (k, v) in d_iiv
                realized[k] = v
            end
            realized_params[i] = realized
        else
            realized_params[i] = snap
        end

        # 4) Simulate with or without IOV, with or without PKPD
        spec_i = ModelSpec(base.kind, base.name * "_i$(i)", params_i, base.doses)

        if pd_spec === nothing
            # PK-only simulation
            if pop.iov === nothing
                individuals[i] = simulate(spec_i, grid, solver)
            else
                # Individual-specific kappas stream, deterministic
                starts = occ_starts::Vector{Float64}
                kappas = sample_iov_kappas(
                    pop.iov.pis, length(starts), pop.iov.seed + UInt64(i)
                )

                pk_params_segments = Vector{typeof(params_i)}(undef, length(starts))
                for occ in 1:length(starts)
                    p_occ, _ = apply_iov(params_i, kappas[occ])
                    pk_params_segments[occ] = p_occ
                end

                individuals[i] = simulate_segmented_pk(
                    spec_i, grid, solver, starts, pk_params_segments
                )
            end
        else
            # Coupled PKPD simulation
            if pop.iov === nothing
                individuals[i] = simulate_pkpd_coupled(spec_i, pd_spec, grid, solver)
            else
                # Individual-specific kappas stream, deterministic
                starts = occ_starts::Vector{Float64}
                kappas = sample_iov_kappas(
                    pop.iov.pis, length(starts), pop.iov.seed + UInt64(i)
                )

                pk_params_segments = Vector{typeof(params_i)}(undef, length(starts))
                for occ in 1:length(starts)
                    p_occ, _ = apply_iov(params_i, kappas[occ])
                    pk_params_segments[occ] = p_occ
                end

                individuals[i] = simulate_segmented_pkpd_coupled(
                    spec_i, pd_spec, grid, solver, starts, pk_params_segments
                )
            end
        end
    end

    summaries = Dict{Symbol,PopulationSummary}()
    if haskey(individuals[1].observations, :conc)
        summaries[:conc] = compute_population_summary(
            individuals, :conc; probs=[0.05, 0.95]
        )
    end

    # Add PD observation summary if pd_spec is provided
    if pd_spec !== nothing
        pd_obs = pd_spec.output_observation
        if haskey(individuals[1].observations, pd_obs)
            summaries[pd_obs] = compute_population_summary(
                individuals, pd_obs; probs=[0.05, 0.95]
            )
        end
    end

    metadata = Dict{String,Any}(
        "n" => n,
        "seed" => seed,
        "iiv_kind" => iiv_kind,
        "omegas" => Dict(String(k) => v for (k, v) in omegas),
        "engine_version" => "0.1.0",
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
        "covariate_model" =>
            pop.covariate_model === nothing ? "none" : pop.covariate_model.name,
        "iov_kind" => pop.iov === nothing ? "none" : string(typeof(pop.iov.kind)),
        "iov_seed" => pop.iov === nothing ? 0 : Int(pop.iov.seed),
        "iov_pis" => if pop.iov === nothing
            Dict{String,Any}()
        else
            Dict(String(k) => v for (k, v) in pop.iov.pis)
        end,
        "occasion_def_mode" =>
            pop.iov === nothing ? "none" : String(pop.iov.occasion_def.mode),
    )

    return PopulationResult(individuals, realized_params, summaries, metadata)
end
