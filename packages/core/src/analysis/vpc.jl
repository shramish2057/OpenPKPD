# Visual Predictive Check (VPC) Computation
# Implements standard VPC and prediction-corrected VPC (pcVPC)

using StableRNGs

export compute_vpc, compute_pcvpc

"""
Compute Visual Predictive Check from observed data and population simulations.

This function computes a VPC by:
1. Binning observed data by time
2. Computing observed percentiles in each bin
3. Running population simulations
4. Computing simulated percentiles with bootstrap CIs

Arguments:
- observed: ObservedData from CDISC or other source
- population_spec: PopulationSpec for simulations
- grid: SimGrid defining simulation time points
- solver: SolverSpec for ODE solving
- config: VPCConfig with analysis settings
- error_spec: Optional ResidualErrorSpec for adding residual error

Returns:
- VPCResult
"""
function compute_vpc(
    observed::ObservedData,
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec;
    config::VPCConfig=VPCConfig(),
    error_spec::Union{Nothing,ResidualErrorSpec}=nothing
)::VPCResult
    rng = StableRNG(config.seed)

    # Extract observed times and values
    obs_times = all_times(observed)
    obs_values = all_observations(observed)

    # Compute bin definitions based on observed data
    bin_defs = compute_bins(obs_times, config.binning)

    if isempty(bin_defs)
        return VPCResult(
            config, VPCBin[], n_subjects(observed), n_observations(observed),
            0, "", config.seed
        )
    end

    # Assign observed data to bins and compute observed percentiles
    obs_binned = assign_to_bins(obs_times, obs_values, bin_defs)

    # Run population simulations
    sim_percentiles_per_bin = _run_simulations_for_vpc(
        population_spec, grid, solver, bin_defs, config, error_spec, rng
    )

    # Build VPCBin results
    vpc_bins = VPCBin[]

    for (bin_def, (bin_id, obs_vals)) in zip(bin_defs, obs_binned)
        percentile_data = VPCPercentileData[]

        for pi_level in config.pi_levels
            # Observed percentile
            obs_pctl = compute_percentile(obs_vals, pi_level)

            # Simulated percentiles (from all simulations)
            sim_pctls = sim_percentiles_per_bin[bin_id][pi_level]

            # Bootstrap CI for simulated percentiles
            lower, median_sim, upper = bootstrap_percentile_ci(
                sim_pctls, 0.5, config.ci_level, config.n_bootstrap, rng
            )

            push!(percentile_data, VPCPercentileData(
                pi_level, obs_pctl, median_sim, lower, upper
            ))
        end

        push!(vpc_bins, VPCBin(
            bin_def.id,
            bin_def.lower,
            bin_def.upper,
            bin_def.midpoint,
            length(obs_vals),
            config.n_simulations,
            percentile_data
        ))
    end

    return VPCResult(
        config,
        vpc_bins,
        n_subjects(observed),
        n_observations(observed),
        config.n_simulations,
        "",
        config.seed
    )
end

"""
Compute prediction-corrected VPC (pcVPC).

In pcVPC, both observed and simulated data are normalized by the
population prediction to reduce variability from the structural model.

Arguments:
- observed: ObservedData from CDISC or other source
- population_spec: PopulationSpec for simulations
- grid: SimGrid defining simulation time points
- solver: SolverSpec for ODE solving
- config: VPCConfig with analysis settings
- error_spec: Optional ResidualErrorSpec for adding residual error

Returns:
- VPCResult with prediction-corrected values
"""
function compute_pcvpc(
    observed::ObservedData,
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec;
    config::VPCConfig=VPCConfig(),
    error_spec::Union{Nothing,ResidualErrorSpec}=nothing
)::VPCResult
    # Enable prediction correction in config
    pc_config = VPCConfig(
        pi_levels=config.pi_levels,
        ci_level=config.ci_level,
        binning=config.binning,
        prediction_corrected=true,
        stratify_by=config.stratify_by,
        lloq=config.lloq,
        n_simulations=config.n_simulations,
        n_bootstrap=config.n_bootstrap,
        seed=config.seed
    )

    return compute_vpc(observed, population_spec, grid, solver;
        config=pc_config, error_spec=error_spec)
end

"""
Run simulations and compute percentiles for each bin.
"""
function _run_simulations_for_vpc(
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec,
    bin_defs::Vector{BinDefinition},
    config::VPCConfig,
    error_spec::Union{Nothing,ResidualErrorSpec},
    rng
)::Dict{Int,Dict{Float64,Vector{Float64}}}
    n_bins = length(bin_defs)

    # Initialize storage: bin_id -> pi_level -> Vector of percentiles from each simulation
    result = Dict{Int,Dict{Float64,Vector{Float64}}}()
    for bin_def in bin_defs
        result[bin_def.id] = Dict{Float64,Vector{Float64}}()
        for pi in config.pi_levels
            result[bin_def.id][pi] = Float64[]
        end
    end

    # Run simulations
    for sim_idx in 1:config.n_simulations
        # Create a new seed for this simulation
        sim_seed = rand(rng, UInt64)

        # Create new IIV spec with the new seed (if IIV is present)
        if population_spec.iiv !== nothing
            new_iiv = IIVSpec(
                population_spec.iiv.kind,
                population_spec.iiv.omegas,
                sim_seed,
                population_spec.iiv.n
            )
            pop_spec_sim = PopulationSpec(
                population_spec.base_model_spec,
                new_iiv,
                population_spec.iov,
                population_spec.covariate_model,
                population_spec.covariates
            )
        else
            pop_spec_sim = population_spec
        end

        # Run population simulation
        pop_result = simulate_population(pop_spec_sim, grid, solver)

        # Extract simulated observations at time points matching bins
        sim_times = Float64[]
        sim_values = Float64[]

        for individual in pop_result.individuals
            append!(sim_times, individual.t)
            append!(sim_values, individual.observations[:conc])
        end

        # Apply residual error if specified
        if error_spec !== nothing
            error_rng = StableRNG(sim_seed + UInt64(1))
            sim_values = apply_residual_error(sim_values, error_spec; rng=error_rng)
        end

        # Assign to bins and compute percentiles
        sim_binned = assign_to_bins(sim_times, sim_values, bin_defs)

        for (bin_id, sim_vals) in sim_binned
            for pi in config.pi_levels
                pctl = compute_percentile(sim_vals, pi)
                if !isnan(pctl)
                    push!(result[bin_id][pi], pctl)
                end
            end
        end
    end

    return result
end

"""
Compute VPC from population simulation result (without separate observed data).

This is useful for simulation-based model evaluation when you have
simulated data and want to compute prediction intervals.

Arguments:
- pop_result: PopulationResult from simulate_population
- config: VPCConfig with analysis settings

Returns:
- VPCResult
"""
function compute_vpc_from_simulation(
    pop_result::PopulationResult;
    config::VPCConfig=VPCConfig()
)::VPCResult
    rng = StableRNG(config.seed)

    # Extract all times and values from simulation
    all_t = Float64[]
    all_v = Float64[]

    for individual in pop_result.individuals
        append!(all_t, individual.t)
        append!(all_v, individual.observations[:conc])
    end

    # Compute bins
    bin_defs = compute_bins(all_t, config.binning)

    if isempty(bin_defs)
        return VPCResult(config, VPCBin[], length(pop_result.individuals), 0, 0, "", config.seed)
    end

    # Assign data to bins
    binned = assign_to_bins(all_t, all_v, bin_defs)

    # Build VPCBin results (simulated only, no observed)
    vpc_bins = VPCBin[]

    for (bin_def, (bin_id, vals)) in zip(bin_defs, binned)
        percentile_data = VPCPercentileData[]

        for pi_level in config.pi_levels
            pctl = compute_percentile(vals, pi_level)

            # Bootstrap CI
            lower, median_val, upper = bootstrap_percentile_ci(
                vals, pi_level, config.ci_level, config.n_bootstrap, rng
            )

            # For simulation-only VPC, observed = simulated
            push!(percentile_data, VPCPercentileData(
                pi_level, pctl, median_val, lower, upper
            ))
        end

        push!(vpc_bins, VPCBin(
            bin_def.id,
            bin_def.lower,
            bin_def.upper,
            bin_def.midpoint,
            length(vals),
            length(vals),
            percentile_data
        ))
    end

    return VPCResult(
        config,
        vpc_bins,
        length(pop_result.individuals),
        length(all_t),
        1,  # Single "simulation" (the actual pop_result)
        "",
        config.seed
    )
end

export compute_vpc_from_simulation
