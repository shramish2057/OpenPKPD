export simulate_pkpd_coupled

function _preset_dose_callback(
    doses::Vector{DoseEvent}, t0::Float64, t1::Float64, target_index::Int
)
    dose_times = Float64[]
    dose_amounts = Float64[]

    for d in doses
        if d.time > t0 && d.time <= t1
            push!(dose_times, d.time)
            push!(dose_amounts, d.amount)
        end
    end

    if isempty(dose_times)
        return nothing
    end

    function affect!(integrator)
        idx = findfirst(==(integrator.t), dose_times)
        if idx === nothing
            error("Internal error: dose time not found for t=$(integrator.t)")
        end
        integrator.u[target_index] += dose_amounts[idx]
    end

    return PresetTimeCallback(dose_times, affect!)
end

"""
General coupled PKPD simulation for:
- PK: any supported ModelSpec{K,P}
- PD: IndirectResponseTurnover

Coupled states:
- first: PK states
- last:  R response

Observations:
- :conc computed by pk_conc
- pd_spec.output_observation maps to response state R
"""
function simulate_pkpd_coupled(
    pk_spec::ModelSpec{K,P},
    pd_spec::PDSpec{IndirectResponseTurnover,IndirectResponseTurnoverParams},
    grid::SimGrid,
    solver::SolverSpec,
) where {K<:ModelKind,P}
    pk_validate(pk_spec)
    validate(pd_spec)
    validate(grid)
    validate(solver)

    kind = pk_spec.kind

    # PK parameter tuple is model-specific
    pkp = pk_param_tuple(pk_spec)

    # PD parameters
    Kin = pd_spec.params.Kin
    Kout = pd_spec.params.Kout
    R0 = pd_spec.params.R0
    Imax = pd_spec.params.Imax
    IC50 = pd_spec.params.IC50

    # Build initial state vector: [PK states..., R]
    u0_pk = pk_u0(pk_spec, grid)
    u0 = vcat(u0_pk, [R0])

    # Combined parameter tuple (NamedTuple), still pure data
    p = merge(pkp, (Kin=Kin, Kout=Kout, Imax=Imax, IC50=IC50))

    n_pk = length(u0_pk)
    r_index = n_pk + 1

    function ode!(du, u, p, t)
        # Split views
        u_pk = @view u[1:n_pk]
        du_pk = @view du[1:n_pk]

        # PK dynamics
        pk_ode!(du_pk, u_pk, p, t, kind)

        # PD dynamics
        R = u[r_index]
        C = pk_conc(u_pk, p, kind)
        I = inhibition(C, p.Imax, p.IC50)
        du[r_index] = p.Kin - p.Kout * (1.0 - I) * R

        return nothing
    end

    prob = ODEProblem(ode!, u0, (grid.t0, grid.t1), p)

    target_index = pk_dose_target_index(kind)
    cb = _preset_dose_callback(pk_spec.doses, grid.t0, grid.t1, target_index)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    # Build state outputs
    state_syms = pk_state_symbols(kind)
    states = Dict{Symbol,Vector{Float64}}()

    for (j, sym) in enumerate(state_syms)
        states[sym] = [u[j] for u in sol.u]
    end

    R = [u[r_index] for u in sol.u]
    states[:R] = R

    # Observations
    C = [pk_conc(@view(u[1:n_pk]), p, kind) for u in sol.u]

    outkey = pd_spec.output_observation
    observations = Dict{Symbol,Vector{Float64}}(:conc => C, outkey => R)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.0",
        "model" => "Coupled PKPD",
        "pk_kind" => string(typeof(kind)),
        "pd_kind" => "IndirectResponseTurnover",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "pk_dose_schedule" => [(d.time, d.amount) for d in pk_spec.doses],
        "pd_output_observation" => String(outkey),
        "pd_params" =>
            Dict("Kin" => Kin, "Kout" => Kout, "R0" => R0, "Imax" => Imax, "IC50" => IC50),
        "deterministic_output_grid" => true,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end
