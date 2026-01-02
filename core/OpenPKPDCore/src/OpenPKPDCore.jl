module OpenPKPDCore

using SciMLBase
using DifferentialEquations

export ModelKind, OneCompIVBolus, OneCompIVBolusParams, DoseEvent, ModelSpec
export SolverSpec, SimGrid, SimResult, simulate

# -------------------------
# Model definitions
# -------------------------

abstract type ModelKind end

"""
One-compartment IV bolus PK model.
"""
struct OneCompIVBolus <: ModelKind end

"""
A single dosing event.

Rules:
- time must be >= 0
- amount must be > 0
"""
struct DoseEvent
    time::Float64
    amount::Float64
end

"""
Typed parameters for OneCompIVBolus.
"""
struct OneCompIVBolusParams
    CL::Float64
    V::Float64
end

"""
Generic model specification container.
"""
struct ModelSpec{K<:ModelKind,P}
    kind::K
    name::String
    params::P
    doses::Vector{DoseEvent}
end

"""
Solver specification is pure data.
"""
struct SolverSpec
    alg::Symbol
    reltol::Float64
    abstol::Float64
    maxiters::Int
end

"""
Simulation time grid definition.
"""
struct SimGrid
    t0::Float64
    t1::Float64
    saveat::Vector{Float64}
end

"""
Simulation result.
"""
struct SimResult
    t::Vector{Float64}
    amount::Vector{Float64}
    conc::Vector{Float64}
    metadata::Dict{String,Any}
end

# -------------------------
# Validation helpers
# -------------------------

function _require_positive(name::String, x::Float64)
    if !(x > 0.0)
        error("Expected positive value for $(name), got $(x)")
    end
    return nothing
end

function validate(spec::ModelSpec{OneCompIVBolus,OneCompIVBolusParams})
    CL = spec.params.CL
    V = spec.params.V

    _require_positive("CL", CL)
    _require_positive("V", V)

    if isempty(spec.doses)
        error("At least one DoseEvent is required")
    end

    for (i, d) in enumerate(spec.doses)
        if d.time < 0.0
            error("DoseEvent time must be >= 0 at index $(i), got $(d.time)")
        end
        _require_positive("DoseEvent amount at index $(i)", d.amount)
    end

    if !issorted([d.time for d in spec.doses])
        error("Dose events must be sorted by time ascending")
    end

    return nothing
end

function validate(grid::SimGrid)
    _require_positive("t1", grid.t1)
    if grid.t0 < 0.0
        error("t0 must be >= 0, got $(grid.t0)")
    end
    if !(grid.t1 > grid.t0)
        error("t1 must be > t0")
    end
    if isempty(grid.saveat)
        error("saveat must not be empty")
    end
    if any(t -> t < grid.t0 || t > grid.t1, grid.saveat)
        error("All saveat values must be within [t0, t1]")
    end
    if !issorted(grid.saveat)
        error("saveat must be sorted ascending")
    end
    return nothing
end

function validate(solver::SolverSpec)
    _require_positive("reltol", solver.reltol)
    _require_positive("abstol", solver.abstol)
    if solver.maxiters < 1
        error("maxiters must be >= 1")
    end
    return nothing
end

# -------------------------
# Solver mapping
# -------------------------

function _solver_alg(alg::Symbol)
    if alg == :Tsit5
        return Tsit5()
    elseif alg == :Rosenbrock23
        return Rosenbrock23()
    else
        error("Unsupported solver alg: $(alg). Supported: :Tsit5, :Rosenbrock23")
    end
end

# -------------------------
# ODE definitions
# -------------------------

function _ode_onecomp_ivbolus!(dA, A, p, t)
    dA[1] = -(p.CL / p.V) * A[1]
    return nothing
end

# -------------------------
# Main simulation entrypoint
# -------------------------

function simulate(
    spec::ModelSpec{OneCompIVBolus,OneCompIVBolusParams}, grid::SimGrid, solver::SolverSpec
)
    validate(spec)
    validate(grid)
    validate(solver)

    CL = spec.params.CL
    V = spec.params.V

    # Apply doses at t0
    A0 = 0.0
    for d in spec.doses
        if d.time == grid.t0
            A0 += d.amount
        end
    end

    p = (CL=CL, V=V)

    u0 = [A0]
    tspan = (grid.t0, grid.t1)

    prob = ODEProblem(_ode_onecomp_ivbolus!, u0, tspan, p)

    dose_times = Float64[]
    dose_amounts = Float64[]

    for d in spec.doses
        if d.time > grid.t0 && d.time <= grid.t1
            push!(dose_times, d.time)
            push!(dose_amounts, d.amount)
        end
    end

    cb = nothing
    if !isempty(dose_times)
        function affect!(integrator)
            idx = findfirst(==(integrator.t), dose_times)
            integrator.u[1] += dose_amounts[idx]
        end
        cb = PresetTimeCallback(dose_times, affect!)
    end

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    A = [u[1] for u in sol.u]
    C = [a / V for a in A]

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.0",
        "model" => "OneCompIVBolus",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount) for d in spec.doses],
        "deterministic_output_grid" => true,
    )

    return SimResult(Vector{Float64}(sol.t), A, C, metadata)
end

end
