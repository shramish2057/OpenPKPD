module OpenPKPDCore

export ModelSpec, SolverSpec, simulate

"""
Model specification is pure data: states, parameters, equations, and observation mapping.
This is a placeholder that will be expanded with strict validation rules.
"""
struct ModelSpec
    name::String
end

"""
Solver specification is pure data: method choice, tolerances, and event semantics.
No implicit defaults in the engine.
"""
struct SolverSpec
    name::String
    reltol::Float64
    abstol::Float64
end

"""
simulate binds a ModelSpec to a SolverSpec and executes a deterministic run.
This is a stub.
"""
function simulate(model::ModelSpec, solver::SolverSpec)
    return (
        model=model,
        solver=solver,
        outputs=Dict{String,Any}(),
        metadata=Dict("engine_version" => "0.1.0", "deterministic" => true),
    )
end

end
