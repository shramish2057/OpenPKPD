export serialize_population_execution, write_population_json

function _serialize_iiv(iiv::Union{Nothing,IIVSpec})
    if iiv === nothing
        return nothing
    end
    return Dict(
        "kind" => string(typeof(iiv.kind)),
        "omegas" => Dict(String(k) => v for (k, v) in iiv.omegas),
        "seed" => Int(iiv.seed),
        "n" => iiv.n,
    )
end

function _serialize_population_spec(pop::PopulationSpec)
    return Dict(
        "base_model_spec" => _serialize_model_spec(pop.base_model_spec),
        "iiv" => _serialize_iiv(pop.iiv),
        "covariates" =>
            [Dict(String(k) => v for (k, v) in c.values) for c in pop.covariates],
    )
end

function _serialize_population_result(popres::PopulationResult)
    return Dict(
        "metadata" => popres.metadata,
        "params" => [Dict(String(k) => v for (k, v) in d) for d in popres.params],
        "summaries" => Dict(
            String(k) => _serialize_population_summary(v) for (k, v) in popres.summaries
        ),
        "individuals" => [_serialize_results(r) for r in popres.individuals],
    )
end

"""
Create a population execution artifact as a Dict.

Stores:
- population_spec
- grid
- solver
- population_result (including per-individual outputs)
- schema version and semantics fingerprint
"""
function serialize_population_execution(;
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec,
    result::PopulationResult,
)
    artifact = Dict(
        "artifact_type" => "population",
        "artifact_schema_version" => ARTIFACT_SCHEMA_VERSION,
        "semantics_fingerprint" => semantics_fingerprint(),
        "population_spec" => _serialize_population_spec(population_spec),
        "grid" => _serialize_grid(grid),
        "solver" => _serialize_solver(solver),
        "population_result" => _serialize_population_result(result),
    )
    return artifact
end

function write_population_json(path::AbstractString; kwargs...)
    artifact = serialize_population_execution(; kwargs...)
    open(path, "w") do io
        JSON.print(io, artifact; indent=2)
    end
    return path
end

function _serialize_population_summary(s::PopulationSummary)
    return Dict(
        "observation" => String(s.observation),
        "probs" => s.probs,
        "mean" => s.mean,
        "median" => s.median,
        "quantiles" => Dict(string(p) => v for (p, v) in s.quantiles),
    )
end
