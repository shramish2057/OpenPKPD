export deserialize_population_execution, replay_population_execution

function _parse_iiv(d)::Union{Nothing,IIVSpec}
    if d === nothing
        return nothing
    end

    kind_str = String(d["kind"])

    # Normalize both qualified and unqualified forms
    if kind_str == "LogNormalIIV" || kind_str == "OpenPKPDCore.LogNormalIIV"
        kind = LogNormalIIV()
    else
        error("Unsupported IIV kind in artifact: $(kind_str)")
    end

    omegas = Dict{Symbol,Float64}()
    for (k, v) in d["omegas"]
        omegas[Symbol(String(k))] = Float64(v)
    end

    seed = UInt64(Int(d["seed"]))
    n = Int(d["n"])

    return IIVSpec(kind, omegas, seed, n)
end

function _parse_covariates(arr)::Vector{IndividualCovariates}
    covs = IndividualCovariates[]
    for item in arr
        d = Dict{Symbol,Float64}()
        for (k, v) in item
            d[Symbol(String(k))] = Float64(v)
        end
        push!(covs, IndividualCovariates(d))
    end
    return covs
end

function _parse_population_spec(d)::PopulationSpec
    base = _parse_model_spec(_to_dict(d["base_model_spec"]))
    iiv = _parse_iiv(_to_dict(d["iiv"]))
    covs = _parse_covariates(d["covariates"])
    return PopulationSpec(base, iiv, covs)
end

"""
Deserialize a population execution artifact into objects.
Returns a NamedTuple:
- population_spec
- grid
- solver
"""
function deserialize_population_execution(artifact::Dict)
    schema = String(artifact["artifact_schema_version"])
    if schema != ARTIFACT_SCHEMA_VERSION
        error(
            "Unsupported artifact schema version: $(schema). Expected: $(ARTIFACT_SCHEMA_VERSION)",
        )
    end

    if haskey(artifact, "artifact_type")
        if String(artifact["artifact_type"]) != "population"
            error("Expected artifact_type=population")
        end
    end

    pop_spec = _parse_population_spec(_to_dict(artifact["population_spec"]))
    grid = _parse_grid(_to_dict(artifact["grid"]))
    solver = _parse_solver(_to_dict(artifact["solver"]))

    return (population_spec=pop_spec, grid=grid, solver=solver)
end

"""
Replay a population artifact by rerunning simulate_population.
Returns PopulationResult.
"""
function replay_population_execution(artifact::Dict)::PopulationResult
    parsed = deserialize_population_execution(artifact)
    return simulate_population(parsed.population_spec, parsed.grid, parsed.solver)
end
