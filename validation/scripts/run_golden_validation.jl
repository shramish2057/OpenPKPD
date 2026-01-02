using Pkg
Pkg.activate("core/OpenPKPDCore")
Pkg.instantiate()

using OpenPKPDCore

function _load_json(path::String)
    return read_execution_json(path)
end

function _require_equal(a, b, label::String)
    if a != b
        error("Mismatch for $(label). Expected $(a), got $(b)")
    end
end

function _compare_vectors(a::Vector{Float64}, b::Vector{Float64}, label::String; rtol=1e-12, atol=1e-12)
    if length(a) != length(b)
        error("Length mismatch for $(label). Expected $(length(a)), got $(length(b))")
    end
    for i in eachindex(a)
        if !isapprox(a[i], b[i]; rtol=rtol, atol=atol)
            error("Value mismatch for $(label) at index $(i). Expected $(a[i]), got $(b[i])")
        end
    end
end

function _compare_dict_of_series(expected::Dict{String, Any}, actual::Dict{Symbol, Vector{Float64}}, label::String)
    exp_keys = sort(collect(keys(expected)))
    act_keys = sort([String(k) for k in keys(actual)])

    _require_equal(exp_keys, act_keys, "$(label) keys")

    for k in exp_keys
        exp_v = [Float64(x) for x in expected[k]]
        act_v = actual[Symbol(k)]
        _compare_vectors(exp_v, act_v, "$(label).$(k)")
    end
end

function validate_one(path::String)
    artifact = _load_json(path)

    schema = String(artifact["artifact_schema_version"])
    if schema != ARTIFACT_SCHEMA_VERSION
        error("Artifact schema version mismatch in $(path). Expected $(ARTIFACT_SCHEMA_VERSION), got $(schema)")
    end

    stored = artifact["result"]
    stored_t = [Float64(x) for x in stored["t"]]
    stored_states = Dict{String, Any}(stored["states"])
    stored_obs = Dict{String, Any}(stored["observations"])
    stored_meta = Dict{String, Any}(stored["metadata"])

    replay = replay_execution(artifact)

    if replay.t != stored_t
        error("Time grid mismatch in $(path). Stored t does not match replay t")
    end

    _compare_dict_of_series(stored_states, replay.states, "states")
    _compare_dict_of_series(stored_obs, replay.observations, "observations")

    if !haskey(stored_meta, "event_semantics_version")
        error("Stored artifact missing event_semantics_version in $(path)")
    end
    if !haskey(stored_meta, "solver_semantics_version")
        error("Stored artifact missing solver_semantics_version in $(path)")
    end

    _require_equal(String(stored_meta["event_semantics_version"]), EVENT_SEMANTICS_VERSION, "event_semantics_version")
    _require_equal(String(stored_meta["solver_semantics_version"]), SOLVER_SEMANTICS_VERSION, "solver_semantics_version")

    return true
end

function main()
    golden_dir = "validation/golden"
    files = filter(f -> endswith(f, ".json"), readdir(golden_dir; join=true))

    if isempty(files)
        error("No golden artifacts found in $(golden_dir)")
    end

    println("Running golden validation on $(length(files)) artifact(s)")

    for f in sort(files)
        println("Validating: $(f)")
        validate_one(f)
    end

    println("Golden validation passed")
end

main()
