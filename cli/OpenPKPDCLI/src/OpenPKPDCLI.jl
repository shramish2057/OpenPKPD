module OpenPKPDCLI

using ArgParse
using JSON
using OpenPKPDCore

function _die(msg::String)
    println(stderr, msg)
    exit(1)
end

function cmd_version()
    println("OpenPKPD " * OpenPKPDCore.OPENPKPD_VERSION)
    println("Event semantics: " * OpenPKPDCore.EVENT_SEMANTICS_VERSION)
    println("Solver semantics: " * OpenPKPDCore.SOLVER_SEMANTICS_VERSION)
    println("Artifact schema: " * OpenPKPDCore.ARTIFACT_SCHEMA_VERSION)
end

function cmd_replay(args)
    path = args["artifact"]
    out = get(args, "out", nothing)

    artifact = OpenPKPDCore.read_execution_json(path)

    atype = get(artifact, "artifact_type", "single")

    if atype == "population"
        res = OpenPKPDCore.replay_population_execution(artifact)
        if out !== nothing
            parsed = OpenPKPDCore.deserialize_population_execution(artifact)
            OpenPKPDCore.write_population_json(
                out;
                population_spec = parsed.population_spec,
                grid = parsed.grid,
                solver = parsed.solver,
                result = res,
            )
        end
        return
    end

    if atype == "sensitivity_single"
        res = OpenPKPDCore.replay_sensitivity_execution(artifact)
        if out !== nothing
            parsed = OpenPKPDCore.deserialize_sensitivity_execution(artifact)
            OpenPKPDCore.write_sensitivity_json(
                out;
                model_spec = parsed.model_spec,
                grid = parsed.grid,
                solver = parsed.solver,
                result = res,
            )
        end
        return
    end

    if atype == "sensitivity_population"
        res = OpenPKPDCore.replay_population_sensitivity_execution(artifact)
        if out !== nothing
            parsed = OpenPKPDCore.deserialize_population_sensitivity_execution(artifact)
            OpenPKPDCore.write_population_sensitivity_json(
                out;
                population_spec = parsed.population_spec,
                grid = parsed.grid,
                solver = parsed.solver,
                result = res,
            )
        end
        return
    end

    # default: single execution
    res = OpenPKPDCore.replay_execution(artifact)
    if out !== nothing
        parsed = OpenPKPDCore.deserialize_execution(artifact)
        OpenPKPDCore.write_execution_json(
            out;
            model_spec = parsed.model_spec,
            grid = parsed.grid,
            solver = parsed.solver,
            result = res,
            pd_spec = parsed.pd_spec,
        )
    end
end

function cmd_validate_golden()
    cli_src = @__DIR__
    repo_root = normpath(joinpath(cli_src, "..", "..", ".."))

    runner = joinpath(repo_root, "validation", "scripts", "run_golden_validation.jl")

    if !isfile(runner)
        _die("Golden validation runner not found: " * runner)
    end

    # Run validation script in a subprocess to avoid Pkg dependency issues
    cmd = `julia $runner`
    run(Cmd(cmd; dir=repo_root))
end

function main()
    s = ArgParseSettings()
    @add_arg_table s begin
        "command"
            required = true
    end

    length(ARGS) ≥ 1 || _die("Missing command")

    cmd = ARGS[1]
    rest = ARGS[2:end]

    if cmd == "version"
        cmd_version()
    elseif cmd == "replay"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--artifact"
                required = true
            "--out"
                required = false
        end
        args = parse_args(rest, rs)
        cmd_replay(args)
    elseif cmd == "validate-golden"
        cmd_validate_golden()
    else
        _die("Unknown command: $cmd")
    end
end

# Run when invoked as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end # module
