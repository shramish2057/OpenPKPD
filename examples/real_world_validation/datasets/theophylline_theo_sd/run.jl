using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using OpenPKPDCore
using JSON

# -------------------------
# CSV Parsing (208.2)
# -------------------------

function parse_maybe_float(s::String)
    ss = strip(s)
    if isempty(ss) || lowercase(ss) in ["na", "nan", "."]
        return nothing
    end
    return parse(Float64, ss)
end

function read_theo_csv(path::String)
    lines = readlines(path)
    isempty(lines) && error("Empty CSV: " * path)

    header = [replace(strip(h), "\"" => "") for h in split(strip(lines[1]), ",")]
    required = ["ID", "TIME", "DV", "AMT", "EVID", "CMT", "WT"]
    header != required && error("CSV header mismatch. Expected: $(required) Got: $(header)")

    rows = []
    for (i, ln) in enumerate(lines[2:end])
        s = strip(ln)
        isempty(s) && continue
        fields = split(s, ",")
        length(fields) != length(required) && error("Bad column count at line $(i+1)")

        # Remove quotes if present
        fields = map(x -> replace(x, "\"" => ""), fields)

        push!(rows, (
            ID = parse(Int, fields[1]),
            TIME = parse(Float64, fields[2]),
            DV = parse_maybe_float(fields[3]),  # 208.2: Handle missing DV
            AMT = parse(Float64, fields[4]),
            EVID = parse(Int, fields[5]),
            CMT = parse(Int, fields[6]),
            WT = parse(Float64, fields[7]),
        ))
    end
    return rows
end

function group_by_id(rows)
    d = Dict{Int,Vector{Any}}()
    for r in rows
        if !haskey(d, r.ID)
            d[r.ID] = Any[]
        end
        push!(d[r.ID], r)
    end
    for (_, v) in d
        sort!(v, by = x -> x.TIME)
    end
    return d
end

# -------------------------
# Model parameters (208.5)
# -------------------------

const KA = 1.59      # 1/hr absorption rate
const CL_REF = 2.75  # L/hr reference clearance
const V_REF = 31.8   # L reference volume
const WT_REF = 70.0  # kg reference weight

function params_fixed()
    return OneCompOralFirstOrderParams(KA, CL_REF, V_REF)
end

function params_wt_scaled(wt::Float64)
    cl = CL_REF * (wt / WT_REF)^0.75
    v = V_REF * (wt / WT_REF)^1.0
    return OneCompOralFirstOrderParams(KA, cl, v)
end

# -------------------------
# Simulation helpers
# -------------------------

function simulate_one(rows_for_id, params, tag::String, dose_mg::Float64, dose_time::Float64)
    obs_rows = filter(r -> r.EVID == 0, rows_for_id)
    obs_times = unique([r.TIME for r in obs_rows])
    sort!(obs_times)

    id = rows_for_id[1].ID
    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "theo_sd_$(tag)_subj_$(id)",
        params,
        [DoseEvent(dose_time, dose_mg)],
    )

    grid = SimGrid(minimum(obs_times), maximum(obs_times), obs_times)
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)
    res = simulate(spec, grid, solver)

    return spec, grid, solver, res
end

function compute_rmse(res, obs_rows)
    pred = res.observations[:conc]
    t_obs = [r.TIME for r in obs_rows]
    dv = [r.DV for r in obs_rows]

    # Map predictions to observation times
    pred_at = Dict{Float64,Float64}()
    for (i, tt) in enumerate(res.t)
        pred_at[tt] = pred[i]
    end
    pred_full = [pred_at[tt] for tt in t_obs]

    # 208.4: RMSE only over non-missing DV
    sse = 0.0
    n = 0
    for i in eachindex(dv)
        if dv[i] === nothing
            continue
        end
        e = pred_full[i] - Float64(dv[i])
        sse += e * e
        n += 1
    end

    return n == 0 ? NaN : sqrt(sse / n)
end

# -------------------------
# Main simulation function
# -------------------------

function simulate_subject(rows_for_id)
    # Find dose event
    dose_rows = filter(r -> r.AMT > 0.0, rows_for_id)
    length(dose_rows) == 0 && error("No dose row found")
    length(dose_rows) > 1 && error("Multiple dose rows; this study expects single-dose")

    wt = dose_rows[1].WT
    amt = dose_rows[1].AMT
    dose_time = dose_rows[1].TIME

    # 208.3: Explicit AMT unit resolution
    # If AMT is small (typical mg/kg scale), convert to mg via WT
    dose_mg = amt < 50.0 ? amt * wt : amt
    dose_rule = amt < 50.0 ? "AMT_as_mg_per_kg_times_WT" : "AMT_as_mg"

    obs_rows = filter(r -> r.EVID == 0, rows_for_id)

    # 208.5: Run both fixed and WT-scaled scenarios
    spec_fixed, grid_fixed, solver_fixed, res_fixed = simulate_one(
        rows_for_id, params_fixed(), "fixed", dose_mg, dose_time
    )
    spec_wt, grid_wt, solver_wt, res_wt = simulate_one(
        rows_for_id, params_wt_scaled(wt), "wt_scaled", dose_mg, dose_time
    )

    # 208.3: Store metadata in artifacts
    res_fixed.metadata["dose_unit_rule"] = dose_rule
    res_fixed.metadata["raw_amt"] = amt
    res_fixed.metadata["wt"] = wt
    res_fixed.metadata["scenario"] = "fixed"

    res_wt.metadata["dose_unit_rule"] = dose_rule
    res_wt.metadata["raw_amt"] = amt
    res_wt.metadata["wt"] = wt
    res_wt.metadata["scenario"] = "wt_scaled"

    # Compute RMSE for both scenarios
    rmse_fixed = compute_rmse(res_fixed, obs_rows)
    rmse_wt_scaled = compute_rmse(res_wt, obs_rows)

    metrics = Dict(
        "id" => rows_for_id[1].ID,
        "wt" => wt,
        "dose_mg" => dose_mg,
        "dose_unit_rule" => dose_rule,
        "rmse_fixed" => rmse_fixed,
        "rmse_wt_scaled" => rmse_wt_scaled,
    )

    return (
        fixed = (spec = spec_fixed, grid = grid_fixed, solver = solver_fixed, res = res_fixed),
        wt_scaled = (spec = spec_wt, grid = grid_wt, solver = solver_wt, res = res_wt),
        metrics = metrics,
    )
end

function main()
    base = "docs/examples/real_world_validation"
    out_dir = joinpath(base, "studies/theophylline_theo_sd/output")
    mkpath(out_dir)

    data_path = joinpath(base, "datasets/theophylline_theo_sd/theo_sd.csv")
    rows = read_theo_csv(data_path)
    byid = group_by_id(rows)

    all_metrics = Any[]

    for id in sort(collect(keys(byid)))
        result = simulate_subject(byid[id])
        push!(all_metrics, result.metrics)

        # Write fixed scenario artifact
        out_fixed = joinpath(out_dir, "subj_$(id)_fixed.json")
        write_execution_json(
            out_fixed;
            model_spec = result.fixed.spec,
            grid = result.fixed.grid,
            solver = result.fixed.solver,
            result = result.fixed.res,
        )
        println("Wrote: $(out_fixed)")

        # Write WT-scaled scenario artifact
        out_wt = joinpath(out_dir, "subj_$(id)_wt_scaled.json")
        write_execution_json(
            out_wt;
            model_spec = result.wt_scaled.spec,
            grid = result.wt_scaled.grid,
            solver = result.wt_scaled.solver,
            result = result.wt_scaled.res,
        )
        println("Wrote: $(out_wt)")
    end

    metrics_path = joinpath(out_dir, "metrics.json")
    open(metrics_path, "w") do io
        JSON.print(io, all_metrics, 2)
    end
    println("Wrote metrics: $(metrics_path)")
end

main()
