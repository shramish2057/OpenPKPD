using Pkg
Pkg.activate("core/OpenPKPDCore")
Pkg.instantiate()

using OpenPKPDCore
using JSON

function write(path, artifact)
    open(path, "w") do io
        JSON.print(io, artifact)
    end
end

function gen_pk_iv_bolus()
    pk = ModelSpec(
        OneCompIVBolus(),
        "golden_pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [
            DoseEvent(0.0, 60.0),
            DoseEvent(0.0, 40.0),   # duplicate time to lock semantics summing
            DoseEvent(12.0, 25.0),
        ],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_oral()
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "golden_pk_oral",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_then_pd_direct_emax()
    pk = ModelSpec(
        OneCompIVBolus(),
        "golden_pk_then_pd",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    pd = PDSpec(
        DirectEmax(),
        "golden_emax",
        DirectEmaxParams(10.0, 40.0, 0.8),
        :conc,
        :effect,
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

function gen_coupled_pkpd_turnover_oral()
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "golden_coupled_oral",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "golden_turnover",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

function main()
    mkpath("validation/golden")

    artifacts = Dict(
        "pk_iv_bolus.json" => gen_pk_iv_bolus(),
        "pk_oral.json" => gen_pk_oral(),
        "pk_then_pd_direct_emax.json" => gen_pk_then_pd_direct_emax(),
        "pkpd_coupled_turnover_oral.json" => gen_coupled_pkpd_turnover_oral(),
    )

    for (fname, art) in artifacts
        path = joinpath("validation/golden", fname)
        write(path, art)
        println("Wrote: $(path)")
    end
end

main()
