using Test
using OpenPKPDCore

function analytic_onecomp_ivbolus_conc(
    t::Float64, doses::Vector{DoseEvent}, CL::Float64, V::Float64
)
    k = CL / V
    c = 0.0
    for d in doses
        if t >= d.time
            dt = t - d.time
            c += (d.amount / V) * exp(-k * dt)
        end
    end
    return c
end

@testset "OneCompIVBolus analytic equivalence" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "1c_iv_bolus",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_ivbolus_conc(t, spec.doses, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "Dose event at nonzero time" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "1c_iv_bolus_delayed",
        OneCompIVBolusParams(3.0, 30.0),
        [DoseEvent(2.0, 120.0)],
    )

    grid = SimGrid(0.0, 12.0, collect(0.0:0.25:12.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        if any(d.time == t for d in spec.doses)
            continue  # skip discontinuity (left-continuous solver output)
        end
        c_ref = analytic_onecomp_ivbolus_conc(t, spec.doses, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "Multiple bolus schedule" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "1c_iv_bolus_multi",
        OneCompIVBolusParams(4.0, 40.0),
        [DoseEvent(0.0, 80.0), DoseEvent(6.0, 50.0), DoseEvent(10.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        if any(d.time == t for d in spec.doses)
            continue
        end
        c_ref = analytic_onecomp_ivbolus_conc(t, spec.doses, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

function analytic_onecomp_oral_first_order_conc(
    t::Float64, doses, Ka::Float64, CL::Float64, V::Float64
)
    k = CL / V
    c = 0.0
    for d in doses
        if t >= d.time
            dt = t - d.time
            # Handle Ka close to k for numerical stability
            if abs(Ka - k) < 1e-12
                # Limit as Ka -> k:
                # C(t) = (Dose/V) * Ka * dt * exp(-k*dt)
                c += (d.amount / V) * Ka * dt * exp(-k * dt)
            else
                c += (d.amount / V) * (Ka / (Ka - k)) * (exp(-k * dt) - exp(-Ka * dt))
            end
        end
    end
    return c
end

@testset "OneCompOralFirstOrder analytic equivalence" begin
    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "1c_oral_fo",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    Ka = spec.params.Ka
    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_oral_first_order_conc(t, spec.doses, Ka, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "OneCompOralFirstOrder multiple doses" begin
    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "1c_oral_fo_multi",
        OneCompOralFirstOrderParams(0.9, 4.0, 40.0),
        [DoseEvent(0.0, 80.0), DoseEvent(8.0, 50.0), DoseEvent(16.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    Ka = spec.params.Ka
    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_oral_first_order_conc(t, spec.doses, Ka, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "State outputs are present and aligned" begin
    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "state_presence",
        OneCompOralFirstOrderParams(1.0, 3.0, 30.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 12.0, collect(0.0:1.0:12.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_gut)
    @test haskey(res.states, :A_central)
    @test haskey(res.observations, :conc)

    for v in values(res.states)
        @test length(v) == length(res.t)
    end
end

function direct_emax(C::Float64, E0::Float64, Emax::Float64, EC50::Float64)
    return E0 + (Emax * C) / (EC50 + C)
end

@testset "PKPD coupling with DirectEmax using IV bolus analytic reference" begin
    pk = ModelSpec(
        OneCompIVBolus(), "pk_iv", OneCompIVBolusParams(5.0, 50.0), [DoseEvent(0.0, 100.0)]
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(DirectEmax(), "pd_emax", DirectEmaxParams(10.0, 40.0, 0.8), :conc, :effect)

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :effect)
    @test length(res.observations[:effect]) == length(res.t)

    CL = pk.params.CL
    V = pk.params.V

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_ivbolus_conc(t, pk.doses, CL, V)
        e_ref = direct_emax(c_ref, pd.params.E0, pd.params.Emax, pd.params.EC50)

        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
        @test isapprox(res.observations[:effect][i], e_ref; rtol=1e-10, atol=1e-12)
    end
end

function analytic_turnover_R(t::Float64, Kin::Float64, Kout::Float64, R0::Float64)
    # dR/dt = Kin - Kout*R
    return R0 * exp(-Kout * t) + (Kin / Kout) * (1.0 - exp(-Kout * t))
end

@testset "IndirectResponseTurnover coupled: Imax=0 matches analytic turnover and PK analytic" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv_for_pd",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(
        IndirectResponseTurnover(),
        "turnover_no_effect",
        IndirectResponseTurnoverParams(
            10.0,  # Kin
            0.5,   # Kout
            15.0,  # R0
            0.0,   # Imax, no drug effect
            1.0,   # IC50, irrelevant here
        ),
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :response)

    CL = pk.params.CL
    V = pk.params.V

    Kin = pd.params.Kin
    Kout = pd.params.Kout
    R0 = pd.params.R0

    for (i, t) in enumerate(res.t)
        if any(d.time == t for d in pk.doses)
            continue  # skip discontinuities (left-continuous solver output)
        end

        c_ref = analytic_onecomp_ivbolus_conc(t, pk.doses, CL, V)
        r_ref = analytic_turnover_R(t, Kin, Kout, R0)

        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
        @test isapprox(res.observations[:response][i], r_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "IndirectResponseTurnover coupled: inhibition raises response above baseline" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv_effect",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    # Baseline at steady state to make interpretation clean
    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "turnover_with_effect",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    r = res.observations[:response]
    @test maximum(r) > Rss
end

@testset "Coupled engine generalization: oral PK with Imax=0 matches analytic PK and analytic turnover" begin
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "pk_oral_for_pd",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(
        IndirectResponseTurnover(),
        "turnover_no_effect_oral",
        IndirectResponseTurnoverParams(
            10.0,  # Kin
            0.5,   # Kout
            15.0,  # R0
            0.0,   # Imax
            1.0,   # IC50
        ),
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    Ka = pk.params.Ka
    CL = pk.params.CL
    V = pk.params.V

    Kin = pd.params.Kin
    Kout = pd.params.Kout
    R0 = pd.params.R0

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_oral_first_order_conc(t, pk.doses, Ka, CL, V)
        r_ref = analytic_turnover_R(t, Kin, Kout, R0)

        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
        @test isapprox(res.observations[:response][i], r_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "Event semantics: duplicate times are summed (IV bolus)" begin
    base_params = OneCompIVBolusParams(5.0, 50.0)
    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pk_dup = ModelSpec(
        OneCompIVBolus(),
        "dup_times",
        base_params,
        [
            DoseEvent(0.0, 60.0),
            DoseEvent(0.0, 40.0),
            DoseEvent(12.0, 10.0),
            DoseEvent(12.0, 15.0),
        ],
    )

    pk_sum = ModelSpec(
        OneCompIVBolus(),
        "summed",
        base_params,
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 25.0)],
    )

    r_dup = simulate(pk_dup, grid, solver)
    r_sum = simulate(pk_sum, grid, solver)

    @test r_dup.metadata["event_semantics_version"] == "1.0.0"
    @test r_sum.metadata["event_semantics_version"] == "1.0.0"

    for i in eachindex(r_dup.t)
        @test isapprox(
            r_dup.observations[:conc][i],
            r_sum.observations[:conc][i];
            rtol=1e-12,
            atol=1e-12,
        )
    end
end

@testset "Event semantics: input ordering does not matter for same-time events" begin
    params = OneCompIVBolusParams(5.0, 50.0)
    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    doses_a = [
        DoseEvent(0.0, 30.0),
        DoseEvent(0.0, 70.0),
        DoseEvent(12.0, 10.0),
        DoseEvent(12.0, 15.0),
    ]
    doses_b = [
        DoseEvent(0.0, 70.0),
        DoseEvent(0.0, 30.0),
        DoseEvent(12.0, 15.0),
        DoseEvent(12.0, 10.0),
    ]

    pk_a = ModelSpec(OneCompIVBolus(), "a", params, doses_a)
    pk_b = ModelSpec(OneCompIVBolus(), "b", params, doses_b)

    r_a = simulate(pk_a, grid, solver)
    r_b = simulate(pk_b, grid, solver)

    for i in eachindex(r_a.t)
        @test isapprox(
            r_a.observations[:conc][i], r_b.observations[:conc][i]; rtol=1e-12, atol=1e-12
        )
    end
end

@testset "Event semantics: duplicate times are summed in coupled PKPD" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_dup_coupled",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 60.0), DoseEvent(0.0, 40.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "pd_coupled",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.0, 1.0),
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test res.metadata["event_semantics_version"] == "1.0.0"
end

@testset "Solver semantics version is present and stable" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "solver_semantics_test",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 12.0, collect(0.0:1.0:12.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.metadata, "solver_semantics_version")
    @test res.metadata["solver_semantics_version"] == "1.0.0"
end

@testset "Execution artifact serialization is complete and self-consistent" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "serialize_test",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 12.0, collect(0.0:1.0:12.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    res = simulate(pk, grid, solver)

    artifact = serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)

    @test artifact["artifact_schema_version"] == "1.0.0"
    @test haskey(artifact, "model_spec")
    @test haskey(artifact, "grid")
    @test haskey(artifact, "solver")
    @test haskey(artifact, "result")

    @test artifact["result"]["metadata"]["event_semantics_version"] == "1.0.0"
    @test artifact["result"]["metadata"]["solver_semantics_version"] == "1.0.0"

    @test artifact["model_spec"]["name"] == "serialize_test"
    @test artifact["solver"]["alg"] == "Tsit5"
end
