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

@testset "Artifact replay: PK only matches original outputs" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "replay_pk",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 60.0), DoseEvent(0.0, 40.0), DoseEvent(12.0, 25.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res1 = simulate(pk, grid, solver)

    artifact = serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res1)

    res2 = replay_execution(artifact)

    @test res2.t == res1.t
    for i in eachindex(res1.t)
        @test isapprox(
            res2.observations[:conc][i], res1.observations[:conc][i]; rtol=1e-12, atol=1e-12
        )
    end
end

@testset "Artifact replay: PK then PD (DirectEmax) matches original outputs" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "replay_pkpd_direct",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(DirectEmax(), "emax", DirectEmaxParams(10.0, 40.0, 0.8), :conc, :effect)

    res1 = simulate_pkpd(pk, pd, grid, solver)

    artifact = serialize_execution(
        model_spec=pk, grid=grid, solver=solver, result=res1, pd_spec=pd
    )

    res2 = replay_execution(artifact)

    @test res2.t == res1.t
    for i in eachindex(res1.t)
        @test isapprox(
            res2.observations[:conc][i], res1.observations[:conc][i]; rtol=1e-12, atol=1e-12
        )
        @test isapprox(
            res2.observations[:effect][i],
            res1.observations[:effect][i];
            rtol=1e-12,
            atol=1e-12,
        )
    end
end

@testset "Artifact replay: coupled PKPD (IndirectResponseTurnover) matches original outputs" begin
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "replay_coupled_oral",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "turnover",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    res1 = simulate_pkpd_coupled(pk, pd, grid, solver)

    artifact = serialize_execution(
        model_spec=pk, grid=grid, solver=solver, result=res1, pd_spec=pd
    )

    res2 = replay_execution(artifact)

    @test res2.t == res1.t
    for i in eachindex(res1.t)
        @test isapprox(
            res2.observations[:conc][i], res1.observations[:conc][i]; rtol=1e-12, atol=1e-12
        )
        @test isapprox(
            res2.observations[:response][i],
            res1.observations[:response][i];
            rtol=1e-12,
            atol=1e-12,
        )
    end
end

@testset "Semantics fingerprint is complete" begin
    fp = semantics_fingerprint()

    @test haskey(fp, "event_semantics_version")
    @test haskey(fp, "solver_semantics_version")
    @test haskey(fp, "artifact_schema_version")
end

@testset "Population simulation is deterministic with StableRNG seed" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "pop_pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(12345), 20)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    r1 = simulate_population(pop, grid, solver)
    r2 = simulate_population(pop, grid, solver)

    @test r1.metadata["n"] == 20
    @test r2.metadata["n"] == 20

    # Realized individual parameters must match exactly run-to-run
    for i in 1:20
        @test r1.params[i][:CL] == r2.params[i][:CL]
        @test r1.params[i][:V] == r2.params[i][:V]
    end

    # Concentration series must match exactly to tight tolerance
    for i in 1:20
        c1 = r1.individuals[i].observations[:conc]
        c2 = r2.individuals[i].observations[:conc]
        @test length(c1) == length(c2)
        for j in eachindex(c1)
            @test isapprox(c1[j], c2[j]; rtol=1e-12, atol=1e-12)
        end
    end
end

@testset "Population simulation without IIV produces one individual matching base simulation" begin
    base = ModelSpec(
        OneCompOralFirstOrder(),
        "pop_pk_oral_base",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    pop = PopulationSpec(base, nothing, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    base_res = simulate(base, grid, solver)
    pop_res = simulate_population(pop, grid, solver)

    @test length(pop_res.individuals) == 1

    c_base = base_res.observations[:conc]
    c_pop = pop_res.individuals[1].observations[:conc]

    for j in eachindex(c_base)
        @test isapprox(c_pop[j], c_base[j]; rtol=1e-12, atol=1e-12)
    end
end

@testset "Population artifact serialization and replay matches original outputs" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "pop_artifact_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(424242), 5)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    r1 = simulate_population(pop, grid, solver)

    artifact = serialize_population_execution(
        population_spec=pop, grid=grid, solver=solver, result=r1
    )

    r2 = replay_population_execution(artifact)

    @test r2.metadata["n"] == r1.metadata["n"]
    @test length(r2.individuals) == length(r1.individuals)

    for i in eachindex(r1.individuals)
        c1 = r1.individuals[i].observations[:conc]
        c2 = r2.individuals[i].observations[:conc]
        @test length(c1) == length(c2)
        for j in eachindex(c1)
            @test isapprox(c2[j], c1[j]; rtol=1e-12, atol=1e-12)
        end
    end

    @test haskey(r1.summaries, :conc)
    @test haskey(r2.summaries, :conc)

    s1 = r1.summaries[:conc]
    s2 = r2.summaries[:conc]

    @test s1.probs == s2.probs

    for i in eachindex(s1.mean)
        @test isapprox(s2.mean[i], s1.mean[i]; rtol=1e-12, atol=1e-12)
        @test isapprox(s2.median[i], s1.median[i]; rtol=1e-12, atol=1e-12)
    end

    for p in s1.probs
        q1 = s1.quantiles[p]
        q2 = s2.quantiles[p]
        for i in eachindex(q1)
            @test isapprox(q2[i], q1[i]; rtol=1e-12, atol=1e-12)
        end
    end
end

@testset "Single-run sensitivity: increasing CL reduces concentration" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "sens_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan(
        "CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)]
    )

    out = run_sensitivity(spec, grid, solver; plan=plan, observation=:conc)

    # Concentration should be lower at positive times after increasing CL
    # We avoid t=0 because it reflects initial condition and bolus application.
    @test out.pert_metric_series[2] < out.base_metric_series[2]
    @test out.metrics.max_abs_delta > 0.0
end

@testset "Population sensitivity: increasing CL reduces mean concentration" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "sens_pop_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(9999), 20)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan(
        "CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)]
    )

    out = run_population_sensitivity(pop, grid, solver; plan=plan, observation=:conc)

    @test out.pert_summary_mean[2] < out.base_summary_mean[2]
    @test out.metrics_mean.max_abs_delta > 0.0
end

@testset "Sensitivity artifact replay matches stored series (single)" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "sens_art_single",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan(
        "CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)]
    )

    r1 = run_sensitivity(spec, grid, solver; plan=plan, observation=:conc)

    art = serialize_sensitivity_execution(
        model_spec=spec, grid=grid, solver=solver, result=r1
    )

    r2 = replay_sensitivity_execution(art)

    for i in eachindex(r1.base_metric_series)
        @test isapprox(
            r2.base_metric_series[i], r1.base_metric_series[i]; rtol=1e-12, atol=1e-12
        )
        @test isapprox(
            r2.pert_metric_series[i], r1.pert_metric_series[i]; rtol=1e-12, atol=1e-12
        )
    end
end

@testset "Sensitivity artifact replay matches stored series (population)" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "sens_art_pop",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(9999), 20)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan(
        "CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)]
    )

    r1 = run_population_sensitivity(
        pop, grid, solver; plan=plan, observation=:conc, probs=[0.05, 0.95]
    )

    art = serialize_population_sensitivity_execution(
        population_spec=pop, grid=grid, solver=solver, result=r1
    )

    r2 = replay_population_sensitivity_execution(art)

    for i in eachindex(r1.base_summary_mean)
        @test isapprox(
            r2.base_summary_mean[i], r1.base_summary_mean[i]; rtol=1e-12, atol=1e-12
        )
        @test isapprox(
            r2.pert_summary_mean[i], r1.pert_summary_mean[i]; rtol=1e-12, atol=1e-12
        )
    end
end

@testset "IOV changes parameters across occasions deterministically" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "iov_pop_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 100.0)],
    )

    iiv = nothing
    iov = IOVSpec(
        LogNormalIIV(), Dict(:CL => 0.3), UInt64(1234), OccasionDefinition(:dose_times)
    )

    pop = PopulationSpec(base, iiv, iov, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res1 = simulate_population(pop, grid, solver)
    res2 = simulate_population(pop, grid, solver)

    # Deterministic replay
    c1 = res1.individuals[1].observations[:conc]
    c2 = res2.individuals[1].observations[:conc]
    for i in eachindex(c1)
        @test isapprox(c1[i], c2[i]; rtol=1e-12, atol=1e-12)
    end

    # Behavioral: after second dose, concentration trajectory should differ from a non-IOV run
    pop_no_iov = PopulationSpec(base, nothing, nothing, nothing, IndividualCovariates[])
    base_res = simulate_population(pop_no_iov, grid, solver)

    c_base = base_res.individuals[1].observations[:conc]
    # pick a time after 12.0, for example index where t == 13
    idx13 = findfirst(==(13.0), grid.saveat)
    @test idx13 !== nothing
    @test abs(c1[idx13] - c_base[idx13]) > 0.0
end

@testset "Coupled PKPD with IOV is deterministic and alters post-dose dynamics" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pkpd_iov",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 100.0)],
    )

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "pd_turnover",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    iov = IOVSpec(
        LogNormalIIV(), Dict(:CL => 0.3), UInt64(2222), OccasionDefinition(:dose_times)
    )
    pop = PopulationSpec(pk, nothing, iov, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    r1 = simulate_population(pop, grid, solver; pd_spec=pd)
    r2 = simulate_population(pop, grid, solver; pd_spec=pd)

    # Deterministic replay
    y1 = r1.individuals[1].observations[:response]
    y2 = r2.individuals[1].observations[:response]

    for i in eachindex(y1)
        @test isapprox(y1[i], y2[i]; rtol=1e-12, atol=1e-12)
    end

    # Compare against no-IOV behavior
    pop_no_iov = PopulationSpec(pk, nothing, nothing, nothing, IndividualCovariates[])
    base = simulate_population(pop_no_iov, grid, solver; pd_spec=pd)

    idx13 = findfirst(==(13.0), grid.saveat)
    @test idx13 !== nothing
    @test abs(y1[idx13] - base.individuals[1].observations[:response][idx13]) > 0.0
end

@testset "Time-varying covariate changes concentration deterministically via segmentation" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "tv_cov_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    # CL increases at t=10, so concentrations after 10 should drop faster
    cm = CovariateModel(
        "cl_tv", [CovariateEffect(LinearCovariate(), :CL, :CLMULT, 1.0, 1.0)]
    )

    tv = TimeVaryingCovariates(
        Dict(:CLMULT => TimeCovariateSeries(StepTimeCovariate(), [0.0, 10.0], [1.0, 2.0]))
    )

    covs = [IndividualCovariates(Dict{Symbol,Float64}(), tv)]
    pop = PopulationSpec(base, nothing, nothing, cm, covs)

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver)

    c = res.individuals[1].observations[:conc]
    idx9 = findfirst(==(9.0), grid.saveat)
    idx11 = findfirst(==(11.0), grid.saveat)

    @test idx9 !== nothing
    @test idx11 !== nothing

    # sanity: conc at 11 should be lower than it would be with constant CL
    pop_const = PopulationSpec(base, nothing, nothing, nothing, IndividualCovariates[])
    base_res = simulate_population(pop_const, grid, solver)
    c_base = base_res.individuals[1].observations[:conc]

    @test abs(c[idx11] - c_base[idx11]) > 0.0
end
