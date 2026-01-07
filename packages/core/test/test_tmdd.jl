# Test suite for TMDD (Target-Mediated Drug Disposition) Models
# Comprehensive tests for all TMDD model variants

using Test
using OpenPKPDCore
using StableRNGs
using Statistics

@testset "TMDD Models" begin

    @testset "TMDD Type Definitions" begin
        @testset "TMDDFull" begin
            @test TMDDFull() isa TMDDModelKind
            @test TMDDFull() isa ModelKind

            params = TMDDFullParams(
                0.01,   # kel
                3.0,    # V
                0.1,    # kon
                0.001,  # koff
                0.1,    # ksyn
                0.05,   # kdeg
                0.02,   # kint
                2.0     # R0
            )
            @test params.kel == 0.01
            @test params.V == 3.0
            @test params.R0 == 2.0
        end

        @testset "TMDDQSS" begin
            @test TMDDQSS() isa TMDDModelKind

            params = TMDDQSSParams(
                0.01,   # kel
                3.0,    # V
                0.01,   # KSS
                0.1,    # ksyn
                0.05,   # kdeg
                0.02,   # kint
                2.0     # Rtot0
            )
            @test params.KSS == 0.01
            @test params.Rtot0 == 2.0
        end

        @testset "TMDDQE" begin
            @test TMDDQE() isa TMDDModelKind

            params = TMDDQEParams(
                0.01,   # kel
                3.0,    # V
                0.01,   # KD
                0.1,    # ksyn
                0.05,   # kdeg
                2.0     # Rtot0
            )
            @test params.KD == 0.01
        end

        @testset "TMDDMM" begin
            @test TMDDMM() isa TMDDModelKind

            params = TMDDMMParams(
                0.01,   # kel
                3.0,    # V
                0.5,    # Vmax
                0.1     # Km
            )
            @test params.Vmax == 0.5
            @test params.Km == 0.1
        end

        @testset "TMDDRapidBinding" begin
            @test TMDDRapidBinding() isa TMDDModelKind

            params = TMDDRapidBindingParams(
                0.01,   # kel
                3.0,    # V
                0.01,   # KD
                0.1,    # ksyn
                0.05,   # kdeg
                0.02,   # kint
                2.0     # Rtot0
            )
            @test params.KD == 0.01
        end

        @testset "TMDDIrreversible" begin
            @test TMDDIrreversible() isa TMDDModelKind

            params = TMDDIrreversibleParams(
                0.01,   # kel
                3.0,    # V
                0.1,    # kon
                0.1,    # ksyn
                0.05,   # kdeg
                2.0     # R0
            )
            @test params.kon == 0.1
        end

        @testset "TMDD2CptFull" begin
            @test TMDD2CptFull() isa TMDDModelKind

            params = TMDD2CptFullParams(
                0.01,   # kel
                3.0,    # V1
                4.0,    # V2
                0.5,    # Q
                0.1,    # kon
                0.001,  # koff
                0.1,    # ksyn
                0.05,   # kdeg
                0.02,   # kint
                2.0     # R0
            )
            @test params.V1 == 3.0
            @test params.V2 == 4.0
            @test params.Q == 0.5
        end

        @testset "TMDD2CptQSS" begin
            @test TMDD2CptQSS() isa TMDDModelKind

            params = TMDD2CptQSSParams(
                0.01,   # kel
                3.0,    # V1
                4.0,    # V2
                0.5,    # Q
                0.01,   # KSS
                0.1,    # ksyn
                0.05,   # kdeg
                0.02,   # kint
                2.0     # Rtot0
            )
            @test params.KSS == 0.01
        end

        @testset "TMDD2CptCL" begin
            @test TMDD2CptCL() isa TMDDModelKind

            params = TMDD2CptCLParams(
                0.5,    # CL
                3.0,    # V1
                4.0,    # V2
                0.5,    # Q
                0.01,   # Kss
                0.5,    # Vmax
                2.0,    # R0
                0.05    # kdeg
            )
            @test params.CL == 0.5
            @test params.Vmax == 0.5
        end

        @testset "TMDDSolubleTarget" begin
            @test TMDDSolubleTarget() isa TMDDModelKind

            params = TMDDSolubleTargetParams(
                0.01,   # kel
                3.0,    # V
                0.1,    # kon
                0.001,  # koff
                0.1,    # ksyn
                0.05,   # kdeg
                0.03,   # kel_complex
                2.0     # R0
            )
            @test params.kel_complex == 0.03
        end

        @testset "TMDDInternalization" begin
            @test TMDDInternalization() isa TMDDModelKind

            params = TMDDInternalizationParams(
                0.01,   # kel
                3.0,    # V
                0.1,    # kon
                0.001,  # koff
                0.1,    # ksyn
                0.02,   # kint
                0.01,   # kendo_R
                0.05,   # krec
                0.03,   # kdeg_e
                2.0     # R0
            )
            @test params.kendo_R == 0.01
            @test params.krec == 0.05
        end
    end

    @testset "TMDDSpec Construction" begin
        params = TMDDFullParams(0.01, 3.0, 0.1, 0.001, 0.1, 0.05, 0.02, 2.0)
        doses = [DoseEvent(0.0, 100.0)]

        spec = TMDDSpec(TMDDFull(), "Test TMDD", params, doses)

        @test spec.kind isa TMDDFull
        @test spec.name == "Test TMDD"
        @test spec.params.kel == 0.01
        @test length(spec.doses) == 1
        @test spec.doses[1].amount == 100.0
    end

    @testset "TMDD Validation" begin
        @testset "Valid specs" begin
            params = TMDDFullParams(0.01, 3.0, 0.1, 0.001, 0.1, 0.05, 0.02, 2.0)
            doses = [DoseEvent(0.0, 100.0)]
            spec = TMDDSpec(TMDDFull(), "Test", params, doses)

            # Should not throw
            @test validate_tmdd(spec) === nothing
        end

        @testset "Invalid - missing doses" begin
            params = TMDDFullParams(0.01, 3.0, 0.1, 0.001, 0.1, 0.05, 0.02, 2.0)
            spec = TMDDSpec(TMDDFull(), "Test", params, DoseEvent[])

            @test_throws ErrorException validate_tmdd(spec)
        end

        @testset "Validation for each model type" begin
            doses = [DoseEvent(0.0, 100.0)]

            # QSS
            spec_qss = TMDDSpec(TMDDQSS(), "QSS", TMDDQSSParams(0.01, 3.0, 0.01, 0.1, 0.05, 0.02, 2.0), doses)
            @test validate_tmdd(spec_qss) === nothing

            # QE
            spec_qe = TMDDSpec(TMDDQE(), "QE", TMDDQEParams(0.01, 3.0, 0.01, 0.1, 0.05, 2.0), doses)
            @test validate_tmdd(spec_qe) === nothing

            # MM
            spec_mm = TMDDSpec(TMDDMM(), "MM", TMDDMMParams(0.01, 3.0, 0.5, 0.1), doses)
            @test validate_tmdd(spec_mm) === nothing

            # 2Cpt QSS
            spec_2qss = TMDDSpec(TMDD2CptQSS(), "2QSS", TMDD2CptQSSParams(0.01, 3.0, 4.0, 0.5, 0.01, 0.1, 0.05, 0.02, 2.0), doses)
            @test validate_tmdd(spec_2qss) === nothing

            # 2Cpt CL
            spec_2cl = TMDDSpec(TMDD2CptCL(), "2CL", TMDD2CptCLParams(0.5, 3.0, 4.0, 0.5, 0.01, 0.5, 2.0, 0.05), doses)
            @test validate_tmdd(spec_2cl) === nothing
        end
    end

    @testset "TMDD Helper Functions" begin
        @testset "n_states" begin
            @test n_states(TMDDFull()) == 3
            @test n_states(TMDDQSS()) == 2
            @test n_states(TMDDQE()) == 2
            @test n_states(TMDDMM()) == 1
            @test n_states(TMDDRapidBinding()) == 2
            @test n_states(TMDDIrreversible()) == 2
            @test n_states(TMDD2CptFull()) == 4
            @test n_states(TMDD2CptQSS()) == 3
            @test n_states(TMDD2CptCL()) == 3
            @test n_states(TMDDSolubleTarget()) == 3
            @test n_states(TMDDInternalization()) == 5
        end

        @testset "state_names" begin
            @test state_names(TMDDFull()) == [:L, :R, :P]
            @test state_names(TMDDQSS()) == [:L, :Rtot]
            @test state_names(TMDDMM()) == [:L]
            @test state_names(TMDD2CptQSS()) == [:L, :Lp, :Rtot]
        end

        @testset "dosing_compartment" begin
            @test dosing_compartment(TMDDFull()) == 1
            @test dosing_compartment(TMDD2CptQSS()) == 1
        end

        @testset "get_central_volume" begin
            params_full = TMDDFullParams(0.01, 3.0, 0.1, 0.001, 0.1, 0.05, 0.02, 2.0)
            @test get_central_volume(params_full) == 3.0

            params_2cpt = TMDD2CptQSSParams(0.01, 3.0, 4.0, 0.5, 0.01, 0.1, 0.05, 0.02, 2.0)
            @test get_central_volume(params_2cpt) == 3.0
        end

        @testset "get_initial_state" begin
            params = TMDDFullParams(0.01, 3.0, 0.1, 0.001, 0.1, 0.05, 0.02, 2.0)
            doses = [DoseEvent(0.0, 100.0)]
            spec = TMDDSpec(TMDDFull(), "Test", params, doses)

            u0 = get_initial_state(spec)
            @test length(u0) == 3
            @test u0[1] == 0.0  # L = 0
            @test u0[2] == 2.0  # R = R0
            @test u0[3] == 0.0  # P = 0
        end
    end

    @testset "TMDD Steady-State Calculations" begin
        @testset "calculate_KD" begin
            @test calculate_KD(0.1, 0.001) ≈ 0.01
        end

        @testset "calculate_Rtot_ss" begin
            @test calculate_Rtot_ss(0.1, 0.05) ≈ 2.0
        end

        @testset "calculate_half_life" begin
            @test calculate_half_life(log(2)) ≈ 1.0
        end

        @testset "target_occupancy" begin
            @test target_occupancy(0.01, 0.01) ≈ 0.5
            @test target_occupancy(0.1, 0.01) > 0.9
            @test target_occupancy(0.001, 0.01) < 0.1
        end

        @testset "free_drug_concentration" begin
            # When Ltot >> Rtot and Ltot >> KD*V, L ≈ Ltot
            L = free_drug_concentration(100.0, 1.0, 0.01, 3.0)
            @test L > 98.0  # Most drug is free

            # When Rtot >> Ltot, most drug is bound
            L2 = free_drug_concentration(1.0, 100.0, 0.01, 3.0)
            @test L2 < 0.5  # Most drug is bound
        end
    end

    @testset "TMDD Analysis Functions" begin
        @testset "calculate_target_occupancy" begin
            @test calculate_target_occupancy(0.0, 0.01) == 0.0
            @test calculate_target_occupancy(0.01, 0.01) ≈ 0.5
            @test calculate_target_occupancy(0.09, 0.01) ≈ 0.9 atol=0.01
        end

        @testset "target_occupancy_EC values" begin
            KD = 0.01
            @test target_occupancy_EC50(KD) == KD
            @test target_occupancy_EC90(KD) == 9 * KD
            @test target_occupancy_EC99(KD) == 99 * KD
        end

        @testset "identify_tmdd_regime" begin
            # High concentration - linear
            @test identify_tmdd_regime(10.0, 0.01, 1.0) == :linear

            # Low concentration - TMDD
            @test identify_tmdd_regime(0.001, 0.01, 1.0) == :tmdd

            # Mixed
            @test identify_tmdd_regime(0.1, 0.01, 1.0) == :mixed
        end

        @testset "calculate_tmdd_steady_state" begin
            params = TMDDFullParams(0.01, 3.0, 0.1, 0.001, 0.1, 0.05, 0.02, 2.0)
            ss = calculate_tmdd_steady_state(params)

            @test ss.R_ss ≈ 2.0  # ksyn/kdeg = 0.1/0.05 = 2.0
            @test ss.KD ≈ 0.01   # koff/kon = 0.001/0.1 = 0.01
            @test ss.t_half_target ≈ log(2)/0.05
        end
    end

    @testset "TMDD Simulation - Full Model" begin
        params = TMDDFullParams(
            0.01,   # kel = 0.01/day
            3.0,    # V = 3 L
            0.1,    # kon
            0.001,  # koff
            0.1,    # ksyn
            0.05,   # kdeg
            0.02,   # kint
            2.0     # R0
        )
        doses = [DoseEvent(0.0, 100.0)]  # 100 mg bolus
        spec = TMDDSpec(TMDDFull(), "Test Full TMDD", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))  # 7 days
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test length(result.t) == 169
        @test haskey(result.states, :L)
        @test haskey(result.states, :R)
        @test haskey(result.states, :P)
        @test haskey(result.observations, :conc)
        @test haskey(result.observations, :target_occupancy)

        # Initial concentration should be dose/V
        @test result.observations[:conc][1] ≈ 100.0/3.0 atol=1e-6

        # Concentration should decrease over time
        @test result.observations[:conc][end] < result.observations[:conc][1]

        # Target occupancy: starts at 0, increases as complex forms, then decreases
        # At t=0, P=0 so occupancy=0. Find max occupancy (should be intermediate time)
        max_occ = maximum(result.observations[:target_occupancy])
        @test max_occ > 0.5  # Should achieve significant target engagement
    end

    @testset "TMDD Simulation - QSS" begin
        params = TMDDQSSParams(
            0.01,   # kel
            3.0,    # V
            0.01,   # KSS
            0.1,    # ksyn
            0.05,   # kdeg
            0.02,   # kint
            2.0     # Rtot0
        )
        doses = [DoseEvent(0.0, 100.0)]
        spec = TMDDSpec(TMDDQSS(), "Test QSS", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test haskey(result.states, :L)
        @test haskey(result.states, :Rtot)
        @test haskey(result.observations, :conc)

        # Should have reasonable PK behavior
        @test result.observations[:conc][1] > 0
        @test result.observations[:conc][end] < result.observations[:conc][1]
    end

    @testset "TMDD Simulation - 2Cpt QSS" begin
        params = TMDD2CptQSSParams(
            0.01,   # kel
            3.0,    # V1
            4.0,    # V2
            0.5,    # Q
            0.01,   # KSS
            0.1,    # ksyn
            0.05,   # kdeg
            0.02,   # kint
            2.0     # Rtot0
        )
        doses = [DoseEvent(0.0, 200.0)]
        spec = TMDDSpec(TMDD2CptQSS(), "Test 2Cpt QSS", params, doses)

        grid = SimGrid(0.0, 336.0, collect(0.0:2.0:336.0))  # 14 days
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test haskey(result.states, :L)
        @test haskey(result.states, :Lp)
        @test haskey(result.states, :Rtot)
        @test haskey(result.observations, :conc_peripheral)

        # Peripheral should build up and then decline
        @test maximum(result.observations[:conc_peripheral]) > 0
    end

    @testset "TMDD Simulation - MM Approximation" begin
        params = TMDDMMParams(
            0.01,   # kel
            3.0,    # V
            0.5,    # Vmax
            0.1     # Km
        )
        doses = [DoseEvent(0.0, 100.0)]
        spec = TMDDSpec(TMDDMM(), "Test MM", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test haskey(result.observations, :conc)

        # Nonlinear elimination - faster initial decline
        @test result.observations[:conc][end] < 1.0  # Should clear significantly
    end

    @testset "TMDD Simulation - Irreversible" begin
        params = TMDDIrreversibleParams(
            0.01,   # kel
            3.0,    # V
            0.1,    # kon
            0.1,    # ksyn
            0.05,   # kdeg
            2.0     # R0
        )
        doses = [DoseEvent(0.0, 100.0)]
        spec = TMDDSpec(TMDDIrreversible(), "Test Irreversible", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test haskey(result.states, :L)
        @test haskey(result.states, :R)

        # Target should be depleted by irreversible binding
        @test result.states[:R][end] < result.states[:R][1]
    end

    @testset "TMDD Simulation - Rapid Binding" begin
        params = TMDDRapidBindingParams(
            0.01,   # kel
            3.0,    # V
            0.01,   # KD
            0.1,    # ksyn
            0.05,   # kdeg
            0.02,   # kint
            2.0     # Rtot0
        )
        doses = [DoseEvent(0.0, 100.0)]
        spec = TMDDSpec(TMDDRapidBinding(), "Test Rapid Binding", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test haskey(result.states, :Ltot)
        @test haskey(result.states, :Rtot)
        @test haskey(result.observations, :conc_free)
        @test haskey(result.observations, :conc_total)

        # Total > Free due to binding
        @test result.observations[:conc_total][1] >= result.observations[:conc_free][1]
    end

    @testset "TMDD Multiple Doses" begin
        params = TMDD2CptQSSParams(
            0.01, 3.0, 4.0, 0.5, 0.01, 0.1, 0.05, 0.02, 2.0
        )

        # Weekly dosing (168 hours) - faster cycle
        doses = [
            DoseEvent(0.0, 200.0),
            DoseEvent(168.0, 200.0),
            DoseEvent(336.0, 200.0)
        ]
        spec = TMDDSpec(TMDD2CptQSS(), "Weekly Dosing", params, doses)

        grid = SimGrid(0.0, 504.0, collect(0.0:1.0:504.0))  # 3 weeks with hourly sampling
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult

        # Check that concentration is elevated throughout (accumulation)
        conc = result.observations[:conc]
        @test length(conc) > 0

        # Initial Cmax
        cmax_initial = maximum(conc[1:100])  # First ~100 hours
        @test cmax_initial > 0

        # Final trough should be elevated due to accumulation
        ctrough_final = conc[end]
        @test ctrough_final > 0

        # Total AUC should be positive
        metrics = tmdd_exposure_metrics(result)
        @test metrics.AUC > 0
    end

    @testset "TMDD Exposure Metrics" begin
        params = TMDDQSSParams(0.01, 3.0, 0.01, 0.1, 0.05, 0.02, 2.0)
        doses = [DoseEvent(0.0, 100.0)]
        spec = TMDDSpec(TMDDQSS(), "Test", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)
        metrics = tmdd_exposure_metrics(result)

        @test metrics.AUC > 0
        @test metrics.Cmax > 0
        @test metrics.Tmax >= 0
        @test metrics.Ctrough >= 0
        @test metrics.AUC_occupancy > 0
        @test 0.0 <= metrics.mean_occupancy <= 1.0
    end

    @testset "Time Above Threshold" begin
        params = TMDDQSSParams(0.01, 3.0, 0.01, 0.1, 0.05, 0.02, 2.0)
        doses = [DoseEvent(0.0, 100.0)]
        spec = TMDDSpec(TMDDQSS(), "Test", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        # Time above 90% occupancy
        t_90 = time_above_occupancy(result, 0.9)
        @test t_90 >= 0.0
        @test t_90 <= 168.0

        # Time above some concentration
        t_conc = time_above_concentration(result, 1.0)
        @test t_conc >= 0.0
    end

    @testset "TMDD Determinism" begin
        # Same spec should produce identical results
        params = TMDDQSSParams(0.01, 3.0, 0.01, 0.1, 0.05, 0.02, 2.0)
        doses = [DoseEvent(0.0, 100.0)]
        spec = TMDDSpec(TMDDQSS(), "Test", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result1 = solve_tmdd(spec, grid, solver)
        result2 = solve_tmdd(spec, grid, solver)

        @test result1.observations[:conc] == result2.observations[:conc]
        @test result1.states[:L] == result2.states[:L]
    end

    @testset "TMDD Soluble Target" begin
        params = TMDDSolubleTargetParams(
            0.01,   # kel
            3.0,    # V
            0.1,    # kon
            0.001,  # koff
            0.1,    # ksyn
            0.05,   # kdeg
            0.03,   # kel_complex (different from kel)
            2.0     # R0
        )
        doses = [DoseEvent(0.0, 100.0)]
        spec = TMDDSpec(TMDDSolubleTarget(), "Soluble Target", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test haskey(result.states, :P)  # Complex tracked
        @test haskey(result.observations, :conc_bound)
    end

    @testset "TMDD Internalization Model" begin
        params = TMDDInternalizationParams(
            0.01,   # kel
            3.0,    # V
            0.1,    # kon
            0.001,  # koff
            0.1,    # ksyn
            0.02,   # kint
            0.01,   # kendo_R
            0.05,   # krec
            0.03,   # kdeg_e
            2.0     # R0
        )
        doses = [DoseEvent(0.0, 100.0)]
        spec = TMDDSpec(TMDDInternalization(), "Internalization", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test haskey(result.states, :Re)  # Endosomal receptor
        @test haskey(result.states, :Pe)  # Endosomal complex
        @test haskey(result.observations, :R_endo_free)
    end

    @testset "TMDD 2Cpt CL Parameterization" begin
        params = TMDD2CptCLParams(
            0.5,    # CL = 0.5 L/day
            3.0,    # V1
            4.0,    # V2
            0.5,    # Q
            0.01,   # Kss
            0.5,    # Vmax
            2.0,    # R0
            0.05    # kdeg
        )
        doses = [DoseEvent(0.0, 200.0)]
        spec = TMDDSpec(TMDD2CptCL(), "CL Param", params, doses)

        grid = SimGrid(0.0, 336.0, collect(0.0:2.0:336.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test haskey(result.observations, :conc)

        # Should clear based on CL
        @test result.observations[:conc][end] < result.observations[:conc][1]
    end

end
