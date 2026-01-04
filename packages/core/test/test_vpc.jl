# Test suite for Visual Predictive Check (VPC)

using Test
using OpenPKPDCore
using StableRNGs

@testset "VPC" begin

    @testset "Binning Strategies" begin
        @testset "QuantileBinning" begin
            strategy = QuantileBinning(5)
            @test strategy.n_bins == 5

            times = collect(1.0:20.0)
            bins = compute_bins(times, strategy)

            @test length(bins) == 5
            @test bins[1].id == 1
            @test bins[5].id == 5

            # All times should be covered
            for t in times
                covered = any(b.lower <= t <= b.upper for b in bins)
                @test covered
            end
        end

        @testset "EqualWidthBinning" begin
            strategy = EqualWidthBinning(4)
            @test strategy.n_bins == 4

            times = [0.0, 1.0, 2.0, 3.0, 4.0, 8.0, 12.0, 24.0]
            bins = compute_bins(times, strategy)

            @test length(bins) == 4
            # Equal width means equal time range per bin
            widths = [b.upper - b.lower for b in bins]
            @test all(isapprox(w, widths[1]; atol=0.01) for w in widths)
        end

        @testset "KMeansBinning" begin
            strategy = KMeansBinning(3; max_iter=50)
            @test strategy.n_bins == 3
            @test strategy.max_iter == 50

            # Times with clear clusters
            times = vcat(fill(1.0, 5), fill(5.0, 5), fill(10.0, 5))
            bins = compute_bins(times, strategy)

            @test length(bins) == 3
        end

        @testset "Empty times" begin
            strategy = QuantileBinning(5)
            bins = compute_bins(Float64[], strategy)
            @test isempty(bins)
        end
    end

    @testset "Percentile Computation" begin
        values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

        @test compute_percentile(values, 0.5) == 5.5  # Median
        @test compute_percentile(values, 0.0) == 1.0
        @test compute_percentile(values, 1.0) == 10.0
        @test isapprox(compute_percentile(values, 0.25), 3.25; atol=0.1)
        @test isapprox(compute_percentile(values, 0.75), 7.75; atol=0.1)

        # Empty vector
        @test isnan(compute_percentile(Float64[], 0.5))

        # Single value
        @test compute_percentile([5.0], 0.5) == 5.0
    end

    @testset "Bootstrap CI" begin
        rng = StableRNG(42)
        values = collect(1.0:100.0)

        lower, median_val, upper = bootstrap_percentile_ci(values, 0.5, 0.95, 500, rng)

        # Median should be around 50
        @test 45 < median_val < 55

        # Lower should be less than median
        @test lower < median_val

        # Upper should be greater than median
        @test upper > median_val

        # CI should contain the true median
        @test lower < 50.5 < upper
    end

    @testset "Assign to Bins" begin
        times = [0.5, 1.5, 2.5, 3.5, 4.5]
        values = [10.0, 20.0, 30.0, 40.0, 50.0]

        bins = [
            BinDefinition(1, 0.0, 2.0, 1.0),
            BinDefinition(2, 2.0, 4.0, 3.0),
            BinDefinition(3, 4.0, 6.0, 5.0)
        ]

        result = assign_to_bins(times, values, bins)

        @test length(result) == 3
        @test result[1][1] == 1  # bin id
        @test result[1][2] == [10.0, 20.0]  # values in bin 1
        @test result[2][2] == [30.0, 40.0]  # values in bin 2
        @test result[3][2] == [50.0]  # values in bin 3
    end

    @testset "VPCConfig" begin
        # Default config
        config1 = VPCConfig()
        @test config1.pi_levels == [0.05, 0.50, 0.95]
        @test config1.ci_level == 0.95
        @test config1.binning isa QuantileBinning
        @test config1.prediction_corrected == false
        @test config1.n_simulations == 200
        @test config1.n_bootstrap == 500

        # Custom config
        config2 = VPCConfig(
            pi_levels=[0.10, 0.50, 0.90],
            ci_level=0.90,
            binning=EqualWidthBinning(8),
            prediction_corrected=true,
            n_simulations=100,
            n_bootstrap=200,
            seed=UInt64(99999)
        )
        @test config2.pi_levels == [0.10, 0.50, 0.90]
        @test config2.ci_level == 0.90
        @test config2.binning isa EqualWidthBinning
        @test config2.prediction_corrected == true
        @test config2.n_simulations == 100
        @test config2.seed == UInt64(99999)
    end

    @testset "VPCResult Accessors" begin
        # Create mock VPC result
        percentiles = [
            VPCPercentileData(0.05, 1.0, 1.1, 0.9, 1.3),
            VPCPercentileData(0.50, 5.0, 5.2, 4.8, 5.6),
            VPCPercentileData(0.95, 10.0, 9.8, 9.5, 10.2)
        ]

        bin1 = VPCBin(1, 0.0, 4.0, 2.0, 10, 100, percentiles)
        bin2 = VPCBin(2, 4.0, 8.0, 6.0, 8, 100, percentiles)

        config = VPCConfig()
        result = VPCResult(config, [bin1, bin2], 10, 18, 100, "", UInt64(12345))

        # Test accessors
        @test bin_midpoints(result) == [2.0, 6.0]
        @test observed_percentile(result, 0.50) == [5.0, 5.0]
        @test simulated_median(result, 0.50) == [5.2, 5.2]
        @test simulated_lower(result, 0.50) == [4.8, 4.8]
        @test simulated_upper(result, 0.50) == [5.6, 5.6]
    end

    @testset "VPC from Population Simulation" begin
        # Create a simple one-compartment model
        model_params = OneCompIVBolusParams(10.0, 50.0)  # CL, V
        doses = [DoseEvent(0.0, 100.0)]
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

        # Create population spec
        iiv_spec = IIVSpec(LogNormalIIV(), Dict(:CL => 0.3, :V => 0.2), UInt64(42), 20)

        pop_spec = PopulationSpec(model_spec, iiv_spec, nothing, nothing, IndividualCovariates[])

        grid = SimGrid(0.0, 24.0, collect(0.0:2.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        # Run population simulation
        pop_result = simulate_population(pop_spec, grid, solver)

        # Compute VPC from simulation
        config = VPCConfig(
            pi_levels=[0.10, 0.50, 0.90],
            binning=QuantileBinning(6),
            n_bootstrap=100,
            seed=UInt64(123)
        )

        vpc_result = compute_vpc_from_simulation(pop_result; config=config)

        @test vpc_result.n_subjects_observed == 20
        @test length(vpc_result.bins) == 6
        @test vpc_result.n_simulations == 1

        # Each bin should have 3 percentiles
        for bin in vpc_result.bins
            @test length(bin.percentiles) == 3
        end

        # Midpoints should be sorted
        midpoints = bin_midpoints(vpc_result)
        @test issorted(midpoints)

        # Median percentile should be reasonable
        median_pctls = observed_percentile(vpc_result, 0.50)
        @test all(!isnan(p) for p in median_pctls)
    end

    @testset "VPC with Observed Data" begin
        # Create observed data
        doses = [DoseEvent(0.0, 100.0)]

        # Create subjects with observations
        subj1 = SubjectData(
            "SUBJ001",
            [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
            [2.0, 1.8, 1.5, 1.2, 0.8, 0.5, 0.2],
            doses
        )
        subj2 = SubjectData(
            "SUBJ002",
            [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
            [2.2, 1.9, 1.6, 1.3, 0.9, 0.6, 0.25],
            doses
        )
        subj3 = SubjectData(
            "SUBJ003",
            [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
            [1.8, 1.7, 1.4, 1.1, 0.7, 0.45, 0.18],
            doses
        )

        observed = ObservedData(
            [subj1, subj2, subj3];
            study_id="TEST001",
            analyte="DRUG1",
            units="ng/mL",
            time_units="h"
        )

        # Create model for simulation
        model_params = OneCompIVBolusParams(10.0, 50.0)  # CL, V
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

        iiv_spec = IIVSpec(LogNormalIIV(), Dict(:CL => 0.3, :V => 0.2), UInt64(42), 50)

        pop_spec = PopulationSpec(model_spec, iiv_spec, nothing, nothing, IndividualCovariates[])

        grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        config = VPCConfig(
            pi_levels=[0.10, 0.50, 0.90],
            binning=QuantileBinning(5),
            n_simulations=20,  # Low for testing
            n_bootstrap=100,   # Minimum required
            seed=UInt64(123)
        )

        # Compute VPC
        vpc_result = compute_vpc(observed, pop_spec, grid, solver; config=config)

        @test vpc_result.n_subjects_observed == 3
        @test vpc_result.n_observations_observed == 21
        @test vpc_result.n_simulations == 20
        @test length(vpc_result.bins) == 5

        # Check that observed percentiles are computed
        for bin in vpc_result.bins
            @test bin.n_observed > 0
            for p in bin.percentiles
                @test !isnan(p.observed) || bin.n_observed == 0
            end
        end
    end

    @testset "VPC with Residual Error" begin
        # Create observed data
        doses = [DoseEvent(0.0, 100.0)]
        subj1 = SubjectData(
            "SUBJ001",
            [0.0, 2.0, 4.0, 8.0, 12.0],
            [2.0, 1.5, 1.0, 0.5, 0.25],
            doses
        )

        observed = ObservedData([subj1])

        # Create model
        model_params = OneCompIVBolusParams(10.0, 50.0)  # CL, V
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

        iiv_spec = IIVSpec(LogNormalIIV(), Dict(:CL => 0.3), UInt64(42), 20)

        pop_spec = PopulationSpec(model_spec, iiv_spec, nothing, nothing, IndividualCovariates[])

        grid = SimGrid(0.0, 12.0, collect(0.0:2.0:12.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        # Proportional error model
        error_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),  # sigma
            :conc,
            UInt64(999)
        )

        config = VPCConfig(
            n_simulations=10,
            n_bootstrap=100,
            seed=UInt64(123)
        )

        vpc_result = compute_vpc(observed, pop_spec, grid, solver;
            config=config, error_spec=error_spec)

        @test vpc_result.n_simulations == 10
        @test length(vpc_result.bins) > 0
    end

    @testset "pcVPC" begin
        # Create observed data
        doses = [DoseEvent(0.0, 100.0)]
        subj1 = SubjectData(
            "SUBJ001",
            [0.0, 2.0, 4.0, 8.0],
            [2.0, 1.5, 1.0, 0.5],
            doses
        )

        observed = ObservedData([subj1])

        # Create model
        model_params = OneCompIVBolusParams(10.0, 50.0)  # CL, V
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

        iiv_spec = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2), UInt64(42), 20)

        pop_spec = PopulationSpec(model_spec, iiv_spec, nothing, nothing, IndividualCovariates[])

        grid = SimGrid(0.0, 8.0, collect(0.0:2.0:8.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        config = VPCConfig(
            n_simulations=10,
            n_bootstrap=100,
            seed=UInt64(123)
        )

        vpc_result = compute_pcvpc(observed, pop_spec, grid, solver; config=config)

        @test vpc_result.config.prediction_corrected == true
        @test length(vpc_result.bins) > 0
    end

end
