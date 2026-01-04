# Test suite for Parameter Estimation (NLME)

using Test
using OpenPKPDCore
using LinearAlgebra
using StableRNGs

@testset "Estimation" begin

    @testset "Estimation Types" begin
        @testset "FOCEIMethod" begin
            method = FOCEIMethod()
            @test method.max_inner_iter == 50
            @test method.inner_tol == 1e-6
            @test method.centered == false

            method2 = FOCEIMethod(max_inner_iter=100, inner_tol=1e-8, centered=true)
            @test method2.max_inner_iter == 100
            @test method2.inner_tol == 1e-8
            @test method2.centered == true
        end

        @testset "SAEMMethod" begin
            method = SAEMMethod()
            @test method.n_burn == 300
            @test method.n_iter == 200
            @test method.n_chains == 3
            @test method.step_size_schedule == :harmonic

            method2 = SAEMMethod(n_burn=100, n_iter=50, step_size_schedule=:constant)
            @test method2.n_burn == 100
            @test method2.step_size_schedule == :constant
        end

        @testset "LaplacianMethod" begin
            method = LaplacianMethod()
            @test method.max_inner_iter == 50
            @test method.inner_tol == 1e-6
        end

        @testset "OmegaStructure" begin
            @test DiagonalOmega() isa OmegaStructure
            @test BlockOmega([2, 1]) isa OmegaStructure
            @test FullOmega() isa OmegaStructure
        end

        @testset "EstimationConfig" begin
            sigma_spec = ResidualErrorSpec(
                ProportionalError(),
                ProportionalErrorParams(0.1),
                :conc,
                UInt64(1)
            )

            config = EstimationConfig(
                FOCEIMethod();
                theta_init=[10.0, 50.0],
                omega_init=diagm([0.09, 0.04]),
                sigma_init=sigma_spec
            )

            @test length(config.theta_init) == 2
            @test size(config.omega_init) == (2, 2)
            @test config.max_iter == 500
            @test config.compute_se == true
        end
    end

    @testset "Utility Functions" begin
        @testset "ensure_positive_definite" begin
            # Already positive definite
            A = [1.0 0.5; 0.5 1.0]
            A_pd = ensure_positive_definite(A)
            @test isposdef(Symmetric(A_pd))

            # Not positive definite
            B = [1.0 2.0; 2.0 1.0]
            B_pd = ensure_positive_definite(B)
            @test isposdef(Symmetric(B_pd))
        end

        @testset "cov_to_corr" begin
            cov = [0.09 0.03; 0.03 0.04]
            corr = cov_to_corr(cov)

            @test corr[1, 1] ≈ 1.0
            @test corr[2, 2] ≈ 1.0
            @test abs(corr[1, 2]) <= 1.0
            @test corr[1, 2] ≈ corr[2, 1]
        end
    end

    @testset "Gradient Utilities" begin
        @testset "gradient_fd" begin
            f(x) = x[1]^2 + x[2]^2
            x = [1.0, 2.0]

            grad = gradient_fd(f, x)
            @test length(grad) == 2
            @test isapprox(grad[1], 2.0; atol=1e-4)
            @test isapprox(grad[2], 4.0; atol=1e-4)
        end

        @testset "hessian_fd" begin
            f(x) = x[1]^2 + x[2]^2
            x = [1.0, 2.0]

            hess = hessian_fd(f, x)
            @test size(hess) == (2, 2)
            @test isapprox(hess[1, 1], 2.0; atol=1e-3)
            @test isapprox(hess[2, 2], 2.0; atol=1e-3)
            @test isapprox(hess[1, 2], 0.0; atol=1e-3)
        end

        @testset "compute_gradient (ForwardDiff)" begin
            f(x) = x[1]^2 + 2*x[1]*x[2] + x[2]^2
            x = [1.0, 2.0]

            grad = compute_gradient(f, x)
            @test isapprox(grad[1], 6.0; atol=1e-10)
            @test isapprox(grad[2], 6.0; atol=1e-10)
        end
    end

    @testset "Diagnostics" begin
        @testset "compute_wres" begin
            obs = [1.0, 2.0, 3.0]
            pred = [1.1, 1.9, 3.1]
            sigma = ResidualErrorSpec(
                ProportionalError(),
                ProportionalErrorParams(0.1),
                :conc,
                UInt64(1)
            )

            wres = compute_wres(obs, pred, sigma)
            @test length(wres) == 3
            @test all(!isnan(w) for w in wres)
        end

        @testset "shrinkage_eta" begin
            # Create mock etas with some variability
            etas = [
                [0.1, 0.05],
                [-0.1, -0.05],
                [0.05, 0.02],
                [-0.05, -0.02]
            ]
            omega = diagm([0.09, 0.04])  # 30%, 20% CV

            shrinkage = shrinkage_eta(etas, omega)
            @test length(shrinkage) == 2
            # Shrinkage should be between 0 and 1 for reasonable data
            @test all(0 <= s <= 1 for s in shrinkage if !isnan(s))
        end

        @testset "shrinkage_epsilon" begin
            # IWRES should be approximately N(0,1) for good model
            iwres = randn(StableRNG(42), 100)

            eps_shrinkage = shrinkage_epsilon(iwres)
            # For N(0,1) data, shrinkage should be near 0
            @test abs(eps_shrinkage) < 0.3
        end

        @testset "vpc_check" begin
            obs = [1.0, 2.0, 3.0, 4.0, 5.0]
            pi_lower = [0.5, 1.5, 2.5, 3.5, 4.5]
            pi_upper = [1.5, 2.5, 3.5, 4.5, 5.5]

            pct = vpc_check(obs, pi_lower, pi_upper)
            @test pct == 100.0  # All should be within PI
        end

        @testset "likelihood_ratio_test" begin
            ofv_full = 100.0
            ofv_reduced = 105.0
            df = 1

            chi_sq, p_value = likelihood_ratio_test(ofv_full, ofv_reduced, df)
            @test chi_sq ≈ 5.0
            @test 0 < p_value < 1
        end

        @testset "residual_summary" begin
            residuals = [0.1, -0.1, 0.2, -0.2, 0.0]
            summary = residual_summary(residuals)

            @test haskey(summary, :mean)
            @test haskey(summary, :sd)
            @test haskey(summary, :median)
            @test isapprox(summary[:mean], 0.0; atol=1e-10)
        end
    end

    @testset "Standard Errors" begin
        @testset "compute_se_from_hessian" begin
            # Positive definite Hessian
            hessian = [4.0 0.0; 0.0 4.0]
            se, cov, success = compute_se_from_hessian(hessian)

            @test success
            @test se !== nothing
            @test isapprox(se[1], 0.5; atol=1e-10)
            @test isapprox(se[2], 0.5; atol=1e-10)
        end

        @testset "compute_ci" begin
            estimates = [10.0, 50.0]
            se = [1.0, 5.0]

            lower, upper = compute_ci(estimates, se; level=0.95)

            @test all(lower .< estimates)
            @test all(upper .> estimates)
            @test isapprox(estimates .- lower, upper .- estimates; atol=1e-10)
        end

        @testset "compute_rse" begin
            estimates = [10.0, 50.0]
            se = [1.0, 5.0]

            rse = compute_rse(estimates, se)
            @test isapprox(rse[1], 10.0; atol=1e-10)
            @test isapprox(rse[2], 10.0; atol=1e-10)
        end
    end

    @testset "Laplacian Estimation" begin
        # Create simple test data
        doses = [DoseEvent(0.0, 100.0)]

        # Create a few subjects with observations
        subj1 = SubjectData(
            "SUBJ001",
            [0.5, 1.0, 2.0, 4.0, 8.0],
            [1.8, 1.6, 1.3, 0.9, 0.4],
            doses
        )
        subj2 = SubjectData(
            "SUBJ002",
            [0.5, 1.0, 2.0, 4.0, 8.0],
            [2.2, 1.9, 1.5, 1.0, 0.5],
            doses
        )
        subj3 = SubjectData(
            "SUBJ003",
            [0.5, 1.0, 2.0, 4.0, 8.0],
            [1.9, 1.7, 1.4, 0.95, 0.45],
            doses
        )

        observed = ObservedData([subj1, subj2, subj3])

        # Model specification
        model_params = OneCompIVBolusParams(10.0, 50.0)  # CL, V
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

        grid = SimGrid(0.0, 12.0, collect(0.0:0.5:12.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        # Estimation config
        sigma_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        config = EstimationConfig(
            LaplacianMethod(max_inner_iter=20, inner_tol=1e-4);
            theta_init=[10.0, 50.0],
            theta_lower=[1.0, 10.0],
            theta_upper=[100.0, 200.0],
            theta_names=[:CL, :V],
            omega_init=diagm([0.09, 0.04]),
            omega_names=[:eta_CL, :eta_V],
            sigma_init=sigma_spec,
            max_iter=10,  # Few iterations for testing
            tol=1e-2,
            compute_se=false,  # Skip SE for faster test
            verbose=false
        )

        # Run estimation
        result = estimate(observed, model_spec, config; grid=grid, solver=solver)

        @test result isa EstimationResult
        @test length(result.theta) == 2
        @test result.theta[1] > 0  # CL should be positive
        @test result.theta[2] > 0  # V should be positive
        @test result.n_iterations > 0
        @test isfinite(result.ofv)
        @test length(result.individuals) == 3
    end

    @testset "Individual Estimates" begin
        # Create test data
        doses = [DoseEvent(0.0, 100.0)]

        subj1 = SubjectData(
            "SUBJ001",
            [0.5, 1.0, 2.0, 4.0],
            [1.8, 1.6, 1.3, 0.9],
            doses
        )

        observed = ObservedData([subj1])

        model_params = OneCompIVBolusParams(10.0, 50.0)
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

        grid = SimGrid(0.0, 8.0, collect(0.0:0.5:8.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sigma_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        config = EstimationConfig(
            LaplacianMethod(max_inner_iter=10, inner_tol=1e-3);
            theta_init=[10.0, 50.0],
            omega_init=diagm([0.09, 0.04]),
            sigma_init=sigma_spec,
            max_iter=5,
            compute_se=false,
            verbose=false
        )

        result = estimate(observed, model_spec, config; grid=grid, solver=solver)

        # Check individual estimates
        @test length(result.individuals) == 1
        ind = result.individuals[1]

        @test ind.subject_id == "SUBJ001"
        @test length(ind.eta) == 2
        @test length(ind.ipred) == 4
        @test length(ind.pred) == 4
        @test length(ind.cwres) == 4
        @test isfinite(ind.ofv_contribution)
    end

    @testset "Model Comparison Metrics" begin
        @testset "compute_aic" begin
            ofv = 100.0
            n_params = 5
            aic = compute_aic(ofv, n_params)
            @test aic == 110.0
        end

        @testset "compute_bic" begin
            ofv = 100.0
            n_params = 5
            n_obs = 100
            bic = compute_bic(ofv, n_params, n_obs)
            @test bic ≈ 100.0 + log(100) * 5
        end
    end

    @testset "Parameter Count" begin
        doses = [DoseEvent(0.0, 100.0)]
        sigma_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        # Test _count_omega_params
        @test OpenPKPDCore._count_omega_params(DiagonalOmega(), 3) == 3
        @test OpenPKPDCore._count_omega_params(FullOmega(), 3) == 6  # 3*4/2
        @test OpenPKPDCore._count_omega_params(BlockOmega([2, 1]), 3) == 4  # 2*3/2 + 1

        # Test _count_sigma_params
        @test OpenPKPDCore._count_sigma_params(sigma_spec) == 1

        combined_spec = ResidualErrorSpec(
            CombinedError(),
            CombinedErrorParams(0.1, 0.05),
            :conc,
            UInt64(1)
        )
        @test OpenPKPDCore._count_sigma_params(combined_spec) == 2
    end

end
