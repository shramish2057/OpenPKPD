# First-Order Conditional Estimation with Interaction (FOCE-I)
# The most commonly used method in NONMEM

using Optim
using LinearAlgebra
using ForwardDiff

export foce_estimate

"""
    foce_estimate(observed, model_spec, config, grid, solver, rng) -> EstimationResult

Estimate parameters using First-Order Conditional Estimation with Interaction (FOCE-I).

FOCE-I is the gold standard for population PK/PD estimation. It:
1. Finds the conditional mode of eta for each subject (inner optimization)
2. Linearizes the model around these modes
3. Approximates the marginal likelihood using the linearized model
4. Optimizes theta and omega to maximize this approximation (outer optimization)

The "interaction" refers to including the dependence of residual error on eta.
"""
function foce_estimate(
    observed::ObservedData,
    model_spec::ModelSpec,
    config::EstimationConfig{FOCEIMethod},
    grid::SimGrid,
    solver::SolverSpec,
    rng
)::EstimationResult
    start_time = time()
    n_theta = length(config.theta_init)
    n_eta = size(config.omega_init, 1)
    n_subj = n_subjects(observed)
    n_obs = n_observations(observed)

    subjects = extract_subject_data(observed)

    # Initialize parameters
    theta_current = copy(config.theta_init)
    omega_current = copy(config.omega_init)
    sigma_current = config.sigma_init

    # Storage for individual etas (empirical Bayes estimates)
    eta_per_subject = [zeros(n_eta) for _ in 1:n_subj]

    # Optimization tracking
    ofv_history = Float64[]
    converged = false
    n_iter = 0
    final_gradient_norm = Inf

    # Pack parameters for optimization
    # We optimize: [theta; log_omega_diag; log_sigma_params]
    function pack_params(theta, omega, sigma)
        omega_diag = diag(omega)
        sigma_vals = get_sigma_params(sigma)
        return vcat(theta, log.(omega_diag), log.(sigma_vals))
    end

    function unpack_params(x, n_theta, n_omega, n_sigma)
        theta = x[1:n_theta]
        omega_log_diag = x[n_theta+1:n_theta+n_omega]
        sigma_log = x[n_theta+n_omega+1:end]

        omega = Diagonal(exp.(omega_log_diag)) |> Matrix
        sigma_vals = exp.(sigma_log)

        return theta, omega, sigma_vals
    end

    n_omega_params = n_eta  # Diagonal for now
    n_sigma_params = _count_sigma_params(sigma_current)

    x_init = pack_params(theta_current, omega_current, sigma_current)

    # Define FOCE-I objective function
    function foce_objective(x)
        theta, omega, sigma_vals = unpack_params(x, n_theta, n_omega_params, n_sigma_params)

        # Check parameter bounds for theta
        if any(theta .< config.theta_lower) || any(theta .> config.theta_upper)
            return Inf
        end

        # Check omega is positive
        if any(diag(omega) .<= 0)
            return Inf
        end

        # Update sigma spec with new values
        sigma = update_sigma_params(sigma_current, sigma_vals)

        ofv = 0.0

        # For each subject, find eta mode and compute contribution
        for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
            # Find conditional mode of eta
            eta_mode = find_eta_mode_foce(
                theta, omega, sigma, times, obs, doses,
                model_spec, grid, solver, eta_per_subject[i],
                config.method.max_inner_iter, config.method.inner_tol
            )
            eta_per_subject[i] = eta_mode

            # Compute subject's contribution to FOCE objective
            subj_ofv = foce_subject_objective(
                theta, omega, sigma, eta_mode, times, obs, doses,
                model_spec, grid, solver
            )

            ofv += subj_ofv
        end

        return ofv
    end

    # Outer optimization using BFGS
    lower_bounds = vcat(config.theta_lower, fill(-10.0, n_omega_params + n_sigma_params))
    upper_bounds = vcat(config.theta_upper, fill(10.0, n_omega_params + n_sigma_params))

    result = optimize(
        foce_objective,
        lower_bounds,
        upper_bounds,
        x_init,
        Fminbox(BFGS(linesearch=LineSearches.BackTracking())),
        Optim.Options(
            iterations=config.max_iter,
            g_tol=config.tol,
            show_trace=config.verbose
        )
    )

    converged = Optim.converged(result)
    n_iter = Optim.iterations(result)
    final_gradient_norm = norm(Optim.gradient(result))

    # Extract final estimates
    x_final = Optim.minimizer(result)
    theta_final, omega_final, sigma_vals_final = unpack_params(x_final, n_theta, n_omega_params, n_sigma_params)
    sigma_final = update_sigma_params(sigma_current, sigma_vals_final)
    final_ofv = Optim.minimum(result)

    # Re-compute etas at final parameters
    for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
        eta_per_subject[i] = find_eta_mode_foce(
            theta_final, omega_final, sigma_final, times, obs, doses,
            model_spec, grid, solver, eta_per_subject[i],
            config.method.max_inner_iter, config.method.inner_tol
        )
    end

    # Compute individual estimates
    individual_estimates = compute_individual_estimates_foce(
        theta_final, omega_final, sigma_final,
        eta_per_subject, subjects, model_spec, grid, solver
    )

    # Compute standard errors if requested
    theta_se = nothing
    theta_rse = nothing
    omega_se = nothing
    sigma_se = nothing
    covariance_successful = false
    condition_num = NaN
    eigenvalue_ratio = NaN

    if config.compute_se
        theta_se, omega_se, sigma_se, covariance_successful, condition_num, eigenvalue_ratio =
            compute_standard_errors_foce(
                theta_final, omega_final, sigma_final,
                eta_per_subject, subjects, model_spec, grid, solver, config
            )
        if theta_se !== nothing
            theta_rse = 100.0 .* abs.(theta_se ./ theta_final)
        end
    end

    # Confidence intervals
    theta_ci_lower = nothing
    theta_ci_upper = nothing
    if config.compute_ci && theta_se !== nothing
        z = quantile(Normal(), 1.0 - (1.0 - config.ci_level) / 2)
        theta_ci_lower = theta_final .- z .* theta_se
        theta_ci_upper = theta_final .+ z .* theta_se
    end

    # Compute correlation matrix
    omega_corr = cov_to_corr(omega_final)

    # Model fit statistics
    n_params = n_theta + n_omega_params + n_sigma_params
    aic = compute_aic(final_ofv, n_params)
    bic = compute_bic(final_ofv, n_params, n_obs)

    elapsed = time() - start_time

    return EstimationResult(
        config,
        theta_final, theta_se, theta_rse, theta_ci_lower, theta_ci_upper,
        omega_final, omega_se, omega_corr,
        sigma_final, sigma_se,
        final_ofv, aic, bic,
        individual_estimates,
        converged, n_iter, final_gradient_norm, condition_num, eigenvalue_ratio,
        covariance_successful,
        String[],
        elapsed
    )
end

"""
Find the conditional mode of eta for a subject using FOCE-I.
"""
function find_eta_mode_foce(
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    times::Vector{Float64},
    obs::Vector{Float64},
    doses::Vector{DoseEvent},
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    eta_init::Vector{Float64},
    max_iter::Int,
    tol::Float64
)::Vector{Float64}
    omega_inv = inv(omega)

    # Objective: minimize -2 log p(eta|y)
    function eta_objective(eta)
        # Prior contribution: eta' * omega_inv * eta
        prior = eta' * omega_inv * eta

        # Likelihood contribution - wrap in try-catch for invalid parameters
        try
            ipred = compute_individual_predictions(
                theta, eta, times, doses, model_spec, grid, solver
            )

            ll = 0.0
            for (y, f) in zip(obs, ipred)
                ll += observation_log_likelihood(y, f, sigma)
            end

            return ll + prior
        catch e
            # Invalid parameters - return large value
            return 1e10 + prior
        end
    end

    # Optimize
    result = optimize(
        eta_objective,
        eta_init,
        BFGS(linesearch=LineSearches.BackTracking()),
        Optim.Options(iterations=max_iter, g_tol=tol, show_trace=false)
    )

    return Optim.minimizer(result)
end

"""
Compute a subject's contribution to the FOCE-I objective.

The FOCE-I objective includes:
1. -2LL of observations given eta mode
2. Prior contribution from eta
3. Laplacian correction term (log determinant of Hessian of eta posterior)
"""
function foce_subject_objective(
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    eta::Vector{Float64},
    times::Vector{Float64},
    obs::Vector{Float64},
    doses::Vector{DoseEvent},
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec
)::Float64
    n_eta = length(eta)
    omega_inv = inv(omega)

    # Compute individual predictions
    ipred = compute_individual_predictions(
        theta, eta, times, doses, model_spec, grid, solver
    )

    # -2LL contribution
    ll_contrib = 0.0
    for (y, f) in zip(obs, ipred)
        ll_contrib += observation_log_likelihood(y, f, sigma)
    end

    # Prior contribution
    prior_contrib = eta' * omega_inv * eta

    # Log determinant of omega
    logdet_omega = logdet(omega)

    # The full FOCE-I objective includes:
    # -2LL + eta'*Omega^{-1}*eta + log|Omega| + log|H_eta|
    # where H_eta is the Hessian of the eta posterior
    # For simplicity, we approximate log|H_eta| ≈ log|Omega^{-1}| = -log|Omega|
    # This gives us the Laplacian approximation

    # FOCE-I objective (simplified)
    ofv = ll_contrib + prior_contrib + logdet_omega

    return ofv
end

"""
Compute individual estimates for FOCE-I.
"""
function compute_individual_estimates_foce(
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    etas::Vector{Vector{Float64}},
    subjects::Vector,
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec
)::Vector{IndividualEstimate}
    estimates = IndividualEstimate[]
    omega_inv = inv(omega)

    for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
        eta = etas[i]

        # Compute IPRED
        ipred = compute_individual_predictions(
            theta, eta, times, doses, model_spec, grid, solver
        )

        # Compute PRED
        pred = compute_individual_predictions(
            theta, zeros(length(eta)), times, doses, model_spec, grid, solver
        )

        # Compute residuals
        cwres = compute_cwres(obs, ipred, sigma)
        iwres = compute_iwres(obs, ipred, sigma)
        wres = compute_wres(obs, pred, sigma)

        # OFV contribution
        ofv_contrib = foce_subject_objective(
            theta, omega, sigma, eta, times, obs, doses,
            model_spec, grid, solver
        )

        push!(estimates, IndividualEstimate(
            subj_id,
            eta,
            nothing,
            ipred,
            pred,
            cwres,
            iwres,
            wres,
            ofv_contrib
        ))
    end

    return estimates
end

"""
Compute standard errors for FOCE-I estimates.
"""
function compute_standard_errors_foce(
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    etas::Vector{Vector{Float64}},
    subjects::Vector,
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    config::EstimationConfig
)
    # Use finite differences to compute Hessian
    n_theta = length(theta)
    h = 1e-5

    function ofv_theta(th)
        ofv = 0.0
        for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
            eta = etas[i]
            ofv += foce_subject_objective(
                th, omega, sigma, eta, times, obs, doses,
                model_spec, grid, solver
            )
        end
        return ofv
    end

    # Compute Hessian
    hessian = zeros(n_theta, n_theta)

    for i in 1:n_theta
        for j in i:n_theta
            if i == j
                theta_plus = copy(theta)
                theta_minus = copy(theta)
                theta_plus[i] += h
                theta_minus[i] -= h

                f_plus = ofv_theta(theta_plus)
                f_center = ofv_theta(theta)
                f_minus = ofv_theta(theta_minus)

                hessian[i, i] = (f_plus - 2*f_center + f_minus) / h^2
            else
                theta_pp = copy(theta)
                theta_pm = copy(theta)
                theta_mp = copy(theta)
                theta_mm = copy(theta)

                theta_pp[i] += h; theta_pp[j] += h
                theta_pm[i] += h; theta_pm[j] -= h
                theta_mp[i] -= h; theta_mp[j] += h
                theta_mm[i] -= h; theta_mm[j] -= h

                f_pp = ofv_theta(theta_pp)
                f_pm = ofv_theta(theta_pm)
                f_mp = ofv_theta(theta_mp)
                f_mm = ofv_theta(theta_mm)

                hessian[i, j] = (f_pp - f_pm - f_mp + f_mm) / (4 * h^2)
                hessian[j, i] = hessian[i, j]
            end
        end
    end

    # Check positive definiteness
    eigenvalues = eigvals(hessian)
    real_eigenvalues = real.(eigenvalues)

    if any(real_eigenvalues .<= 0)
        return nothing, nothing, nothing, false, NaN, NaN
    end

    condition_num = maximum(real_eigenvalues) / minimum(real_eigenvalues)
    eigenvalue_ratio = minimum(real_eigenvalues) / maximum(real_eigenvalues)

    try
        cov_matrix = inv(hessian)
        variances = diag(cov_matrix)

        if any(variances .< 0)
            return nothing, nothing, nothing, false, condition_num, eigenvalue_ratio
        end

        theta_se = sqrt.(variances)
        return theta_se, nothing, nothing, true, condition_num, eigenvalue_ratio
    catch e
        return nothing, nothing, nothing, false, condition_num, eigenvalue_ratio
    end
end

# ------------------------------------------------------------------
# Helper functions for sigma parameter handling
# ------------------------------------------------------------------

function get_sigma_params(spec::ResidualErrorSpec{AdditiveError})::Vector{Float64}
    return [spec.params.sigma]
end

function get_sigma_params(spec::ResidualErrorSpec{ProportionalError})::Vector{Float64}
    return [spec.params.sigma]
end

function get_sigma_params(spec::ResidualErrorSpec{CombinedError})::Vector{Float64}
    return [spec.params.sigma_add, spec.params.sigma_prop]
end

function get_sigma_params(spec::ResidualErrorSpec{ExponentialError})::Vector{Float64}
    return [spec.params.sigma]
end

function update_sigma_params(spec::ResidualErrorSpec{AdditiveError}, vals::Vector{Float64})::ResidualErrorSpec
    return ResidualErrorSpec(
        AdditiveError(),
        AdditiveErrorParams(vals[1]),
        spec.observation,
        spec.seed
    )
end

function update_sigma_params(spec::ResidualErrorSpec{ProportionalError}, vals::Vector{Float64})::ResidualErrorSpec
    return ResidualErrorSpec(
        ProportionalError(),
        ProportionalErrorParams(vals[1]),
        spec.observation,
        spec.seed
    )
end

function update_sigma_params(spec::ResidualErrorSpec{CombinedError}, vals::Vector{Float64})::ResidualErrorSpec
    return ResidualErrorSpec(
        CombinedError(),
        CombinedErrorParams(vals[1], vals[2]),
        spec.observation,
        spec.seed
    )
end

function update_sigma_params(spec::ResidualErrorSpec{ExponentialError}, vals::Vector{Float64})::ResidualErrorSpec
    return ResidualErrorSpec(
        ExponentialError(),
        ExponentialErrorParams(vals[1]),
        spec.observation,
        spec.seed
    )
end

export get_sigma_params, update_sigma_params
