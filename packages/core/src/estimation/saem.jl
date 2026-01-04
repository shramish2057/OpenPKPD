# Stochastic Approximation Expectation Maximization (SAEM)
# More robust than FOCE-I for complex models

using LinearAlgebra
using StableRNGs
using Distributions

export saem_estimate

"""
    saem_estimate(observed, model_spec, config, grid, solver, rng) -> EstimationResult

Estimate parameters using Stochastic Approximation Expectation Maximization (SAEM).

SAEM is often more robust than FOCE-I for:
- Models with high inter-individual variability
- Complex nonlinear models
- Sparse data

The algorithm alternates between:
- E-step: MCMC sampling of eta from posterior
- M-step: Stochastic approximation update of theta, omega, sigma
"""
function saem_estimate(
    observed::ObservedData,
    model_spec::ModelSpec,
    config::EstimationConfig{SAEMMethod},
    grid::SimGrid,
    solver::SolverSpec,
    rng
)::EstimationResult
    start_time = time()
    n_theta = length(config.theta_init)
    n_eta = size(config.omega_init, 1)
    n_subj = n_subjects(observed)
    n_obs_total = n_observations(observed)

    subjects = extract_subject_data(observed)

    # Initialize parameters
    theta_current = copy(config.theta_init)
    omega_current = copy(config.omega_init)
    sigma_current = config.sigma_init

    # Initialize sufficient statistics
    S_theta = zeros(n_theta)
    S_omega = zeros(n_eta, n_eta)
    n_sigma = _count_sigma_params(sigma_current)
    S_sigma = zeros(n_sigma)

    # Initialize eta chains for each subject
    eta_chains = [[zeros(n_eta) for _ in 1:config.method.n_chains] for _ in 1:n_subj]

    # MCMC proposal variance (adaptive)
    proposal_sd = 0.1 * ones(n_eta)

    # Tracking
    ofv_history = Float64[]
    theta_history = Vector{Float64}[]
    n_total_iter = config.method.n_burn + config.method.n_iter

    if config.verbose
        println("SAEM: Starting $(config.method.n_burn) burn-in + $(config.method.n_iter) main iterations")
    end

    # Main SAEM loop
    for iter in 1:n_total_iter
        is_burn_in = iter <= config.method.n_burn

        # Step size (decreases after burn-in for stochastic approximation)
        if config.method.step_size_schedule == :harmonic
            if is_burn_in
                step_size = 1.0  # Full step during burn-in
            else
                step_size = 1.0 / (iter - config.method.n_burn)
            end
        else  # :constant
            step_size = is_burn_in ? 1.0 : 0.1
        end

        # E-step: Sample eta for each subject using MCMC
        sampled_etas = Vector{Vector{Float64}}(undef, n_subj)

        for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
            # Run MCMC chains
            eta_samples = mcmc_sample_eta(
                theta_current, omega_current, sigma_current,
                times, obs, doses, model_spec, grid, solver,
                eta_chains[i], proposal_sd, 10, rng  # 10 MCMC steps
            )

            # Update chains
            eta_chains[i] = eta_samples

            # Use last sample from first chain for statistics
            sampled_etas[i] = eta_samples[1]
        end

        # M-step: Update sufficient statistics

        # Sufficient statistic for omega: sum of eta * eta'
        S_omega_new = zeros(n_eta, n_eta)
        for eta in sampled_etas
            S_omega_new += eta * eta'
        end
        S_omega_new /= n_subj

        # Update omega sufficient statistic with stochastic approximation
        S_omega = (1 - step_size) * S_omega + step_size * S_omega_new

        # Update omega estimate
        omega_current = ensure_positive_definite(S_omega)

        # If diagonal structure, zero out off-diagonals
        if config.omega_structure isa DiagonalOmega
            omega_current = Diagonal(diag(omega_current)) |> Matrix
        end

        # Update theta using closed-form or optimization
        theta_new = update_theta_saem(
            theta_current, sampled_etas, subjects, sigma_current,
            model_spec, grid, solver, config, step_size
        )
        theta_current = theta_new

        # Update sigma (for proportional error, estimate from residuals)
        sigma_vals_new = update_sigma_saem(
            theta_current, sampled_etas, subjects, sigma_current,
            model_spec, grid, solver
        )

        sigma_vals_current = get_sigma_params(sigma_current)
        sigma_vals_updated = (1 - step_size) .* sigma_vals_current .+ step_size .* sigma_vals_new
        sigma_current = update_sigma_params(sigma_current, sigma_vals_updated)

        # Compute current OFV for monitoring
        if iter % 50 == 0 || iter == n_total_iter
            current_ofv = compute_saem_ofv(
                theta_current, omega_current, sigma_current,
                sampled_etas, subjects, model_spec, grid, solver
            )
            push!(ofv_history, current_ofv)
            push!(theta_history, copy(theta_current))

            if config.verbose
                phase = is_burn_in ? "burn-in" : "main"
                println("  Iter $iter ($phase): OFV=$(round(current_ofv, digits=3)), theta=$(round.(theta_current, digits=3))")
            end
        end
    end

    # Final OFV computation
    final_etas = [eta_chains[i][1] for i in 1:n_subj]
    final_ofv = compute_saem_ofv(
        theta_current, omega_current, sigma_current,
        final_etas, subjects, model_spec, grid, solver
    )

    # Compute individual estimates
    individual_estimates = compute_individual_estimates_saem(
        theta_current, omega_current, sigma_current,
        final_etas, subjects, model_spec, grid, solver
    )

    # Check convergence based on OFV stabilization
    converged = length(ofv_history) >= 2 &&
                abs(ofv_history[end] - ofv_history[end-1]) < config.tol * abs(ofv_history[end])

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
            compute_standard_errors_saem(
                theta_current, omega_current, sigma_current,
                final_etas, subjects, model_spec, grid, solver, config
            )
        if theta_se !== nothing
            theta_rse = 100.0 .* abs.(theta_se ./ theta_current)
        end
    end

    # Confidence intervals
    theta_ci_lower = nothing
    theta_ci_upper = nothing
    if config.compute_ci && theta_se !== nothing
        z = quantile(Normal(), 1.0 - (1.0 - config.ci_level) / 2)
        theta_ci_lower = theta_current .- z .* theta_se
        theta_ci_upper = theta_current .+ z .* theta_se
    end

    # Compute correlation matrix
    omega_corr = cov_to_corr(omega_current)

    # Model fit statistics
    n_params = n_theta + _count_omega_params(config.omega_structure, n_eta) +
               _count_sigma_params(sigma_current)
    aic = compute_aic(final_ofv, n_params)
    bic = compute_bic(final_ofv, n_params, n_obs_total)

    elapsed = time() - start_time

    return EstimationResult(
        config,
        theta_current, theta_se, theta_rse, theta_ci_lower, theta_ci_upper,
        omega_current, omega_se, omega_corr,
        sigma_current, sigma_se,
        final_ofv, aic, bic,
        individual_estimates,
        converged, n_total_iter, NaN, condition_num, eigenvalue_ratio,
        covariance_successful,
        String[],
        elapsed
    )
end

"""
MCMC sampling of eta using Metropolis-Hastings.
"""
function mcmc_sample_eta(
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    times::Vector{Float64},
    obs::Vector{Float64},
    doses::Vector{DoseEvent},
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    current_chains::Vector{Vector{Float64}},
    proposal_sd::Vector{Float64},
    n_steps::Int,
    rng
)::Vector{Vector{Float64}}
    n_chains = length(current_chains)
    n_eta = length(current_chains[1])
    omega_inv = inv(omega)

    # Log posterior for eta
    function log_posterior(eta)
        # Prior: N(0, omega)
        log_prior = -0.5 * eta' * omega_inv * eta - 0.5 * logdet(omega)

        # Likelihood
        ipred = compute_individual_predictions(
            theta, eta, times, doses, model_spec, grid, solver
        )

        log_lik = 0.0
        for (y, f) in zip(obs, ipred)
            var_res = residual_variance(f, sigma)
            if var_res > 0
                log_lik += -0.5 * (log(2π) + log(var_res) + (y - f)^2 / var_res)
            else
                log_lik = -Inf
                break
            end
        end

        return log_prior + log_lik
    end

    new_chains = Vector{Vector{Float64}}(undef, n_chains)

    for c in 1:n_chains
        eta_current = copy(current_chains[c])
        log_p_current = log_posterior(eta_current)

        for step in 1:n_steps
            # Propose new eta
            eta_proposed = eta_current .+ randn(rng, n_eta) .* proposal_sd
            log_p_proposed = log_posterior(eta_proposed)

            # Accept/reject
            log_alpha = log_p_proposed - log_p_current
            if log(rand(rng)) < log_alpha
                eta_current = eta_proposed
                log_p_current = log_p_proposed
            end
        end

        new_chains[c] = eta_current
    end

    return new_chains
end

"""
Update theta in the M-step of SAEM.
"""
function update_theta_saem(
    theta_current::Vector{Float64},
    etas::Vector{Vector{Float64}},
    subjects::Vector,
    sigma::ResidualErrorSpec,
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    config::EstimationConfig,
    step_size::Float64
)::Vector{Float64}
    # Use a few optimization steps to update theta
    function theta_objective(theta)
        # Check bounds
        if any(theta .< config.theta_lower) || any(theta .> config.theta_upper)
            return Inf
        end

        ll = 0.0
        for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
            eta = etas[i]
            ipred = compute_individual_predictions(
                theta, eta, times, doses, model_spec, grid, solver
            )

            for (y, f) in zip(obs, ipred)
                ll += observation_log_likelihood(y, f, sigma)
            end
        end

        return ll
    end

    # Take a bounded optimization step
    result = optimize(
        theta_objective,
        config.theta_lower,
        config.theta_upper,
        theta_current,
        Fminbox(BFGS(linesearch=LineSearches.BackTracking())),
        Optim.Options(iterations=10, show_trace=false)
    )

    theta_new = Optim.minimizer(result)

    # Apply stochastic approximation
    return (1 - step_size) .* theta_current .+ step_size .* theta_new
end

"""
Update sigma parameters in the M-step of SAEM.
"""
function update_sigma_saem(
    theta::Vector{Float64},
    etas::Vector{Vector{Float64}},
    subjects::Vector,
    sigma::ResidualErrorSpec,
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec
)::Vector{Float64}
    # Estimate sigma from residuals
    residuals_sq = Float64[]
    predictions = Float64[]

    for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
        eta = etas[i]
        ipred = compute_individual_predictions(
            theta, eta, times, doses, model_spec, grid, solver
        )

        for (y, f) in zip(obs, ipred)
            push!(residuals_sq, (y - f)^2)
            push!(predictions, f)
        end
    end

    # Estimate sigma based on error model type
    return estimate_sigma_from_residuals(residuals_sq, predictions, sigma)
end

function estimate_sigma_from_residuals(
    residuals_sq::Vector{Float64},
    predictions::Vector{Float64},
    spec::ResidualErrorSpec{AdditiveError}
)::Vector{Float64}
    # Additive: var(res) = sigma^2
    sigma = sqrt(mean(residuals_sq))
    return [max(sigma, 0.001)]
end

function estimate_sigma_from_residuals(
    residuals_sq::Vector{Float64},
    predictions::Vector{Float64},
    spec::ResidualErrorSpec{ProportionalError}
)::Vector{Float64}
    # Proportional: var(res) = (sigma * F)^2
    # sigma = sqrt(mean(res^2 / F^2))
    weighted_sq = [r / max(f^2, 1e-10) for (r, f) in zip(residuals_sq, predictions)]
    sigma = sqrt(mean(weighted_sq))
    return [max(sigma, 0.001)]
end

function estimate_sigma_from_residuals(
    residuals_sq::Vector{Float64},
    predictions::Vector{Float64},
    spec::ResidualErrorSpec{CombinedError}
)::Vector{Float64}
    # Combined: more complex, use simple heuristic
    sigma_add = sqrt(mean(residuals_sq)) * 0.3
    sigma_prop = sqrt(mean([r / max(f^2, 1e-10) for (r, f) in zip(residuals_sq, predictions)])) * 0.7
    return [max(sigma_add, 0.001), max(sigma_prop, 0.001)]
end

function estimate_sigma_from_residuals(
    residuals_sq::Vector{Float64},
    predictions::Vector{Float64},
    spec::ResidualErrorSpec{ExponentialError}
)::Vector{Float64}
    # Exponential: on log scale
    sigma = sqrt(mean(residuals_sq))
    return [max(sigma, 0.001)]
end

"""
Compute SAEM objective function value.
"""
function compute_saem_ofv(
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    etas::Vector{Vector{Float64}},
    subjects::Vector,
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec
)::Float64
    n_subj = length(subjects)
    omega_inv = inv(omega)
    ofv = 0.0

    for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
        eta = etas[i]

        # Prior contribution
        ofv += eta' * omega_inv * eta

        # Likelihood contribution
        ipred = compute_individual_predictions(
            theta, eta, times, doses, model_spec, grid, solver
        )

        for (y, f) in zip(obs, ipred)
            ofv += observation_log_likelihood(y, f, sigma)
        end
    end

    # Add log determinant term
    ofv += n_subj * logdet(omega)

    return ofv
end

"""
Compute individual estimates for SAEM.
"""
function compute_individual_estimates_saem(
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
        ofv_contrib = eta' * omega_inv * eta
        for (y, f) in zip(obs, ipred)
            ofv_contrib += observation_log_likelihood(y, f, sigma)
        end

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
Compute standard errors for SAEM (using Louis' method or Fisher information).
"""
function compute_standard_errors_saem(
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
    # Use numerical Hessian (similar to Laplacian/FOCE)
    n_theta = length(theta)
    h = 1e-5

    function ofv_theta(th)
        return compute_saem_ofv(th, omega, sigma, etas, subjects, model_spec, grid, solver)
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

# Helper function
function mean(x::Vector{Float64})::Float64
    return sum(x) / length(x)
end
