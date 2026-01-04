# Standard Error Computation for Parameter Estimates
# Implements various methods for computing SEs and confidence intervals

using LinearAlgebra
using Distributions

export compute_se_from_hessian, compute_se_sandwich, compute_ci

"""
Compute standard errors from the Hessian (observed Fisher information).

Arguments:
- hessian: Hessian matrix of the objective function at the minimum
- method: :inverse (standard) or :robust (use eigenvalue regularization)

Returns:
- (se, cov_matrix, success) tuple
"""
function compute_se_from_hessian(
    hessian::Matrix{Float64};
    method::Symbol=:inverse,
    min_eigenvalue::Float64=1e-8
)::Tuple{Union{Nothing,Vector{Float64}},Union{Nothing,Matrix{Float64}},Bool}
    n = size(hessian, 1)

    # Check symmetry
    if !issymmetric(hessian)
        hessian = (hessian + hessian') / 2
    end

    # Check positive definiteness
    eigenvalues = eigvals(Symmetric(hessian))

    if any(real.(eigenvalues) .<= 0)
        if method == :robust
            # Regularize eigenvalues
            eigvecs = eigvecs(Symmetric(hessian))
            eigenvalues = max.(real.(eigenvalues), min_eigenvalue)
            hessian = eigvecs * Diagonal(eigenvalues) * eigvecs'
        else
            return (nothing, nothing, false)
        end
    end

    # Invert Hessian
    try
        cov_matrix = inv(hessian)
        variances = diag(cov_matrix)

        if any(variances .< 0)
            return (nothing, nothing, false)
        end

        se = sqrt.(variances)
        return (se, cov_matrix, true)
    catch e
        return (nothing, nothing, false)
    end
end

"""
Compute sandwich (robust) standard errors.

The sandwich estimator is: cov = H^{-1} * G * H^{-1}
where H is the Hessian and G is the outer product of gradients.

Arguments:
- hessian: Hessian matrix
- gradients: Matrix of individual gradient contributions (n_obs x n_params)

Returns:
- (se, cov_matrix, success) tuple
"""
function compute_se_sandwich(
    hessian::Matrix{Float64},
    gradients::Matrix{Float64}
)::Tuple{Union{Nothing,Vector{Float64}},Union{Nothing,Matrix{Float64}},Bool}
    n_params = size(hessian, 1)

    # Compute H^{-1}
    try
        H_inv = inv(hessian)
    catch
        return (nothing, nothing, false)
    end

    # Compute meat matrix: G = sum of g_i * g_i'
    G = zeros(n_params, n_params)
    for i in 1:size(gradients, 1)
        g = gradients[i, :]
        G += g * g'
    end

    # Sandwich: H^{-1} * G * H^{-1}
    try
        cov_matrix = H_inv * G * H_inv
        variances = diag(cov_matrix)

        if any(variances .< 0)
            return (nothing, nothing, false)
        end

        se = sqrt.(variances)
        return (se, cov_matrix, true)
    catch e
        return (nothing, nothing, false)
    end
end

"""
Compute confidence intervals from estimates and standard errors.

Arguments:
- estimates: Parameter estimates
- se: Standard errors
- level: Confidence level (default: 0.95)
- method: :wald (normal approximation) or :profile (profile likelihood)

Returns:
- (lower, upper) bounds vectors
"""
function compute_ci(
    estimates::Vector{Float64},
    se::Vector{Float64};
    level::Float64=0.95,
    method::Symbol=:wald
)::Tuple{Vector{Float64},Vector{Float64}}
    if method == :wald
        z = quantile(Normal(), 1.0 - (1.0 - level) / 2)
        lower = estimates .- z .* se
        upper = estimates .+ z .* se
        return (lower, upper)
    else
        error("Method $method not implemented")
    end
end

"""
Compute relative standard errors (coefficient of variation).

Arguments:
- estimates: Parameter estimates
- se: Standard errors

Returns:
- RSE as percentages
"""
function compute_rse(
    estimates::Vector{Float64},
    se::Vector{Float64}
)::Vector{Float64}
    return 100.0 .* abs.(se ./ estimates)
end

export compute_rse

"""
Compute correlation matrix from covariance matrix.
"""
function cov_to_corr_matrix(cov::Matrix{Float64})::Matrix{Float64}
    d = sqrt.(diag(cov))
    n = length(d)
    corr = zeros(n, n)

    for i in 1:n
        for j in 1:n
            if d[i] > 0 && d[j] > 0
                corr[i, j] = cov[i, j] / (d[i] * d[j])
            elseif i == j
                corr[i, j] = 1.0
            end
        end
    end

    return corr
end

export cov_to_corr_matrix

"""
Bootstrap standard errors.

Arguments:
- data: Vector of observations
- statistic: Function that computes the statistic from data
- n_bootstrap: Number of bootstrap samples
- rng: Random number generator

Returns:
- Bootstrap SE estimate
"""
function bootstrap_se(
    data::Vector{Float64},
    statistic::Function,
    n_bootstrap::Int,
    rng
)::Float64
    n = length(data)
    bootstrap_stats = Float64[]

    for _ in 1:n_bootstrap
        indices = rand(rng, 1:n, n)
        sample = data[indices]
        push!(bootstrap_stats, statistic(sample))
    end

    return std(bootstrap_stats)
end

export bootstrap_se

# Simple std function
function std(x::Vector{Float64})::Float64
    n = length(x)
    if n < 2
        return 0.0
    end
    m = sum(x) / n
    return sqrt(sum((x .- m).^2) / (n - 1))
end
