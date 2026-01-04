# Estimation Types for NLME Parameter Estimation
# Defines methods, configurations, and results

using LinearAlgebra

export EstimationMethod, FOCEIMethod, SAEMMethod, LaplacianMethod
export EstimationConfig, IndividualEstimate, EstimationResult
export OmegaStructure, DiagonalOmega, BlockOmega

# ------------------------------------------------------------------
# Estimation Methods
# ------------------------------------------------------------------

"""
Abstract type for estimation methods.
"""
abstract type EstimationMethod end

"""
First-Order Conditional Estimation with Interaction (FOCE-I).

This is the most commonly used method in NONMEM. It linearizes the model
around the conditional mode of the random effects (eta).

Fields:
- max_inner_iter: Maximum iterations for finding eta mode (default: 50)
- inner_tol: Convergence tolerance for inner optimization (default: 1e-6)
- centered: Use centered differences for FO approximation (default: false)
"""
struct FOCEIMethod <: EstimationMethod
    max_inner_iter::Int
    inner_tol::Float64
    centered::Bool

    function FOCEIMethod(;
        max_inner_iter::Int=50,
        inner_tol::Float64=1e-6,
        centered::Bool=false
    )
        new(max_inner_iter, inner_tol, centered)
    end
end

"""
Stochastic Approximation Expectation Maximization (SAEM).

SAEM uses MCMC sampling to approximate the E-step and stochastic
approximation for the M-step. Often more robust than FOCE-I for
complex models.

Fields:
- n_burn: Number of burn-in iterations (default: 300)
- n_iter: Number of main iterations after burn-in (default: 200)
- n_chains: Number of MCMC chains per subject (default: 3)
- step_size_schedule: Step size decay schedule (:harmonic or :constant)
"""
struct SAEMMethod <: EstimationMethod
    n_burn::Int
    n_iter::Int
    n_chains::Int
    step_size_schedule::Symbol

    function SAEMMethod(;
        n_burn::Int=300,
        n_iter::Int=200,
        n_chains::Int=3,
        step_size_schedule::Symbol=:harmonic
    )
        @assert step_size_schedule in (:harmonic, :constant) "Invalid step size schedule"
        new(n_burn, n_iter, n_chains, step_size_schedule)
    end
end

"""
Laplacian Estimation Method.

The simplest approximation method. Uses the mode of the posterior
for random effects and the Laplace approximation for the marginal
likelihood. Less accurate than FOCE-I but computationally faster.
"""
struct LaplacianMethod <: EstimationMethod
    max_inner_iter::Int
    inner_tol::Float64

    function LaplacianMethod(;
        max_inner_iter::Int=50,
        inner_tol::Float64=1e-6
    )
        new(max_inner_iter, inner_tol)
    end
end

# ------------------------------------------------------------------
# Omega Structure
# ------------------------------------------------------------------

"""
Abstract type for omega matrix structure.
"""
abstract type OmegaStructure end

"""
Diagonal omega matrix - no correlations between random effects.
"""
struct DiagonalOmega <: OmegaStructure end

"""
Block diagonal omega matrix - allows correlations within blocks.

Fields:
- block_sizes: Vector of block sizes (e.g., [2, 1] for a 2x2 block and 1x1 block)
"""
struct BlockOmega <: OmegaStructure
    block_sizes::Vector{Int}
end

"""
Full omega matrix - all random effects correlated.
"""
struct FullOmega <: OmegaStructure end

export FullOmega

# ------------------------------------------------------------------
# Estimation Configuration
# ------------------------------------------------------------------

"""
Configuration for parameter estimation.

Fields:
- method: Estimation method (FOCEIMethod, SAEMMethod, or LaplacianMethod)
- theta_init: Initial values for fixed effects (population parameters)
- theta_lower: Lower bounds for theta
- theta_upper: Upper bounds for theta
- theta_names: Names of theta parameters
- omega_init: Initial omega matrix (variance of random effects)
- omega_structure: Structure of omega matrix
- omega_names: Names of eta parameters
- sigma_init: Initial residual error specification
- max_iter: Maximum outer iterations
- tol: Convergence tolerance
- compute_se: Whether to compute standard errors
- compute_ci: Whether to compute confidence intervals
- ci_level: Confidence level (default: 0.95)
- verbose: Print progress during estimation
- seed: Random seed for reproducibility
"""
struct EstimationConfig{M<:EstimationMethod}
    method::M
    theta_init::Vector{Float64}
    theta_lower::Vector{Float64}
    theta_upper::Vector{Float64}
    theta_names::Vector{Symbol}
    omega_init::Matrix{Float64}
    omega_structure::OmegaStructure
    omega_names::Vector{Symbol}
    sigma_init::ResidualErrorSpec
    max_iter::Int
    tol::Float64
    compute_se::Bool
    compute_ci::Bool
    ci_level::Float64
    verbose::Bool
    seed::UInt64

    function EstimationConfig(
        method::M;
        theta_init::Vector{Float64},
        theta_lower::Vector{Float64}=fill(-Inf, length(theta_init)),
        theta_upper::Vector{Float64}=fill(Inf, length(theta_init)),
        theta_names::Vector{Symbol}=Symbol[],
        omega_init::Matrix{Float64},
        omega_structure::OmegaStructure=DiagonalOmega(),
        omega_names::Vector{Symbol}=Symbol[],
        sigma_init::ResidualErrorSpec,
        max_iter::Int=500,
        tol::Float64=1e-4,
        compute_se::Bool=true,
        compute_ci::Bool=false,
        ci_level::Float64=0.95,
        verbose::Bool=false,
        seed::UInt64=UInt64(12345)
    ) where {M<:EstimationMethod}
        n_theta = length(theta_init)
        n_eta = size(omega_init, 1)

        # Validate dimensions
        @assert length(theta_lower) == n_theta "theta_lower length mismatch"
        @assert length(theta_upper) == n_theta "theta_upper length mismatch"
        @assert size(omega_init, 1) == size(omega_init, 2) "omega_init must be square"
        @assert issymmetric(omega_init) "omega_init must be symmetric"
        @assert all(theta_lower .<= theta_init .<= theta_upper) "theta_init out of bounds"
        @assert 0.0 < ci_level < 1.0 "ci_level must be in (0, 1)"

        # Default names if not provided
        if isempty(theta_names)
            theta_names = [Symbol("theta_$i") for i in 1:n_theta]
        end
        if isempty(omega_names)
            omega_names = [Symbol("eta_$i") for i in 1:n_eta]
        end

        new{M}(
            method, theta_init, theta_lower, theta_upper, theta_names,
            omega_init, omega_structure, omega_names, sigma_init,
            max_iter, tol, compute_se, compute_ci, ci_level, verbose, seed
        )
    end
end

# ------------------------------------------------------------------
# Individual Estimates
# ------------------------------------------------------------------

"""
Estimated quantities for a single individual.

Fields:
- subject_id: Subject identifier
- eta: Estimated random effects (empirical Bayes estimates)
- eta_se: Standard errors of eta (if computed)
- ipred: Individual predictions at observation times
- pred: Population predictions at observation times
- cwres: Conditional weighted residuals
- iwres: Individual weighted residuals
- wres: Weighted residuals (population)
- ofv_contribution: Subject's contribution to objective function
"""
struct IndividualEstimate
    subject_id::String
    eta::Vector{Float64}
    eta_se::Union{Nothing,Vector{Float64}}
    ipred::Vector{Float64}
    pred::Vector{Float64}
    cwres::Vector{Float64}
    iwres::Vector{Float64}
    wres::Vector{Float64}
    ofv_contribution::Float64
end

# ------------------------------------------------------------------
# Estimation Result
# ------------------------------------------------------------------

"""
Result of parameter estimation.

Fields:
- config: EstimationConfig used
- theta: Estimated fixed effects
- theta_se: Standard errors of theta
- theta_rse: Relative standard errors (CV%)
- theta_ci_lower: Lower confidence bounds
- theta_ci_upper: Upper confidence bounds
- omega: Estimated omega matrix
- omega_se: Standard errors of omega elements
- omega_corr: Correlation matrix derived from omega
- sigma: Estimated residual error specification
- sigma_se: Standard errors of sigma parameters
- ofv: Objective function value (-2LL)
- aic: Akaike Information Criterion
- bic: Bayesian Information Criterion
- individuals: Per-subject estimates
- convergence: Did optimization converge?
- n_iterations: Number of iterations used
- gradient_norm: Final gradient norm
- condition_number: Condition number of Hessian
- eigenvalue_ratio: Ratio of smallest to largest eigenvalue
- covariance_step_successful: Did SE computation succeed?
- messages: Any warnings or messages
- runtime_seconds: Total runtime
"""
struct EstimationResult{M<:EstimationMethod}
    config::EstimationConfig{M}

    # Fixed effects
    theta::Vector{Float64}
    theta_se::Union{Nothing,Vector{Float64}}
    theta_rse::Union{Nothing,Vector{Float64}}
    theta_ci_lower::Union{Nothing,Vector{Float64}}
    theta_ci_upper::Union{Nothing,Vector{Float64}}

    # Random effects variance
    omega::Matrix{Float64}
    omega_se::Union{Nothing,Matrix{Float64}}
    omega_corr::Matrix{Float64}

    # Residual error
    sigma::ResidualErrorSpec
    sigma_se::Union{Nothing,Vector{Float64}}

    # Model fit statistics
    ofv::Float64
    aic::Float64
    bic::Float64

    # Individual estimates
    individuals::Vector{IndividualEstimate}

    # Convergence information
    convergence::Bool
    n_iterations::Int
    gradient_norm::Float64
    condition_number::Float64
    eigenvalue_ratio::Float64
    covariance_step_successful::Bool

    # Messages and runtime
    messages::Vector{String}
    runtime_seconds::Float64
end

# ------------------------------------------------------------------
# Utility Functions for Results
# ------------------------------------------------------------------

"""
Get number of estimated parameters (for AIC/BIC calculation).
"""
function n_parameters(result::EstimationResult)::Int
    n_theta = length(result.theta)
    n_omega = _count_omega_params(result.config.omega_structure, size(result.omega, 1))
    n_sigma = _count_sigma_params(result.sigma)
    return n_theta + n_omega + n_sigma
end

function _count_omega_params(::DiagonalOmega, n::Int)::Int
    return n  # Just the diagonal elements
end

function _count_omega_params(structure::BlockOmega, n::Int)::Int
    total = 0
    for block_size in structure.block_sizes
        total += block_size * (block_size + 1) ÷ 2
    end
    return total
end

function _count_omega_params(::FullOmega, n::Int)::Int
    return n * (n + 1) ÷ 2  # Lower triangle including diagonal
end

function _count_sigma_params(spec::ResidualErrorSpec{AdditiveError})::Int
    return 1  # sigma
end

function _count_sigma_params(spec::ResidualErrorSpec{ProportionalError})::Int
    return 1  # sigma
end

function _count_sigma_params(spec::ResidualErrorSpec{CombinedError})::Int
    return 2  # sigma_add, sigma_prop
end

function _count_sigma_params(spec::ResidualErrorSpec{ExponentialError})::Int
    return 1  # sigma
end

"""
Get number of observations across all subjects.
"""
function n_observations(result::EstimationResult)::Int
    return sum(length(ind.ipred) for ind in result.individuals)
end

"""
Get number of subjects.
"""
function n_subjects(result::EstimationResult)::Int
    return length(result.individuals)
end

export n_parameters, n_observations, n_subjects

# ------------------------------------------------------------------
# Pretty Printing
# ------------------------------------------------------------------

function Base.show(io::IO, result::EstimationResult)
    println(io, "EstimationResult")
    println(io, "  Method: ", typeof(result.config.method))
    println(io, "  Convergence: ", result.convergence ? "Yes" : "No")
    println(io, "  OFV: ", round(result.ofv, digits=3))
    println(io, "  AIC: ", round(result.aic, digits=3))
    println(io, "  BIC: ", round(result.bic, digits=3))
    println(io, "  Subjects: ", n_subjects(result))
    println(io, "  Observations: ", n_observations(result))
    println(io, "  Iterations: ", result.n_iterations)
    if result.theta_se !== nothing
        println(io, "\n  Fixed Effects (theta):")
        for (i, (name, val, se)) in enumerate(zip(
            result.config.theta_names, result.theta, result.theta_se
        ))
            rse = result.theta_rse !== nothing ? result.theta_rse[i] : NaN
            println(io, "    $name: $(round(val, digits=4)) (SE=$(round(se, digits=4)), RSE=$(round(rse, digits=1))%)")
        end
    end
end
