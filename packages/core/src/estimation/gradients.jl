# Gradient Computation Utilities for Parameter Estimation
# Provides finite difference and automatic differentiation support

using ForwardDiff
using LinearAlgebra

export compute_gradient, compute_hessian, gradient_fd, hessian_fd

"""
Compute gradient of a function using finite differences.

Arguments:
- f: Objective function f(x) -> Float64
- x: Point at which to compute gradient
- h: Step size for finite differences (default: 1e-6)

Returns:
- Vector of partial derivatives
"""
function gradient_fd(
    f::Function,
    x::Vector{Float64};
    h::Float64=1e-6
)::Vector{Float64}
    n = length(x)
    grad = zeros(n)

    for i in 1:n
        x_plus = copy(x)
        x_minus = copy(x)
        x_plus[i] += h
        x_minus[i] -= h

        # Central difference
        grad[i] = (f(x_plus) - f(x_minus)) / (2 * h)
    end

    return grad
end

"""
Compute Hessian of a function using finite differences.

Arguments:
- f: Objective function f(x) -> Float64
- x: Point at which to compute Hessian
- h: Step size for finite differences (default: 1e-5)

Returns:
- Symmetric Hessian matrix
"""
function hessian_fd(
    f::Function,
    x::Vector{Float64};
    h::Float64=1e-5
)::Matrix{Float64}
    n = length(x)
    hess = zeros(n, n)
    f_center = f(x)

    for i in 1:n
        for j in i:n
            if i == j
                # Diagonal: second derivative
                x_plus = copy(x)
                x_minus = copy(x)
                x_plus[i] += h
                x_minus[i] -= h

                hess[i, i] = (f(x_plus) - 2 * f_center + f(x_minus)) / h^2
            else
                # Off-diagonal: mixed partial
                x_pp = copy(x)
                x_pm = copy(x)
                x_mp = copy(x)
                x_mm = copy(x)

                x_pp[i] += h; x_pp[j] += h
                x_pm[i] += h; x_pm[j] -= h
                x_mp[i] -= h; x_mp[j] += h
                x_mm[i] -= h; x_mm[j] -= h

                hess[i, j] = (f(x_pp) - f(x_pm) - f(x_mp) + f(x_mm)) / (4 * h^2)
                hess[j, i] = hess[i, j]
            end
        end
    end

    return hess
end

"""
Compute gradient using ForwardDiff automatic differentiation.

Arguments:
- f: Objective function f(x) -> Float64
- x: Point at which to compute gradient

Returns:
- Vector of partial derivatives
"""
function compute_gradient(
    f::Function,
    x::Vector{Float64}
)::Vector{Float64}
    return ForwardDiff.gradient(f, x)
end

"""
Compute Hessian using ForwardDiff automatic differentiation.

Arguments:
- f: Objective function f(x) -> Float64
- x: Point at which to compute Hessian

Returns:
- Hessian matrix
"""
function compute_hessian(
    f::Function,
    x::Vector{Float64}
)::Matrix{Float64}
    return ForwardDiff.hessian(f, x)
end

"""
Check gradient accuracy by comparing AD and FD.
Useful for debugging.

Returns maximum relative difference between AD and FD gradients.
"""
function check_gradient(
    f::Function,
    x::Vector{Float64};
    h::Float64=1e-6,
    verbose::Bool=false
)::Float64
    grad_ad = compute_gradient(f, x)
    grad_fd = gradient_fd(f, x; h=h)

    if verbose
        println("AD gradient: ", grad_ad)
        println("FD gradient: ", grad_fd)
        println("Difference:  ", grad_ad .- grad_fd)
    end

    # Relative difference
    rel_diff = abs.(grad_ad .- grad_fd) ./ max.(abs.(grad_ad), 1.0)
    return maximum(rel_diff)
end

export check_gradient

"""
Compute Jacobian of a vector-valued function.

Arguments:
- f: Function f(x) -> Vector{Float64}
- x: Point at which to compute Jacobian

Returns:
- Jacobian matrix (m x n where m = length(f(x)), n = length(x))
"""
function compute_jacobian(
    f::Function,
    x::Vector{Float64}
)::Matrix{Float64}
    return ForwardDiff.jacobian(f, x)
end

export compute_jacobian

"""
Compute numerical condition number of a matrix.
"""
function condition_number(A::Matrix{Float64})::Float64
    eigenvalues = eigvals(A)
    real_eigenvalues = real.(eigenvalues)

    max_eig = maximum(abs.(real_eigenvalues))
    min_eig = minimum(abs.(real_eigenvalues))

    if min_eig < 1e-12
        return Inf
    end

    return max_eig / min_eig
end

export condition_number

"""
Check if matrix is positive definite.
"""
function is_positive_definite(A::Matrix{Float64})::Bool
    try
        cholesky(Symmetric(A))
        return true
    catch
        return false
    end
end

export is_positive_definite
