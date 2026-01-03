export PerturbationKind, RelativePerturbation, AbsolutePerturbation, Perturbation
export PerturbationPlan, SensitivityMetric

abstract type PerturbationKind end

"""
Relative perturbation of a parameter: new = base * (1 + delta)
delta is a fraction, e.g. 0.1 means +10 percent
"""
struct RelativePerturbation <: PerturbationKind end

"""
Absolute perturbation of a parameter: new = base + delta
"""
struct AbsolutePerturbation <: PerturbationKind end

struct Perturbation{K<:PerturbationKind}
    kind::K
    param::Symbol
    delta::Float64
end

"""
A named set of perturbations.
"""
struct PerturbationPlan
    name::String
    perturbations::Vector{Perturbation}
end

"""
Simple sensitivity metrics for time series comparison.
"""
struct SensitivityMetric
    max_abs_delta::Float64
    max_rel_delta::Float64
    l2_norm_delta::Float64
end
