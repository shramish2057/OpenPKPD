export apply_covariates

function _require_has_covariate(covs::IndividualCovariates, name::Symbol)
    if !haskey(covs.values, name)
        error("Missing covariate $(name) in IndividualCovariates")
    end
    return covs.values[name]
end

function _apply_effect(theta::Float64, cov::Float64, eff::CovariateEffect{LinearCovariate})
    return theta * (1.0 + eff.beta * (cov - eff.ref))
end

function _apply_effect(theta::Float64, cov::Float64, eff::CovariateEffect{PowerCovariate})
    return theta * (cov / eff.ref) ^ eff.beta
end

function _apply_effect(theta::Float64, cov::Float64, eff::CovariateEffect{ExpCovariate})
    return theta * exp(eff.beta * (cov - eff.ref))
end

"""
Apply covariate effects to a typed parameter struct.

Order:
- effects are applied sequentially in the order provided
- multiple effects can target the same parameter (explicit order)
"""
function apply_covariates(params, cov_model::CovariateModel, covs::IndividualCovariates)
    T = typeof(params)
    fn = fieldnames(T)

    # Start from base param values
    vals = Dict{Symbol,Float64}()
    for f in fn
        vals[f] = Float64(getfield(params, f))
    end

    for eff in cov_model.effects
        if !(eff.param in fn)
            error("Covariate effect targets unknown param $(eff.param) for $(T)")
        end

        cov = _require_has_covariate(covs, eff.covariate)
        θ = vals[eff.param]

        newθ = _apply_effect(θ, cov, eff)

        if !isfinite(newθ)
            error("Non-finite covariate transformed value for $(eff.param)")
        end

        vals[eff.param] = newθ
    end

    # Rebuild typed params
    return T((vals[f] for f in fn)...), vals
end
