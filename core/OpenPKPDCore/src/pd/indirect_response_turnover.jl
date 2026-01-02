export validate, inhibition

function validate(spec::PDSpec{IndirectResponseTurnover,IndirectResponseTurnoverParams})
    p = spec.params

    _require_positive("Kin", p.Kin)
    _require_positive("Kout", p.Kout)
    _require_positive("R0", p.R0)

    if p.Imax < 0.0 || p.Imax > 1.0
        error("Imax must be within [0, 1], got $(p.Imax)")
    end

    _require_positive("IC50", p.IC50)

    return nothing
end

"""
Inhibition term I(C) in [0, Imax].
"""
function inhibition(C::Float64, Imax::Float64, IC50::Float64)
    if C <= 0.0
        return 0.0
    end
    return (Imax * C) / (IC50 + C)
end
