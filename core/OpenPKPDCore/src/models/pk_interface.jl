export pk_validate, pk_param_tuple, pk_state_symbols, pk_u0, pk_ode!, pk_conc
export pk_dose_target_index

"""
Internal PK interface.

Each supported PK model kind must implement these methods:
- pk_validate(spec)
- pk_param_tuple(spec)
- pk_state_symbols(kind)
- pk_u0(spec, grid)
- pk_ode!(du, u, p, t, kind)
- pk_conc(u, p, kind)
- pk_dose_target_index(kind)
"""

# -------------------------
# Generic fallbacks
# -------------------------

function pk_validate(spec::ModelSpec)
    validate(spec)
end

# -------------------------
# OneCompIVBolus
# -------------------------

function pk_param_tuple(spec::ModelSpec{OneCompIVBolus,OneCompIVBolusParams})
    return (CL=spec.params.CL, V=spec.params.V)
end

pk_state_symbols(::OneCompIVBolus) = [:A_central]

pk_dose_target_index(::OneCompIVBolus) = 1

function pk_u0(spec::ModelSpec{OneCompIVBolus,OneCompIVBolusParams}, grid::SimGrid)
    A0 = 0.0
    for d in spec.doses
        if d.time == grid.t0
            A0 += d.amount
        end
    end
    return [A0]
end

function pk_ode!(du, u, p, t, ::OneCompIVBolus)
    A = u[1]
    du[1] = -(p.CL / p.V) * A
    return nothing
end

function pk_conc(u, p, ::OneCompIVBolus)
    return u[1] / p.V
end

# -------------------------
# OneCompOralFirstOrder
# -------------------------

function pk_param_tuple(spec::ModelSpec{OneCompOralFirstOrder,OneCompOralFirstOrderParams})
    return (Ka=spec.params.Ka, CL=spec.params.CL, V=spec.params.V)
end

pk_state_symbols(::OneCompOralFirstOrder) = [:A_gut, :A_central]

pk_dose_target_index(::OneCompOralFirstOrder) = 1

function pk_u0(
    spec::ModelSpec{OneCompOralFirstOrder,OneCompOralFirstOrderParams}, grid::SimGrid
)
    Agut0 = 0.0
    for d in spec.doses
        if d.time == grid.t0
            Agut0 += d.amount
        end
    end
    return [Agut0, 0.0]
end

function pk_ode!(du, u, p, t, ::OneCompOralFirstOrder)
    Agut = u[1]
    Acent = u[2]

    du[1] = -p.Ka * Agut
    du[2] = p.Ka * Agut - (p.CL / p.V) * Acent

    return nothing
end

function pk_conc(u, p, ::OneCompOralFirstOrder)
    return u[2] / p.V
end
