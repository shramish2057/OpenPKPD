# TMDD Simulation/Solve Functions
#
# High-level simulation interface for TMDD models.
# Integrates with OpenPKPD's solver and callback infrastructure.

using DifferentialEquations
using SciMLBase

export solve_tmdd, TMDDSimResult

# ============================================================================
# TMDD Simulation Result
# ============================================================================

"""
Result container for TMDD simulation.

Contains:
- t: Time points
- states: All state variables by name
- observations: Derived quantities (concentrations, occupancy, etc.)
- metadata: Additional information about the simulation
"""
struct TMDDSimResult
    t::Vector{Float64}
    states::Dict{Symbol,Vector{Float64}}
    observations::Dict{Symbol,Vector{Float64}}
    metadata::Dict{String,Any}
end

# ============================================================================
# Dosing Callbacks
# ============================================================================

"""
Create callbacks for TMDD dosing events.
"""
function _create_tmdd_callbacks(spec::TMDDSpec, dose_cpt::Int)
    doses = spec.doses
    if isempty(doses)
        return nothing
    end

    callbacks = VectorContinuousCallback[]

    for dose in doses
        if is_bolus(dose)
            # Bolus: instantaneous dose at dose.time
            condition = (u, t, integrator) -> t - dose.time
            affect! = function(integrator)
                integrator.u[dose_cpt] += dose.amount
            end
            push!(callbacks, ContinuousCallback(condition, affect!))
        else
            # Infusion: rate = dose.amount / dose.duration over [dose.time, dose.time + dose.duration]
            rate = dose.amount / dose.duration
            t_start = dose.time
            t_end = dose.time + dose.duration

            # We'll handle infusions via a piecewise ODE modification
            # For simplicity, use DiscreteCallback to add the dose gradually
            # This is a simplified approach; full infusion handling would require
            # modifying the ODE during the infusion period

            # Start infusion
            condition_start = (u, t, integrator) -> t - t_start
            affect_start! = function(integrator)
                # Store infusion rate in parameters (if supported)
                # For now, we add half at start and half at end as approximation
                # TODO: Implement proper zero-order infusion
            end

            # For IV bolus approximation of infusion
            condition_end = (u, t, integrator) -> t - t_end
            affect_end! = function(integrator)
                integrator.u[dose_cpt] += dose.amount
            end
            push!(callbacks, ContinuousCallback(condition_end, affect_end!))
        end
    end

    return CallbackSet(callbacks...)
end

"""
Create discrete callbacks for bolus dosing (more efficient than continuous).
"""
function _create_bolus_callbacks(doses::Vector{DoseEvent}, dose_cpt::Int)
    bolus_doses = filter(is_bolus, doses)
    if isempty(bolus_doses)
        return nothing
    end

    dose_times = [d.time for d in bolus_doses]
    dose_amounts = [d.amount for d in bolus_doses]

    affect! = function(integrator)
        idx = findfirst(t -> isapprox(t, integrator.t, atol=1e-10), dose_times)
        if idx !== nothing
            integrator.u[dose_cpt] += dose_amounts[idx]
        end
    end

    return PresetTimeCallback(dose_times, affect!)
end

# ============================================================================
# Main Solve Function
# ============================================================================

"""
    solve_tmdd(spec::TMDDSpec, grid::SimGrid, solver::SolverSpec) -> TMDDSimResult

Simulate a TMDD model.

# Arguments
- `spec::TMDDSpec`: TMDD model specification with parameters and dosing
- `grid::SimGrid`: Simulation time grid (t0, t1, saveat)
- `solver::SolverSpec`: ODE solver configuration

# Returns
- `TMDDSimResult`: Simulation results with states and derived observations

# Example
```julia
spec = TMDDSpec(
    TMDD2CptQSS(),
    "mAb PK",
    TMDD2CptQSSParams(
        kel=0.01, V1=3.0, V2=4.0, Q=0.5,
        KSS=0.05, ksyn=0.1, kdeg=0.1, kint=0.05, Rtot0=1.0
    ),
    [DoseEvent(0.0, 200.0)]
)

grid = SimGrid(0.0, 672.0, collect(0.0:1.0:672.0))  # 28 days
solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

result = solve_tmdd(spec, grid, solver)
```
"""
function solve_tmdd(spec::TMDDSpec{K,P}, grid::SimGrid, solver::SolverSpec) where {K<:TMDDModelKind,P}
    # Validate
    validate_tmdd(spec)

    # Get model components
    ode_fn! = get_ode_function(spec.kind)
    u0 = get_initial_state(spec)
    dose_cpt = dosing_compartment(spec.kind)
    names = state_names(spec.kind)

    # Apply first dose at t=0 if applicable
    first_dose_at_zero = !isempty(spec.doses) && spec.doses[1].time == 0.0 && is_bolus(spec.doses[1])
    if first_dose_at_zero
        u0[dose_cpt] += spec.doses[1].amount
    end

    # Create callbacks for remaining doses
    remaining_doses = first_dose_at_zero ? spec.doses[2:end] : spec.doses
    callbacks = _create_bolus_callbacks(remaining_doses, dose_cpt)

    # Get solver algorithm
    alg = _get_algorithm(solver.alg)

    # Setup ODE problem
    tspan = (grid.t0, grid.t1)
    prob = ODEProblem(ode_fn!, u0, tspan, spec.params)

    # Solve
    sol = solve(
        prob,
        alg;
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=callbacks
    )

    # Extract results
    t = sol.t
    n_times = length(t)
    n_states_val = n_states(spec.kind)

    # Build states dictionary
    states = Dict{Symbol,Vector{Float64}}()
    for (i, name) in enumerate(names)
        states[name] = [sol.u[j][i] for j in 1:n_times]
    end

    # Compute derived observations
    observations = _compute_tmdd_observations(spec, states, t)

    # Metadata
    metadata = Dict{String,Any}(
        "model_type" => string(K),
        "model_name" => spec.name,
        "n_doses" => length(spec.doses),
        "solver" => string(solver.alg),
        "retcode" => string(sol.retcode)
    )

    return TMDDSimResult(t, states, observations, metadata)
end

"""
Get DifferentialEquations algorithm from symbol.
"""
function _get_algorithm(alg::Symbol)
    if alg == :Tsit5
        return Tsit5()
    elseif alg == :Rodas4
        return Rodas4()
    elseif alg == :Rodas5
        return Rodas5()
    elseif alg == :CVODE_BDF
        # Requires Sundials
        return Rodas4()  # Fallback
    elseif alg == :Rosenbrock23
        return Rosenbrock23()
    elseif alg == :TRBDF2
        return TRBDF2()
    else
        return Tsit5()  # Default
    end
end

# ============================================================================
# Observation Computation
# ============================================================================

"""
Compute derived observations from TMDD states.
"""
function _compute_tmdd_observations(spec::TMDDSpec{TMDDFull,TMDDFullParams}, states, t)
    p = spec.params
    V = p.V

    L = states[:L]
    R = states[:R]
    P = states[:P]

    observations = Dict{Symbol,Vector{Float64}}()

    # Concentrations
    observations[:conc] = L ./ V  # Free drug concentration
    observations[:conc_free] = L ./ V
    observations[:conc_bound] = P ./ V  # Complex concentration
    observations[:conc_total] = (L .+ P) ./ V  # Total drug concentration

    # Target
    observations[:R_free] = R ./ V  # Free target concentration
    observations[:R_total] = (R .+ P) ./ V  # Total target concentration

    # Target occupancy
    Rtot = R .+ P
    observations[:target_occupancy] = [Rtot[i] > 0 ? P[i] / Rtot[i] : 0.0 for i in 1:length(t)]

    return observations
end

function _compute_tmdd_observations(spec::TMDDSpec{TMDDQSS,TMDDQSSParams}, states, t)
    p = spec.params
    V = p.V
    KSS = p.KSS

    L = states[:L]
    Rtot = states[:Rtot]

    observations = Dict{Symbol,Vector{Float64}}()

    # QSS complex calculation
    KSS_V = KSS * V
    RL = [L[i] * Rtot[i] / (KSS_V + L[i]) for i in 1:length(t)]

    # Concentrations
    observations[:conc] = L ./ V
    observations[:conc_free] = L ./ V
    observations[:conc_bound] = RL ./ V
    observations[:conc_total] = (L .+ RL) ./ V

    # Target
    R_free = Rtot .- RL
    observations[:R_free] = R_free ./ V
    observations[:R_total] = Rtot ./ V

    # Target occupancy
    observations[:target_occupancy] = [Rtot[i] > 0 ? RL[i] / Rtot[i] : 0.0 for i in 1:length(t)]

    return observations
end

function _compute_tmdd_observations(spec::TMDDSpec{TMDDQE,TMDDQEParams}, states, t)
    p = spec.params
    V = p.V
    KD = p.KD

    L = states[:L]
    Rtot = states[:Rtot]

    observations = Dict{Symbol,Vector{Float64}}()

    # QE complex calculation
    KD_V = KD * V
    RL = [L[i] * Rtot[i] / (KD_V + L[i]) for i in 1:length(t)]

    # Concentrations
    observations[:conc] = L ./ V
    observations[:conc_free] = L ./ V
    observations[:conc_bound] = RL ./ V
    observations[:conc_total] = (L .+ RL) ./ V

    # Target
    R_free = Rtot .- RL
    observations[:R_free] = R_free ./ V
    observations[:R_total] = Rtot ./ V

    # Target occupancy
    observations[:target_occupancy] = [Rtot[i] > 0 ? RL[i] / Rtot[i] : 0.0 for i in 1:length(t)]

    return observations
end

function _compute_tmdd_observations(spec::TMDDSpec{TMDDMM,TMDDMMParams}, states, t)
    p = spec.params
    V = p.V

    L = states[:L]

    observations = Dict{Symbol,Vector{Float64}}()
    observations[:conc] = L ./ V
    observations[:conc_free] = L ./ V
    observations[:conc_total] = L ./ V

    return observations
end

function _compute_tmdd_observations(spec::TMDDSpec{TMDDRapidBinding,TMDDRapidBindingParams}, states, t)
    p = spec.params
    V = p.V
    KD = p.KD

    Ltot = states[:Ltot]
    Rtot = states[:Rtot]

    observations = Dict{Symbol,Vector{Float64}}()

    # Calculate free drug from quadratic
    KD_V = KD * V
    L_free = Vector{Float64}(undef, length(t))
    for i in 1:length(t)
        a = Ltot[i] - Rtot[i] - KD_V
        discriminant = a^2 + 4 * KD_V * Ltot[i]
        L_free[i] = max(0.0, 0.5 * (a + sqrt(max(0.0, discriminant))))
    end

    P = Ltot .- L_free

    observations[:conc] = L_free ./ V
    observations[:conc_free] = L_free ./ V
    observations[:conc_bound] = P ./ V
    observations[:conc_total] = Ltot ./ V

    # Target
    R_free = Rtot .- P
    observations[:R_free] = R_free ./ V
    observations[:R_total] = Rtot ./ V

    # Target occupancy
    observations[:target_occupancy] = [Rtot[i] > 0 ? P[i] / Rtot[i] : 0.0 for i in 1:length(t)]

    return observations
end

function _compute_tmdd_observations(spec::TMDDSpec{TMDDIrreversible,TMDDIrreversibleParams}, states, t)
    p = spec.params
    V = p.V

    L = states[:L]
    R = states[:R]

    observations = Dict{Symbol,Vector{Float64}}()

    observations[:conc] = L ./ V
    observations[:conc_free] = L ./ V
    observations[:R_free] = R ./ V

    return observations
end

function _compute_tmdd_observations(spec::TMDDSpec{TMDD2CptFull,TMDD2CptFullParams}, states, t)
    p = spec.params
    V1 = p.V1
    V2 = p.V2

    L = states[:L]
    Lp = states[:Lp]
    R = states[:R]
    P = states[:P]

    observations = Dict{Symbol,Vector{Float64}}()

    # Central concentrations
    observations[:conc] = L ./ V1
    observations[:conc_free] = L ./ V1
    observations[:conc_bound] = P ./ V1
    observations[:conc_total] = (L .+ P) ./ V1

    # Peripheral
    observations[:conc_peripheral] = Lp ./ V2

    # Target
    observations[:R_free] = R ./ V1
    observations[:R_total] = (R .+ P) ./ V1

    # Target occupancy
    Rtot = R .+ P
    observations[:target_occupancy] = [Rtot[i] > 0 ? P[i] / Rtot[i] : 0.0 for i in 1:length(t)]

    return observations
end

function _compute_tmdd_observations(spec::TMDDSpec{TMDD2CptQSS,TMDD2CptQSSParams}, states, t)
    p = spec.params
    V1 = p.V1
    V2 = p.V2
    KSS = p.KSS

    L = states[:L]
    Lp = states[:Lp]
    Rtot = states[:Rtot]

    observations = Dict{Symbol,Vector{Float64}}()

    # QSS complex
    KSS_V1 = KSS * V1
    RL = [L[i] * Rtot[i] / (KSS_V1 + L[i]) for i in 1:length(t)]

    # Central concentrations
    observations[:conc] = L ./ V1
    observations[:conc_free] = L ./ V1
    observations[:conc_bound] = RL ./ V1
    observations[:conc_total] = (L .+ RL) ./ V1

    # Peripheral
    observations[:conc_peripheral] = Lp ./ V2

    # Target
    R_free = Rtot .- RL
    observations[:R_free] = R_free ./ V1
    observations[:R_total] = Rtot ./ V1

    # Target occupancy
    observations[:target_occupancy] = [Rtot[i] > 0 ? RL[i] / Rtot[i] : 0.0 for i in 1:length(t)]

    return observations
end

function _compute_tmdd_observations(spec::TMDDSpec{TMDD2CptCL,TMDD2CptCLParams}, states, t)
    p = spec.params
    V1 = p.V1
    V2 = p.V2
    Kss = p.Kss

    L = states[:L]
    Lp = states[:Lp]
    Rtot = states[:Rtot]

    observations = Dict{Symbol,Vector{Float64}}()

    # QSS complex
    Kss_V1 = Kss * V1
    RL = [L[i] * Rtot[i] / (Kss_V1 + L[i]) for i in 1:length(t)]

    # Central concentrations
    observations[:conc] = L ./ V1
    observations[:conc_free] = L ./ V1
    observations[:conc_bound] = RL ./ V1
    observations[:conc_total] = (L .+ RL) ./ V1

    # Peripheral
    observations[:conc_peripheral] = Lp ./ V2

    # Target
    R_free = Rtot .- RL
    observations[:R_free] = R_free ./ V1
    observations[:R_total] = Rtot ./ V1

    # Target occupancy
    observations[:target_occupancy] = [Rtot[i] > 0 ? RL[i] / Rtot[i] : 0.0 for i in 1:length(t)]

    return observations
end

function _compute_tmdd_observations(spec::TMDDSpec{TMDDSolubleTarget,TMDDSolubleTargetParams}, states, t)
    p = spec.params
    V = p.V

    L = states[:L]
    R = states[:R]
    P = states[:P]

    observations = Dict{Symbol,Vector{Float64}}()

    # Concentrations
    observations[:conc] = L ./ V
    observations[:conc_free] = L ./ V
    observations[:conc_bound] = P ./ V
    observations[:conc_total] = (L .+ P) ./ V

    # Target
    observations[:R_free] = R ./ V
    observations[:R_total] = (R .+ P) ./ V

    # Target occupancy
    Rtot = R .+ P
    observations[:target_occupancy] = [Rtot[i] > 0 ? P[i] / Rtot[i] : 0.0 for i in 1:length(t)]

    return observations
end

function _compute_tmdd_observations(spec::TMDDSpec{TMDDInternalization,TMDDInternalizationParams}, states, t)
    p = spec.params
    V = p.V

    L = states[:L]
    R = states[:R]
    P = states[:P]
    Re = states[:Re]
    Pe = states[:Pe]

    observations = Dict{Symbol,Vector{Float64}}()

    # Concentrations
    observations[:conc] = L ./ V
    observations[:conc_free] = L ./ V
    observations[:conc_bound_surface] = P ./ V
    observations[:conc_total] = (L .+ P) ./ V

    # Surface target
    observations[:R_surface_free] = R ./ V
    observations[:R_surface_total] = (R .+ P) ./ V

    # Endosomal
    observations[:R_endo_free] = Re ./ V
    observations[:R_endo_bound] = Pe ./ V

    # Total target (all compartments)
    R_total = R .+ P .+ Re .+ Pe
    observations[:R_total] = R_total ./ V

    # Surface target occupancy
    R_surface_total = R .+ P
    observations[:target_occupancy] = [R_surface_total[i] > 0 ? P[i] / R_surface_total[i] : 0.0 for i in 1:length(t)]

    return observations
end
