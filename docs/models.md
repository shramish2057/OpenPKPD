# PK/PD Models Reference

OpenPKPD provides validated pharmacokinetic (PK) and pharmacodynamic (PD) models for drug concentration and effect simulation.

## Pharmacokinetic Models

### One-Compartment IV Bolus

The simplest PK model representing instantaneous drug administration into a single well-mixed compartment.

**Model Kind**: `OneCompIVBolus`

**Parameters**: `OneCompIVBolusParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `CL` | CL | volume/time | Clearance | CL > 0 |
| `V` | V | volume | Volume of distribution | V > 0 |

**State Variables**:

| State | Symbol | Description |
|-------|--------|-------------|
| `A_central` | A | Amount in central compartment |

**Differential Equation**:

$$\frac{dA}{dt} = -\frac{CL}{V} \cdot A$$

**Observation**:

$$C_{conc} = \frac{A}{V}$$

**Dose Target**: Central compartment (index 1)

**Example**:

```julia
using OpenPKPDCore

# Parameters: CL = 5 L/h, V = 50 L
params = OneCompIVBolusParams(5.0, 50.0)

# Single 100 mg dose at t=0
doses = [DoseEvent(0.0, 100.0)]

# Create model specification
spec = ModelSpec(OneCompIVBolus(), "iv_bolus_example", params, doses)

# Simulation grid: 0 to 24 hours
grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))

# High-precision solver
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run simulation
result = simulate(spec, grid, solver)

# Initial concentration: dose/V = 100/50 = 2 mg/L
println("C(0) = ", result.observations[:conc][1])  # 2.0

# Half-life: t½ = 0.693 * V / CL = 0.693 * 50 / 5 = 6.93 h
```

**Multiple Doses**:

```julia
# Multiple doses: 100 mg at t=0, 50 mg at t=12
doses = [
    DoseEvent(0.0, 100.0),
    DoseEvent(12.0, 50.0)
]

spec = ModelSpec(OneCompIVBolus(), "multi_dose", params, doses)
result = simulate(spec, grid, solver)
```

---

### One-Compartment Oral First-Order

One-compartment model with first-order absorption from a depot (gut) compartment.

**Model Kind**: `OneCompOralFirstOrder`

**Parameters**: `OneCompOralFirstOrderParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `Ka` | Ka | 1/time | Absorption rate constant | Ka > 0 |
| `CL` | CL | volume/time | Clearance | CL > 0 |
| `V` | V | volume | Volume of distribution | V > 0 |

**State Variables**:

| State | Symbol | Description |
|-------|--------|-------------|
| `A_gut` | A_gut | Amount in absorption depot |
| `A_central` | A_central | Amount in central compartment |

**Differential Equations**:

$$\frac{dA_{gut}}{dt} = -K_a \cdot A_{gut}$$

$$\frac{dA_{central}}{dt} = K_a \cdot A_{gut} - \frac{CL}{V} \cdot A_{central}$$

**Observation**:

$$C_{conc} = \frac{A_{central}}{V}$$

**Dose Target**: Gut compartment (index 1)

**Example**:

```julia
using OpenPKPDCore

# Parameters: Ka = 1.5 /h, CL = 5 L/h, V = 50 L
params = OneCompOralFirstOrderParams(1.5, 5.0, 50.0)

# 200 mg oral dose at t=0
doses = [DoseEvent(0.0, 200.0)]

spec = ModelSpec(OneCompOralFirstOrder(), "oral_example", params, doses)
grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

result = simulate(spec, grid, solver)

# Concentration starts at 0 and rises to Cmax before declining
println("C(0) = ", result.observations[:conc][1])   # ~0
println("Times: ", result.t)
println("Concentrations: ", result.observations[:conc])
```

**Pharmacokinetic Metrics**:

```julia
# Calculate Cmax and AUC
t = result.t
c = result.observations[:conc]

cmax_val = cmax(t, c)
auc_val = auc_trapezoid(t, c)

println("Cmax = ", cmax_val)
println("AUC(0-24) = ", auc_val)
```

---

## Pharmacodynamic Models

### Direct Emax Model

A direct effect model where drug effect is an instantaneous function of concentration.

**Model Kind**: `DirectEmax`

**Parameters**: `DirectEmaxParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `E0` | E0 | effect units | Baseline effect | - |
| `Emax` | Emax | effect units | Maximum effect | Emax > 0 |
| `EC50` | EC50 | concentration | Concentration at 50% Emax | EC50 > 0 |

**Effect Equation**:

$$E(C) = E_0 + \frac{E_{max} \cdot C}{EC_{50} + C}$$

**Example - Sequential PKPD**:

```julia
using OpenPKPDCore

# PK model
pk_params = OneCompIVBolusParams(5.0, 50.0)
pk_spec = ModelSpec(OneCompIVBolus(), "pk", pk_params, [DoseEvent(0.0, 100.0)])

# PD model: E0=10, Emax=90, EC50=1.0 mg/L
pd_params = DirectEmaxParams(10.0, 90.0, 1.0)
pd_spec = PDSpec(DirectEmax(), "pd", pd_params, :conc, :effect)

grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run sequential PKPD (PK solved first, then PD evaluated)
result = simulate_pkpd(pk_spec, pd_spec, grid, solver)

println("Concentrations: ", result.observations[:conc])
println("Effects: ", result.observations[:effect])
```

---

### Indirect Response Turnover Model

An indirect response model where drug inhibits the production or elimination of a response variable.

**Model Kind**: `IndirectResponseTurnover`

**Parameters**: `IndirectResponseTurnoverParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `Kin` | Kin | response/time | Zero-order input rate | Kin > 0 |
| `Kout` | Kout | 1/time | First-order output rate | Kout > 0 |
| `R0` | R0 | response units | Baseline response | R0 > 0 |
| `Imax` | Imax | - | Maximum inhibition | 0 ≤ Imax ≤ 1 |
| `IC50` | IC50 | concentration | Concentration at 50% Imax | IC50 > 0 |

**State Variable**:

| State | Symbol | Description |
|-------|--------|-------------|
| `R` | R | Response |

**Differential Equation**:

$$\frac{dR}{dt} = K_{in} - K_{out} \cdot (1 - I(C)) \cdot R$$

**Inhibition Function**:

$$I(C) = \frac{I_{max} \cdot C}{IC_{50} + C}$$

**Note**: At steady state without drug: $R_0 = K_{in} / K_{out}$

**Example - Coupled PKPD**:

```julia
using OpenPKPDCore

# PK model
pk_params = OneCompIVBolusParams(5.0, 50.0)
pk_spec = ModelSpec(OneCompIVBolus(), "pk", pk_params, [DoseEvent(0.0, 100.0)])

# PD model: Kin=10, Kout=0.5, R0=20, Imax=0.8, IC50=1.0
pd_params = IndirectResponseTurnoverParams(10.0, 0.5, 20.0, 0.8, 1.0)
pd_spec = PDSpec(IndirectResponseTurnover(), "indirect_pd", pd_params, :conc, :response)

grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run coupled PKPD (simultaneous PK-PD ODE system)
result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)

println("Concentrations: ", result.observations[:conc])
println("Response: ", result.observations[:response])
```

---

## Dose Events

Doses are specified using the `DoseEvent` struct.

**Structure**:

```julia
DoseEvent(time::Float64, amount::Float64)
```

| Field | Description |
|-------|-------------|
| `time` | Time of dose administration |
| `amount` | Dose amount (mass units) |

**Dose Handling**:

- Doses at `t0` are added to the initial condition
- Doses in `(t0, t1]` are handled via callbacks
- Multiple doses at the same time are summed
- Doses outside `[t0, t1]` are ignored

**Example**:

```julia
# QD dosing for 7 days
doses = [DoseEvent(Float64(d * 24), 100.0) for d in 0:6]

# BID dosing
doses = [
    DoseEvent(0.0, 50.0),
    DoseEvent(12.0, 50.0),
    DoseEvent(24.0, 50.0),
    DoseEvent(36.0, 50.0)
]
```

---

## Simulation Grid

The `SimGrid` struct defines the time domain for simulation.

**Structure**:

```julia
SimGrid(t0::Float64, t1::Float64, saveat::Vector{Float64})
```

| Field | Description | Constraints |
|-------|-------------|-------------|
| `t0` | Start time | - |
| `t1` | End time | t1 > t0 |
| `saveat` | Output time points | Sorted, within [t0, t1] |

**Example**:

```julia
# Hourly output from 0 to 24 hours
grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))

# Fine resolution around Cmax
grid = SimGrid(0.0, 24.0, collect(0.0:0.1:24.0))

# Sparse sampling
grid = SimGrid(0.0, 168.0, [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0, 72.0, 168.0])
```

---

## Solver Specification

The `SolverSpec` struct configures the ODE solver.

**Structure**:

```julia
SolverSpec(alg::Symbol, reltol::Float64, abstol::Float64, maxiters::Int)
```

| Field | Description | Recommended Values |
|-------|-------------|-------------------|
| `alg` | Solver algorithm | `:Tsit5` (non-stiff), `:Rosenbrock23` (stiff) |
| `reltol` | Relative tolerance | 1e-6 to 1e-10 |
| `abstol` | Absolute tolerance | 1e-9 to 1e-12 |
| `maxiters` | Maximum iterations | 10000 to 10000000 |

**Golden Standard Settings**:

```julia
# High precision for reproducibility
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)
```

---

## Simulation Results

The `SimResult` struct contains simulation output.

**Structure**:

```julia
struct SimResult
    t::Vector{Float64}
    states::Dict{Symbol, Vector{Float64}}
    observations::Dict{Symbol, Vector{Float64}}
    metadata::Dict{String, Any}
end
```

| Field | Description |
|-------|-------------|
| `t` | Time points (matches `saveat`) |
| `states` | State variable trajectories |
| `observations` | Observable outputs (`:conc`, `:effect`, etc.) |
| `metadata` | Simulation metadata |

**Metadata Contents**:

- `model`: Model name
- `solver_alg`: Algorithm used
- `reltol`, `abstol`: Tolerances
- `dose_schedule`: Applied doses
- `engine_version`: OpenPKPD version
- `event_semantics_version`: Event handling version
- `solver_semantics_version`: Solver semantics version

---

## Exposure Metrics

OpenPKPD provides functions for common pharmacokinetic metrics.

### Maximum Concentration (Cmax)

```julia
cmax_value = cmax(t, c)
```

Returns the maximum value in the concentration vector.

### Area Under the Curve (AUC)

```julia
auc_value = auc_trapezoid(t, c)
```

Computes AUC using the trapezoidal rule. Requires sorted time points.

**Example**:

```julia
result = simulate(spec, grid, solver)

t = result.t
c = result.observations[:conc]

println("Cmax = ", cmax(t, c))
println("AUC = ", auc_trapezoid(t, c))
```
