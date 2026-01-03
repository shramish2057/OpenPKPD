# PK/PD Models Reference

OpenPKPD provides a comprehensive library of validated pharmacokinetic (PK) and pharmacodynamic (PD) models for drug concentration and effect simulation.

## Model Summary

### Pharmacokinetic Models

| Model | Compartments | Parameters | Route | Use Case |
|-------|--------------|------------|-------|----------|
| `OneCompIVBolus` | 1 | CL, V | IV bolus | Simple IV drugs |
| `OneCompOralFirstOrder` | 2 | Ka, CL, V | Oral | Simple oral drugs |
| `TwoCompIVBolus` | 2 | CL, V1, Q, V2 | IV bolus | Distribution kinetics |
| `TwoCompOral` | 3 | Ka, CL, V1, Q, V2 | Oral | Oral with distribution |
| `ThreeCompIVBolus` | 3 | CL, V1, Q2, V2, Q3, V3 | IV bolus | Deep tissue distribution |
| `TransitAbsorption` | N+1 | N, Ktr, Ka, CL, V | Oral | Delayed absorption |
| `MichaelisMentenElimination` | 1 | Vmax, Km, V | IV bolus | Saturable elimination |

### Pharmacodynamic Models

| Model | Type | Parameters | Mechanism |
|-------|------|------------|-----------|
| `DirectEmax` | Direct | E0, Emax, EC50 | Hyperbolic response |
| `SigmoidEmax` | Direct | E0, Emax, EC50, gamma | Hill equation |
| `IndirectResponseTurnover` | Indirect | Kin, Kout, R0, Imax, IC50 | Inhibition of output |
| `BiophaseEquilibration` | Effect compartment | ke0, E0, Emax, EC50 | PK-PD hysteresis |

---

## Pharmacokinetic Models

### One-Compartment IV Bolus

The simplest PK model representing instantaneous drug administration into a single well-mixed compartment.

**Model Kind**: `OneCompIVBolus`

**Parameters**: `OneCompIVBolusParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `CL` | CL | L/h | Clearance | CL > 0 |
| `V` | V | L | Volume of distribution | V > 0 |

**State Variables**:

| State | Symbol | Description |
|-------|--------|-------------|
| `A_central` | A | Amount in central compartment |

**Differential Equation**:

$$\frac{dA}{dt} = -\frac{CL}{V} \cdot A$$

**Observation**: $C = \frac{A}{V}$

**Half-life**: $t_{1/2} = \frac{0.693 \cdot V}{CL}$

**Example (Julia)**:

```julia
using OpenPKPDCore

params = OneCompIVBolusParams(5.0, 50.0)  # CL=5 L/h, V=50 L
doses = [DoseEvent(0.0, 100.0)]
spec = ModelSpec(OneCompIVBolus(), "iv_bolus", params, doses)

grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

result = simulate(spec, grid, solver)
```

**Example (Python)**:

```python
import openpkpd

openpkpd.init_julia()

result = openpkpd.simulate_pk_iv_bolus(
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=24.0,
    saveat=[float(t) for t in range(25)]
)
```

---

### One-Compartment Oral First-Order

One-compartment model with first-order absorption from a depot (gut) compartment.

**Model Kind**: `OneCompOralFirstOrder`

**Parameters**: `OneCompOralFirstOrderParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `Ka` | Ka | 1/h | Absorption rate constant | Ka > 0 |
| `CL` | CL | L/h | Clearance | CL > 0 |
| `V` | V | L | Volume of distribution | V > 0 |

**State Variables**:

| State | Symbol | Description |
|-------|--------|-------------|
| `A_gut` | A_gut | Amount in absorption depot |
| `A_central` | A_central | Amount in central compartment |

**Differential Equations**:

$$\frac{dA_{gut}}{dt} = -K_a \cdot A_{gut}$$

$$\frac{dA_{central}}{dt} = K_a \cdot A_{gut} - \frac{CL}{V} \cdot A_{central}$$

**Observation**: $C = \frac{A_{central}}{V}$

**Dose Target**: Gut compartment

**Tmax (approximate)**: $t_{max} \approx \frac{\ln(K_a) - \ln(k_{el})}{K_a - k_{el}}$ where $k_{el} = CL/V$

**Example (Julia)**:

```julia
params = OneCompOralFirstOrderParams(1.5, 5.0, 50.0)  # Ka=1.5/h, CL=5, V=50
doses = [DoseEvent(0.0, 200.0)]
spec = ModelSpec(OneCompOralFirstOrder(), "oral", params, doses)
result = simulate(spec, grid, solver)
```

**Example (Python)**:

```python
result = openpkpd.simulate_pk_oral_first_order(
    ka=1.5,
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 200.0}],
    t0=0.0,
    t1=24.0
)
```

---

### Two-Compartment IV Bolus

Two-compartment model with central and peripheral compartments, exhibiting bi-exponential decline.

**Model Kind**: `TwoCompIVBolus`

**Parameters**: `TwoCompIVBolusParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `CL` | CL | L/h | Clearance from central | CL > 0 |
| `V1` | V1 | L | Central volume | V1 > 0 |
| `Q` | Q | L/h | Inter-compartmental clearance | Q > 0 |
| `V2` | V2 | L | Peripheral volume | V2 > 0 |

**Micro-rate Constants**:

- $k_{10} = CL/V_1$ (elimination)
- $k_{12} = Q/V_1$ (central to peripheral)
- $k_{21} = Q/V_2$ (peripheral to central)

**State Variables**:

| State | Symbol | Description |
|-------|--------|-------------|
| `A_central` | A1 | Amount in central compartment |
| `A_peripheral` | A2 | Amount in peripheral compartment |

**Differential Equations**:

$$\frac{dA_1}{dt} = -k_{10} \cdot A_1 - k_{12} \cdot A_1 + k_{21} \cdot A_2$$

$$\frac{dA_2}{dt} = k_{12} \cdot A_1 - k_{21} \cdot A_2$$

**Observation**: $C = \frac{A_1}{V_1}$

**Bi-exponential Profile**:

$$C(t) = A \cdot e^{-\alpha t} + B \cdot e^{-\beta t}$$

Where $\alpha$ (distribution) > $\beta$ (elimination)

**Example (Julia)**:

```julia
params = TwoCompIVBolusParams(10.0, 20.0, 15.0, 50.0)
# CL=10 L/h, V1=20 L, Q=15 L/h, V2=50 L
doses = [DoseEvent(0.0, 500.0)]
spec = ModelSpec(TwoCompIVBolus(), "twocomp_iv", params, doses)
result = simulate(spec, grid, solver)
```

**Example (Python)**:

```python
result = openpkpd.simulate_pk_twocomp_iv_bolus(
    cl=10.0,
    v1=20.0,
    q=15.0,
    v2=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0,
    t1=48.0
)
```

---

### Two-Compartment Oral

Two-compartment model with first-order oral absorption.

**Model Kind**: `TwoCompOral`

**Parameters**: `TwoCompOralParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `Ka` | Ka | 1/h | Absorption rate constant | Ka > 0 |
| `CL` | CL | L/h | Clearance from central | CL > 0 |
| `V1` | V1 | L | Central volume | V1 > 0 |
| `Q` | Q | L/h | Inter-compartmental clearance | Q > 0 |
| `V2` | V2 | L | Peripheral volume | V2 > 0 |

**State Variables**:

| State | Symbol | Description |
|-------|--------|-------------|
| `A_gut` | A_gut | Amount in gut compartment |
| `A_central` | A1 | Amount in central compartment |
| `A_peripheral` | A2 | Amount in peripheral compartment |

**Differential Equations**:

$$\frac{dA_{gut}}{dt} = -K_a \cdot A_{gut}$$

$$\frac{dA_1}{dt} = K_a \cdot A_{gut} - \frac{CL}{V_1} \cdot A_1 - \frac{Q}{V_1} \cdot A_1 + \frac{Q}{V_2} \cdot A_2$$

$$\frac{dA_2}{dt} = \frac{Q}{V_1} \cdot A_1 - \frac{Q}{V_2} \cdot A_2$$

**Example (Python)**:

```python
result = openpkpd.simulate_pk_twocomp_oral(
    ka=1.2,
    cl=8.0,
    v1=25.0,
    q=12.0,
    v2=60.0,
    doses=[{"time": 0.0, "amount": 400.0}],
    t0=0.0,
    t1=72.0
)
```

---

### Three-Compartment IV Bolus

Three-compartment mammillary model with central, shallow peripheral, and deep peripheral compartments.

**Model Kind**: `ThreeCompIVBolus`

**Parameters**: `ThreeCompIVBolusParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `CL` | CL | L/h | Clearance from central | CL > 0 |
| `V1` | V1 | L | Central volume | V1 > 0 |
| `Q2` | Q2 | L/h | Clearance to shallow peripheral | Q2 > 0 |
| `V2` | V2 | L | Shallow peripheral volume | V2 > 0 |
| `Q3` | Q3 | L/h | Clearance to deep peripheral | Q3 > 0 |
| `V3` | V3 | L | Deep peripheral volume | V3 > 0 |

**State Variables**:

| State | Symbol | Description |
|-------|--------|-------------|
| `A_central` | A1 | Central compartment |
| `A_periph1` | A2 | Shallow peripheral (rapid equilibration) |
| `A_periph2` | A3 | Deep peripheral (slow equilibration) |

**Micro-rate Constants**:

- $k_{10} = CL/V_1$
- $k_{12} = Q_2/V_1$, $k_{21} = Q_2/V_2$
- $k_{13} = Q_3/V_1$, $k_{31} = Q_3/V_3$

**Tri-exponential Profile**:

$$C(t) = A \cdot e^{-\alpha t} + B \cdot e^{-\beta t} + C \cdot e^{-\gamma t}$$

Where $\alpha$ > $\beta$ > $\gamma$

**Use Cases**: Lipophilic drugs, anesthetics, drugs with deep tissue binding

**Example (Python)**:

```python
result = openpkpd.simulate_pk_threecomp_iv_bolus(
    cl=5.0,
    v1=10.0,
    q2=20.0,
    v2=30.0,
    q3=5.0,
    v3=100.0,
    doses=[{"time": 0.0, "amount": 1000.0}],
    t0=0.0,
    t1=168.0
)
```

---

### Transit Absorption

Models delayed and complex oral absorption using a chain of transit compartments.

**Model Kind**: `TransitAbsorption`

**Parameters**: `TransitAbsorptionParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `N` | N | - | Number of transit compartments | 1 ≤ N ≤ 20 |
| `Ktr` | Ktr | 1/h | Transit rate constant | Ktr > 0 |
| `Ka` | Ka | 1/h | Absorption rate from last transit | Ka > 0 |
| `CL` | CL | L/h | Clearance | CL > 0 |
| `V` | V | L | Volume of distribution | V > 0 |

**State Variables**: N+1 total

- Transit_1 through Transit_N
- A_central

**Differential Equations**:

$$\frac{dTransit_1}{dt} = -K_{tr} \cdot Transit_1$$ (receives dose)

$$\frac{dTransit_i}{dt} = K_{tr} \cdot Transit_{i-1} - K_{tr} \cdot Transit_i$$ (i > 1)

$$\frac{dA_{central}}{dt} = K_a \cdot Transit_N - \frac{CL}{V} \cdot A_{central}$$

**Mean Transit Time (MTT)**: $(N + 1) / K_{tr}$

**Absorption Profile**: Gamma-like distribution (delayed peak, allows for lag time)

**Use Cases**:
- Controlled-release formulations
- Enteric-coated tablets
- Drugs with complex GI transit
- When simple lag time is insufficient

**Example (Python)**:

```python
result = openpkpd.simulate_pk_transit_absorption(
    n_transit=5,
    ktr=0.5,
    ka=2.0,
    cl=10.0,
    v=70.0,
    doses=[{"time": 0.0, "amount": 300.0}],
    t0=0.0,
    t1=24.0
)
```

---

### Michaelis-Menten Elimination

Models saturable (capacity-limited) elimination kinetics.

**Model Kind**: `MichaelisMentenElimination`

**Parameters**: `MichaelisMentenEliminationParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `Vmax` | Vmax | mg/h | Maximum elimination rate | Vmax > 0 |
| `Km` | Km | mg/L | Concentration at half Vmax | Km > 0 |
| `V` | V | L | Volume of distribution | V > 0 |

**Differential Equation**:

$$\frac{dA}{dt} = -\frac{V_{max} \cdot A}{K_m \cdot V + A}$$

Or in terms of concentration:

$$\frac{dC}{dt} = -\frac{V_{max} \cdot C}{K_m + C}$$

**Kinetic Behavior**:

| Condition | Kinetics | Effective CL |
|-----------|----------|--------------|
| C << Km | First-order (linear) | CL ≈ Vmax/Km |
| C >> Km | Zero-order (saturable) | Rate ≈ Vmax |
| C ≈ Km | Mixed-order | Nonlinear |

**Characteristics**:
- Dose-dependent half-life (increases with dose)
- Disproportionate increase in AUC with dose
- Time to steady state depends on dose

**Clinical Examples**: Phenytoin, ethanol, high-dose aspirin, some biologics

**Example (Python)**:

```python
result = openpkpd.simulate_pk_michaelis_menten(
    vmax=500.0,    # mg/h
    km=10.0,       # mg/L
    v=50.0,        # L
    doses=[{"time": 0.0, "amount": 1000.0}],
    t0=0.0,
    t1=48.0
)
```

---

## Pharmacodynamic Models

### Direct Emax Model

A direct effect model where drug effect is an instantaneous function of concentration.

**Model Kind**: `DirectEmax`

**Parameters**: `DirectEmaxParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `E0` | E0 | effect units | Baseline effect | Any real |
| `Emax` | Emax | effect units | Maximum effect | Emax > 0 |
| `EC50` | EC50 | mg/L | Concentration at 50% Emax | EC50 > 0 |

**Effect Equation**:

$$E(C) = E_0 + \frac{E_{max} \cdot C}{EC_{50} + C}$$

**Properties**:
- At C = 0: E = E0
- At C = EC50: E = E0 + Emax/2
- As C → ∞: E → E0 + Emax

**Use Cases**: Rapid equilibrium effects, receptor binding without delays

**Example (Python)**:

```python
result = openpkpd.simulate_pkpd_direct_emax(
    # PK parameters
    cl=5.0, v=50.0,
    # PD parameters
    e0=0.0, emax=100.0, ec50=2.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0
)
```

---

### Sigmoid Emax (Hill Equation)

Extension of Emax model with Hill coefficient for steepness control.

**Model Kind**: `SigmoidEmax`

**Parameters**: `SigmoidEmaxParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `E0` | E0 | effect units | Baseline effect | Any real |
| `Emax` | Emax | effect units | Maximum effect | Any (negative for inhibition) |
| `EC50` | EC50 | mg/L | Concentration at 50% Emax | EC50 > 0 |
| `gamma` | γ | - | Hill coefficient (steepness) | γ > 0 |

**Effect Equation**:

$$E(C) = E_0 + \frac{E_{max} \cdot C^{\gamma}}{EC_{50}^{\gamma} + C^{\gamma}}$$

**Hill Coefficient Interpretation**:

| γ Value | Response Shape | Example |
|---------|----------------|---------|
| γ = 1 | Hyperbolic (standard Emax) | Many drugs |
| γ < 1 | Gradual, shallow curve | Some enzymes |
| γ > 1 | Steep, switch-like | Neuromuscular blockers |
| γ = 2-4 | Moderate steepness | Ion channels |
| γ > 4 | Near all-or-none | Cooperative binding |

**Example (Python)**:

```python
result = openpkpd.simulate_pkpd_sigmoid_emax(
    cl=5.0, v=50.0,
    e0=10.0, emax=90.0, ec50=1.5, gamma=2.5,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0
)
```

---

### Indirect Response Turnover Model

Models drug effects on turnover of endogenous substances.

**Model Kind**: `IndirectResponseTurnover`

**Parameters**: `IndirectResponseTurnoverParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `Kin` | Kin | units/h | Zero-order input rate | Kin > 0 |
| `Kout` | Kout | 1/h | First-order output rate | Kout > 0 |
| `R0` | R0 | units | Baseline response | R0 > 0 |
| `Imax` | Imax | - | Maximum inhibition | 0 ≤ Imax ≤ 1 |
| `IC50` | IC50 | mg/L | Concentration at 50% Imax | IC50 > 0 |

**Note**: At baseline, R0 = Kin/Kout

**Inhibition Function**:

$$I(C) = \frac{I_{max} \cdot C}{IC_{50} + C}$$

**Differential Equation** (Type IV - inhibition of Kout):

$$\frac{dR}{dt} = K_{in} - K_{out} \cdot (1 - I(C)) \cdot R$$

**Indirect Response Types**:

| Type | Mechanism | Drug Effect |
|------|-----------|-------------|
| I | Inhibit Kin | Decreases R (e.g., warfarin on clotting factors) |
| II | Inhibit Kout | Increases R (e.g., corticosteroids on cortisol) |
| III | Stimulate Kin | Increases R (e.g., EPO on RBCs) |
| IV | Stimulate Kout | Decreases R |

**Characteristics**:
- Delayed response relative to concentration
- Counter-clockwise hysteresis in effect-concentration plots
- Time to maximum effect often occurs after Cmax

**Example (Python)**:

```python
result = openpkpd.simulate_pkpd_indirect_response(
    cl=5.0, v=50.0,
    kin=10.0, kout=0.5, r0=20.0, imax=0.8, ic50=1.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=72.0
)
```

---

### Biophase Equilibration (Effect Compartment)

Models temporal delay between plasma and effect site concentrations.

**Model Kind**: `BiophaseEquilibration`

**Parameters**: `BiophaseEquilibrationParams`

| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| `ke0` | ke0 | 1/h | Effect site equilibration rate | ke0 > 0 |
| `E0` | E0 | effect units | Baseline effect | Any real |
| `Emax` | Emax | effect units | Maximum effect | Emax > 0 |
| `EC50` | EC50 | mg/L | Effect site conc at 50% Emax | EC50 > 0 |

**Effect Compartment ODE**:

$$\frac{dC_e}{dt} = k_{e0} \cdot (C_p - C_e)$$

Where $C_p$ = plasma concentration, $C_e$ = effect site concentration

**Effect Equation** (based on Ce):

$$E(C_e) = E_0 + \frac{E_{max} \cdot C_e}{EC_{50} + C_e}$$

**Equilibration Kinetics**:
- Half-life of equilibration: $t_{1/2,ke0} = \ln(2) / k_{e0}$
- Time to 90% equilibration: $t_{90} \approx 2.3 / k_{e0}$

**Use Cases**:
- CNS-active drugs (anesthetics, sedatives)
- Drugs with hysteresis in effect-concentration plots
- Neuromuscular blocking agents

**Example (Python)**:

```python
result = openpkpd.simulate_pkpd_biophase_equilibration(
    cl=5.0, v=50.0,
    ke0=0.5,  # Equilibration half-life ~1.4 hours
    e0=0.0, emax=100.0, ec50=2.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0
)
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

**Python format**:

```python
doses = [
    {"time": 0.0, "amount": 100.0},
    {"time": 24.0, "amount": 100.0},
]
```

**Dose Handling Rules**:

- Doses at `t0` are added to the initial condition
- Doses in `(t0, t1]` are handled via callbacks
- Multiple doses at the same time are summed
- Doses outside `[t0, t1]` are ignored

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

**Examples**:

```julia
# Hourly output
grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))

# Fine resolution
grid = SimGrid(0.0, 24.0, collect(0.0:0.1:24.0))

# Sparse PK sampling
grid = SimGrid(0.0, 168.0, [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0, 72.0, 168.0])
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

**Golden Standard Settings** (for reproducibility):

```julia
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)
```

---

## See Also

- [Population Simulation](population.md) - IIV, IOV, and covariate modeling
- [NCA Reference](nca.md) - Non-compartmental analysis
- [Trial Simulation](trial.md) - Clinical trial simulation
- [Python Bindings](python.md) - Complete Python API reference
