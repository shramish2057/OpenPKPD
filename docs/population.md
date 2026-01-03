# Population Simulation

OpenPKPD supports population pharmacokinetic/pharmacodynamic simulations with multiple sources of variability.

## Overview

Population simulations model drug behavior across multiple individuals, accounting for:

- **Inter-Individual Variability (IIV)**: Differences between individuals
- **Inter-Occasion Variability (IOV)**: Differences within an individual across occasions
- **Covariate Effects**: Patient characteristics affecting parameters
- **Time-Varying Covariates**: Dynamic patient characteristics

## Population Specification

The `PopulationSpec` struct defines a population simulation.

```julia
struct PopulationSpec{MS}
    base_model_spec::MS         # Typical value model
    iiv::Union{Nothing, IIVSpec}           # Inter-individual variability
    iov::Union{Nothing, IOVSpec}           # Inter-occasion variability
    covariate_model::Union{Nothing, CovariateModel}  # Covariate model
    covariates::Vector{IndividualCovariates}         # Per-individual covariates
end
```

---

## Inter-Individual Variability (IIV)

IIV represents random differences in PK/PD parameters between individuals.

### Log-Normal Distribution

OpenPKPD implements log-normal IIV, the standard in pharmacometrics:

$$\theta_i = \theta_{pop} \cdot e^{\eta_i}$$

where $\eta_i \sim N(0, \omega^2)$

### IIV Specification

```julia
struct IIVSpec{K<:RandomEffectKind}
    kind::K                      # LogNormalIIV()
    omegas::Dict{Symbol, Float64}  # Parameter-specific omegas
    seed::UInt64                 # RNG seed for reproducibility
    n::Int                       # Number of individuals
end
```

| Field | Description |
|-------|-------------|
| `kind` | Distribution type (`LogNormalIIV()`) |
| `omegas` | Standard deviations for each parameter |
| `seed` | Deterministic random seed |
| `n` | Population size |

### Example: Basic IIV

=== "Julia"

    ```julia
    using OpenPKPDCore

    # Base model (typical values)
    base_params = OneCompIVBolusParams(5.0, 50.0)  # CL=5, V=50
    base_spec = ModelSpec(
        OneCompIVBolus(),
        "pop_iv",
        base_params,
        [DoseEvent(0.0, 100.0)]
    )

    # IIV specification: 30% CV on CL, 20% CV on V
    # omega = sqrt(ln(1 + CV^2)) ≈ CV for small CV
    iiv = IIVSpec(
        LogNormalIIV(),
        Dict(:CL => 0.3, :V => 0.2),  # omegas
        UInt64(12345),                 # seed
        100                            # n individuals
    )

    # Population specification (no IOV, no covariates)
    pop = PopulationSpec(base_spec, iiv, nothing, nothing, [])

    # Run simulation
    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

    result = simulate_population(pop, grid, solver)

    # Access results
    println("Number of individuals: ", length(result.individuals))
    println("First individual CL: ", result.params[1][:CL])
    println("Mean concentration at t=0: ", result.summaries[:conc].mean[1])
    ```

=== "Python"

    ```python
    import openpkpd

    openpkpd.init_julia()

    # Run population simulation with IIV
    result = openpkpd.simulate_population_iv_bolus(
        cl=5.0,              # Typical CL
        v=50.0,              # Typical V
        doses=[{"time": 0.0, "amount": 100.0}],
        t0=0.0,
        t1=24.0,
        saveat=[t * 0.5 for t in range(49)],
        n=100,               # Number of individuals
        seed=12345,          # RNG seed for reproducibility
        omegas={             # IIV standard deviations (30% on CL, 20% on V)
            "CL": 0.3,
            "V": 0.2
        }
    )

    # Access results
    print("Number of individuals:", len(result["individuals"]))
    print("First individual CL:", result["params"][0]["CL"])
    print("Mean concentration at t=0:", result["summaries"]["conc"]["mean"][0])
    ```

### Population Results

The `PopulationResult` struct contains:

```julia
struct PopulationResult
    individuals::Vector{SimResult}              # Per-individual results
    params::Vector{Dict{Symbol, Float64}}       # Realized parameters
    summaries::Dict{Symbol, PopulationSummary}  # Summary statistics
    metadata::Dict{String, Any}                 # Metadata
end
```

### Population Summary

```julia
struct PopulationSummary
    observation::Symbol                        # Which observation
    probs::Vector{Float64}                     # Quantile probabilities
    mean::Vector{Float64}                      # Mean at each time
    median::Vector{Float64}                    # Median at each time
    quantiles::Dict{Float64, Vector{Float64}}  # Quantiles at each time
end
```

**Example - Accessing Summaries**:

```julia
summary = result.summaries[:conc]

println("Mean concentrations: ", summary.mean)
println("Median concentrations: ", summary.median)
println("5th percentile: ", summary.quantiles[0.05])
println("95th percentile: ", summary.quantiles[0.95])
```

---

## Inter-Occasion Variability (IOV)

IOV represents random differences in parameters between dosing occasions within an individual.

### IOV Transform

$$\theta_{occasion} = \theta_{base} \cdot e^{\kappa}$$

where $\kappa \sim N(0, \pi^2)$

### Occasion Definition

Occasions are automatically derived from dose times:

- `t0` starts Occasion 1
- Each unique dose time after `t0` starts a new occasion

```julia
struct OccasionDefinition
    mode::Symbol  # Currently only :dose_times
end
```

### IOV Specification

```julia
struct IOVSpec{K<:RandomEffectKind}
    kind::K                       # LogNormalIIV()
    pis::Dict{Symbol, Float64}    # Parameter-specific standard deviations
    seed::UInt64                  # IOV RNG seed (separate from IIV)
    occasion_def::OccasionDefinition
end
```

### Example: IIV + IOV

```julia
using OpenPKPDCore

# Base model with multiple doses (defining occasions)
base_params = OneCompIVBolusParams(5.0, 50.0)
doses = [
    DoseEvent(0.0, 100.0),   # Occasion 1
    DoseEvent(24.0, 100.0),  # Occasion 2
    DoseEvent(48.0, 100.0)   # Occasion 3
]
base_spec = ModelSpec(OneCompIVBolus(), "iov_example", base_params, doses)

# IIV: 30% on CL
iiv = IIVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.3),
    UInt64(12345),
    50
)

# IOV: 15% on CL between occasions
iov = IOVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.15),
    UInt64(67890),  # Different seed than IIV
    OccasionDefinition(:dose_times)
)

# Population with both IIV and IOV
pop = PopulationSpec(base_spec, iiv, iov, nothing, [])

grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

result = simulate_population(pop, grid, solver)
```

### IOV Implementation Details

- IOV uses **segmented simulation**: separate ODE solves per occasion
- **State continuity** is maintained at occasion boundaries
- Parameters are re-sampled at each occasion start
- IOV seed is separate from IIV seed for independent sampling

---

## Covariate Effects

Covariates model the influence of patient characteristics on PK/PD parameters.

### Covariate Effect Types

#### Linear Covariate

$$\theta_i = \theta_{pop} \cdot (1 + \beta \cdot (cov - ref))$$

```julia
LinearCovariate()
```

#### Power Covariate

$$\theta_i = \theta_{pop} \cdot \left(\frac{cov}{ref}\right)^\beta$$

```julia
PowerCovariate()
```

#### Exponential Covariate

$$\theta_i = \theta_{pop} \cdot e^{\beta \cdot (cov - ref)}$$

```julia
ExpCovariate()
```

### Covariate Effect Structure

```julia
struct CovariateEffect{K<:CovariateEffectKind}
    kind::K           # LinearCovariate, PowerCovariate, or ExpCovariate
    param::Symbol     # Parameter being modified
    covariate::Symbol # Covariate variable
    beta::Float64     # Effect coefficient
    ref::Float64      # Reference covariate value
end
```

### Covariate Model

```julia
struct CovariateModel
    name::String
    effects::Vector{CovariateEffect}  # Applied in order
end
```

### Individual Covariates

```julia
struct IndividualCovariates
    values::Dict{Symbol, Float64}                    # Static covariates
    time_varying::Union{Nothing, TimeVaryingCovariates}  # Dynamic covariates
end
```

### Example: Weight-Based Covariate

```julia
using OpenPKPDCore

# Base model
base_params = OneCompIVBolusParams(5.0, 50.0)
base_spec = ModelSpec(
    OneCompIVBolus(),
    "cov_example",
    base_params,
    [DoseEvent(0.0, 100.0)]
)

# Covariate model: CL scales with weight^0.75 (allometric)
cov_model = CovariateModel(
    "weight_model",
    [
        CovariateEffect(PowerCovariate(), :CL, :WT, 0.75, 70.0),  # CL ~ (WT/70)^0.75
        CovariateEffect(LinearCovariate(), :V, :WT, 0.01, 70.0)   # V ~ 1 + 0.01*(WT-70)
    ]
)

# Individual covariates (3 patients with different weights)
covariates = [
    IndividualCovariates(Dict(:WT => 50.0), nothing),
    IndividualCovariates(Dict(:WT => 70.0), nothing),
    IndividualCovariates(Dict(:WT => 100.0), nothing)
]

# Population without IIV (covariate effects only)
pop = PopulationSpec(base_spec, nothing, nothing, cov_model, covariates)

grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

result = simulate_population(pop, grid, solver)

# Each individual has different parameters based on weight
for (i, p) in enumerate(result.params)
    println("Individual $i: CL=$(p[:CL]), V=$(p[:V])")
end
```

### Example: Combined IIV + Covariates

```julia
# IIV specification
iiv = IIVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.3, :V => 0.2),
    UInt64(12345),
    100
)

# Covariate model
cov_model = CovariateModel(
    "full_model",
    [
        CovariateEffect(PowerCovariate(), :CL, :WT, 0.75, 70.0),
        CovariateEffect(LinearCovariate(), :CL, :AGE, -0.01, 40.0)
    ]
)

# Generate random covariates for 100 individuals
using Random
Random.seed!(999)
covariates = [
    IndividualCovariates(
        Dict(:WT => 50.0 + 40.0 * rand(), :AGE => 20.0 + 60.0 * rand()),
        nothing
    )
    for _ in 1:100
]

# Population with IIV and covariates
pop = PopulationSpec(base_spec, iiv, nothing, cov_model, covariates)
result = simulate_population(pop, grid, solver)
```

---

## Time-Varying Covariates

Time-varying covariates represent dynamic patient characteristics that change during the simulation.

### Time Covariate Types

#### Step Time Covariate

Value held constant until the next knot point.

```julia
StepTimeCovariate()
```

#### Linear Time Covariate

Linear interpolation between knot points.

```julia
LinearTimeCovariate()
```

### Time Covariate Series

```julia
struct TimeCovariateSeries{K<:TimeCovariateKind}
    kind::K                    # StepTimeCovariate or LinearTimeCovariate
    times::Vector{Float64}     # Sorted, unique time points
    values::Vector{Float64}    # Values at each time point
end
```

### Time-Varying Covariates Container

```julia
struct TimeVaryingCovariates
    series::Dict{Symbol, Any}  # Symbol => TimeCovariateSeries
end
```

### Example: Enzyme Induction

```julia
using OpenPKPDCore

# Base model
base_params = OneCompIVBolusParams(5.0, 50.0)
doses = [DoseEvent(0.0, 100.0), DoseEvent(24.0, 100.0), DoseEvent(48.0, 100.0)]
base_spec = ModelSpec(OneCompIVBolus(), "enzyme_induction", base_params, doses)

# Covariate model: CL modified by enzyme activity
cov_model = CovariateModel(
    "enzyme_model",
    [CovariateEffect(LinearCovariate(), :CL, :ENZYME, 1.0, 1.0)]
)

# Time-varying enzyme activity (increases over time due to induction)
enzyme_series = TimeCovariateSeries(
    LinearTimeCovariate(),
    [0.0, 24.0, 48.0, 72.0],  # Times
    [1.0, 1.2, 1.5, 1.8]      # Enzyme activity (increasing)
)

tvc = TimeVaryingCovariates(Dict(:ENZYME => enzyme_series))

# Individual with time-varying covariate
covariates = [IndividualCovariates(Dict{Symbol,Float64}(), tvc)]

# Population (single individual with time-varying covariate)
pop = PopulationSpec(base_spec, nothing, nothing, cov_model, covariates)

grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

result = simulate_population(pop, grid, solver)
```

### Time-Varying Implementation Details

- Simulation is **segmented** at covariate boundary times
- Parameters are updated at each segment start
- State continuity is maintained across segments
- Both step and linear interpolation supported

---

## Population PKPD Simulation

Population simulations can include PD models.

### Example: Population PKPD with IIV

```julia
using OpenPKPDCore

# Base PK model
pk_params = OneCompIVBolusParams(5.0, 50.0)
pk_spec = ModelSpec(OneCompIVBolus(), "pop_pkpd", pk_params, [DoseEvent(0.0, 100.0)])

# PD model
pd_params = IndirectResponseTurnoverParams(10.0, 0.5, 20.0, 0.8, 1.0)
pd_spec = PDSpec(IndirectResponseTurnover(), "indirect", pd_params, :conc, :response)

# IIV on PK parameters only
iiv = IIVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.3, :V => 0.2),
    UInt64(12345),
    100
)

pop = PopulationSpec(pk_spec, iiv, nothing, nothing, [])

grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Include pd_spec in population simulation
result = simulate_population(pop, grid, solver; pd_spec=pd_spec)

# Summaries available for both PK and PD
println("Conc summary: ", result.summaries[:conc].mean)
println("Response summary: ", result.summaries[:response].mean)
```

---

## Segmented Simulation

For IOV and time-varying covariates, OpenPKPD uses segmented simulation.

### How It Works

1. Identify segment boundaries (occasion starts, covariate change times)
2. Run separate ODE solve for each segment
3. Use end state of segment N as initial state for segment N+1
4. Update parameters at segment boundaries
5. Merge results into single output aligned with `saveat` grid

### Direct Segmented API

```julia
# Manual segmented simulation (advanced use)
segment_starts = [0.0, 24.0, 48.0]
params_per_segment = [
    OneCompIVBolusParams(5.0, 50.0),   # Segment 1
    OneCompIVBolusParams(6.0, 50.0),   # Segment 2 (higher CL)
    OneCompIVBolusParams(7.0, 50.0)    # Segment 3 (even higher CL)
]

result = simulate_segmented_pk(
    base_spec,
    grid,
    solver,
    segment_starts,
    params_per_segment
)
```

---

## Complete Example

```julia
using OpenPKPDCore

# 1. Define base model
base_params = OneCompIVBolusParams(5.0, 50.0)
doses = [
    DoseEvent(0.0, 100.0),
    DoseEvent(24.0, 100.0),
    DoseEvent(48.0, 100.0)
]
base_spec = ModelSpec(OneCompIVBolus(), "complete_example", base_params, doses)

# 2. Define IIV
iiv = IIVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.3, :V => 0.2),
    UInt64(12345),
    100
)

# 3. Define IOV
iov = IOVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.15),
    UInt64(67890),
    OccasionDefinition(:dose_times)
)

# 4. Define covariate model
cov_model = CovariateModel(
    "full_cov_model",
    [
        CovariateEffect(PowerCovariate(), :CL, :WT, 0.75, 70.0),
        CovariateEffect(PowerCovariate(), :V, :WT, 1.0, 70.0)
    ]
)

# 5. Generate individual covariates
using Random
Random.seed!(999)
covariates = [
    IndividualCovariates(Dict(:WT => 50.0 + 50.0 * rand()), nothing)
    for _ in 1:100
]

# 6. Create population spec
pop = PopulationSpec(base_spec, iiv, iov, cov_model, covariates)

# 7. Run simulation
grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

result = simulate_population(pop, grid, solver)

# 8. Analyze results
summary = result.summaries[:conc]
println("Population mean Cmax: ", maximum(summary.mean))
println("Population 90% CI: ",
    minimum(summary.quantiles[0.05]), " - ",
    maximum(summary.quantiles[0.95]))
```
