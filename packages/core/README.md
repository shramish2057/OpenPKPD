# OpenPKPDCore

**Julia simulation engine for pharmacokinetics and pharmacodynamics**

OpenPKPDCore is the computational heart of OpenPKPD, providing high-performance ODE-based simulation, population modeling, parameter estimation, and model diagnostics.

## Installation

```julia
using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using OpenPKPDCore
```

## Quick Start

### Single Subject Simulation

```julia
using OpenPKPDCore

# Define model specification
spec = ModelSpec(
    OneCompIVBolus(),
    "example",
    OneCompIVBolusParams(5.0, 50.0),  # CL=5 L/h, V=50 L
    [DoseEvent(0.0, 100.0)]            # 100 mg bolus at t=0
)

# Define simulation grid and solver
grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run simulation
result = simulate(spec, grid, solver)

# Access results
println(result.t)                    # Time points
println(result.observations[:conc])  # Concentrations
println(result.states[:A_central])   # Amount in central compartment
```

### IV Infusion

```julia
# 100 mg infused over 1 hour (duration > 0)
doses = [DoseEvent(0.0, 100.0, 1.0)]  # time, amount, duration

spec = ModelSpec(
    OneCompIVBolus(),
    "infusion",
    OneCompIVBolusParams(5.0, 50.0),
    doses
)

result = simulate(spec, grid, solver)
```

### Population Simulation

```julia
# Define IIV specification
iiv = IIVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.3, :V => 0.2),  # Omega values (CV on log scale)
    UInt64(12345),                 # Seed for reproducibility
    100                            # Number of subjects
)

# Create population spec
pop_spec = PopulationSpec(spec, iiv, nothing, nothing, [])

# Run population simulation
pop_result = simulate_population(pop_spec, grid, solver)

# Access summaries
println(pop_result.summaries[:conc].mean)
println(pop_result.summaries[:conc].quantiles["0.05"])
println(pop_result.summaries[:conc].quantiles["0.95"])
```

## Available Models

### Pharmacokinetic Models

| Model | Parameters | Description |
|-------|------------|-------------|
| `OneCompIVBolus` | CL, V | One-compartment IV bolus/infusion |
| `OneCompOralFirstOrder` | Ka, CL, V | One-compartment oral absorption |
| `TwoCompIVBolus` | CL, V1, Q, V2 | Two-compartment IV |
| `TwoCompOral` | Ka, CL, V1, Q, V2 | Two-compartment oral |
| `ThreeCompIVBolus` | CL, V1, Q2, V2, Q3, V3 | Three-compartment IV |
| `TransitAbsorption` | N, Ktr, Ka, CL, V | Transit compartment absorption |
| `MichaelisMentenElimination` | Vmax, Km, V | Saturable elimination |

### Pharmacodynamic Models

| Model | Parameters | Description |
|-------|------------|-------------|
| `DirectEmax` | E0, Emax, EC50 | Direct effect model |
| `SigmoidEmax` | E0, Emax, EC50, gamma | Hill equation |
| `BiophaseEquilibration` | ke0, E0, Emax, EC50 | Effect compartment |
| `IndirectResponseTurnover` | Kin, Kout, R0, Imax, IC50 | Indirect response |

## Parameter Estimation (NLME)

OpenPKPDCore provides three estimation methods:

### FOCE-I (First-Order Conditional Estimation with Interaction)

```julia
config = EstimationConfig(
    FOCEIMethod(max_inner_iter=100, inner_tol=1e-6, centered=false),
    theta_init=[5.0, 50.0],      # Initial CL, V
    theta_lower=[0.1, 1.0],
    theta_upper=[100.0, 500.0],
    omega_init=diagm([0.09, 0.04]),
    omega_structure=:diagonal,
    sigma_init=ResidualErrorSpec(ProportionalError(), (sigma=0.1,), :conc, UInt64(0)),
    max_iter=500,
    compute_se=true
)

result = estimate(observed_data, model_spec, config; grid=grid, solver=solver)

println("Theta: ", result.theta)
println("Theta SE: ", result.theta_se)
println("OFV: ", result.ofv)
println("AIC: ", result.aic)
```

### SAEM (Stochastic Approximation EM)

```julia
config = EstimationConfig(
    SAEMMethod(n_burn=200, n_iter=300, n_chains=3),
    # ... same initialization
)

result = estimate(observed_data, model_spec, config; grid=grid, solver=solver)
```

### Laplacian Approximation

```julia
config = EstimationConfig(
    LaplacianMethod(),
    # ... same initialization
)
```

## Residual Error Models

```julia
# Additive error: Y = F + eps, eps ~ N(0, sigma^2)
error_spec = ResidualErrorSpec(AdditiveError(), (sigma=0.5,), :conc, seed)

# Proportional error: Y = F * (1 + eps), eps ~ N(0, sigma^2)
error_spec = ResidualErrorSpec(ProportionalError(), (sigma=0.1,), :conc, seed)

# Combined error: Y = F + F*eps1 + eps2
error_spec = ResidualErrorSpec(CombinedError(), (sigma_add=0.5, sigma_prop=0.1), :conc, seed)

# Exponential error: Y = F * exp(eps)
error_spec = ResidualErrorSpec(ExponentialError(), (sigma=0.1,), :conc, seed)
```

## Visual Predictive Check (VPC)

```julia
vpc_config = VPCConfig(
    pi_levels=[0.05, 0.50, 0.95],
    ci_level=0.95,
    binning=QuantileBinning(n_bins=10),
    prediction_corrected=false,
    stratify_by=Symbol[],
    lloq=nothing,
    n_bootstrap=500,
    seed=UInt64(12345)
)

vpc_result = compute_vpc(observed_data, pop_spec, grid, solver; config=vpc_config)
```

## CDISC Data Import

```julia
# Read CDISC domains
dataset = read_cdisc_csv("pc.csv", "ex.csv", "dm.csv")

# Or read SAS transport files
dataset = read_cdisc_xpt("pc.xpt", "ex.xpt", "dm.xpt")

# Validate dataset
warnings = validate_cdisc_dataset(dataset)

# Convert to OpenPKPD format
pop_spec, observed = cdisc_to_population(dataset, model_spec)
```

## NONMEM/Monolix Import

```julia
# Parse NONMEM control file
nmctl = parse_nonmem_control("run001.ctl")
model_spec, pop_spec, mapping = convert_nonmem_to_openpkpd(nmctl)

# Parse Monolix project
mlx = parse_monolix_project("project.mlxtran")
model_spec, pop_spec, mapping = convert_monolix_to_openpkpd(mlx)
```

## Covariate Models

```julia
# Static covariates
cov_model = CovariateModel([
    CovariateEffect(:CL, :WT, PowerCovariate(0.75, 70.0)),    # Allometric
    CovariateEffect(:V, :WT, PowerCovariate(1.0, 70.0)),
    CovariateEffect(:CL, :CRCL, LinearCovariate(0.5, 100.0)), # Renal function
])

covariates = [
    SubjectCovariates("001", Dict(:WT => 80.0, :CRCL => 90.0)),
    SubjectCovariates("002", Dict(:WT => 65.0, :CRCL => 120.0)),
]

pop_spec = PopulationSpec(base_spec, iiv, nothing, cov_model, covariates)

# Time-varying covariates
tv_cov = TimeVaryingCovariate(:WT, [0.0, 24.0, 48.0], [70.0, 71.0, 72.0], :linear)
```

## Sensitivity Analysis

```julia
# Single-subject sensitivity
result = run_sensitivity_single(
    spec, grid, solver,
    :CL, 0.1,  # Parameter and perturbation fraction
    :conc      # Observation to analyze
)

println(result.metrics.max_abs_delta)
println(result.metrics.max_rel_delta)

# Population sensitivity
pop_result = run_sensitivity_population(
    pop_spec, grid, solver,
    :CL, 0.1, :conc
)
```

## Serialization

```julia
# Serialize to artifact
artifact = serialize_single_artifact(spec, grid, solver, result)
write_artifact("simulation.json", artifact)

# Deserialize and replay
loaded = read_artifact("simulation.json")
replayed = replay_artifact(loaded)

# Population artifacts
pop_artifact = serialize_population_artifact(pop_spec, grid, solver, pop_result)
```

## Key Types

| Type | Description |
|------|-------------|
| `ModelSpec` | Complete model specification (kind, params, doses) |
| `SimGrid` | Time domain (t0, t1, saveat points) |
| `SolverSpec` | ODE solver configuration |
| `SimResult` | Single simulation output |
| `PopulationSpec` | Population model (IIV, IOV, covariates) |
| `PopulationResult` | Population simulation output with summaries |
| `DoseEvent` | Dose administration (time, amount, duration) |
| `EstimationConfig` | NLME estimation settings |
| `EstimationResult` | Fitted parameters with diagnostics |
| `VPCConfig` | VPC computation settings |
| `VPCResult` | VPC bins and statistics |

## Module Structure

```
OpenPKPDCore/
├── src/
│   ├── models/           # PK/PD model implementations
│   ├── engine/           # solve.jl, population.jl, infusion.jl
│   ├── specs/            # Type definitions
│   ├── estimation/       # FOCE, SAEM, Laplacian
│   ├── import/           # NONMEM/Monolix parsers
│   ├── data/             # CDISC handling
│   ├── analysis/         # VPC, sensitivity
│   └── serialization/    # JSON I/O
└── test/                 # Test suite
```

## Testing

```bash
julia --project=packages/core -e 'using Pkg; Pkg.test()'
```

## License

MIT License - see repository root for details.
