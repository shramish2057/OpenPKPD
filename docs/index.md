# OpenPKPD Documentation

**OpenPKPD** is a transparent, validated pharmacokinetics and pharmacodynamics (PK/PD) modeling infrastructure built for reproducibility and scientific rigor.

## Key Features

- **Validated PK/PD Models**: One-compartment IV bolus, oral first-order absorption, direct Emax, and indirect response models
- **Population Simulation**: Inter-individual variability (IIV), inter-occasion variability (IOV), and covariate effects
- **Time-Varying Covariates**: Step and linear interpolation for dynamic patient characteristics
- **Sensitivity Analysis**: Single-subject and population-level parameter sensitivity
- **Reproducible Artifacts**: JSON-serialized execution artifacts with semantic versioning
- **Multi-Language Support**: Julia core with Python bindings and CLI interface

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/openpkpd/openpkpd.git
cd openpkpd

# Install Julia dependencies
julia --project=core/OpenPKPDCore -e 'using Pkg; Pkg.instantiate()'

# (Optional) Install Python bindings
cd python
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

### Your First Simulation (Julia)

```julia
using OpenPKPDCore

# Define a one-compartment IV bolus model
spec = ModelSpec(
    OneCompIVBolus(),
    "quickstart",
    OneCompIVBolusParams(5.0, 50.0),  # CL=5 L/h, V=50 L
    [DoseEvent(0.0, 100.0)]            # 100 mg at t=0
)

# Define simulation grid and solver
grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
solver = SolverSpec(:Tsit5, 1e-6, 1e-9, 10000)

# Run simulation
result = simulate(spec, grid, solver)

# Access results
println("Time points: ", result.t)
println("Concentrations: ", result.observations[:conc])
```

### Your First Simulation (Python)

```python
import openpkpd

# Initialize Julia bridge
openpkpd.init_julia()

# Run IV bolus simulation
result = openpkpd.simulate_pk_iv_bolus(
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=24.0,
    saveat=[float(t) for t in range(25)]
)

print("Concentrations:", result["observations"]["conc"])
```

### Using the CLI

```bash
# Check version
./bin/openpkpd version

# Replay an artifact
./bin/openpkpd replay --artifact validation/golden/pk_iv_bolus.json

# Validate golden artifacts
./bin/openpkpd validate-golden
```

## Core Concepts

### Models

OpenPKPD provides a comprehensive library of validated pharmacokinetic and pharmacodynamic models:

#### Pharmacokinetic Models

| Model | Description | Key Use Cases |
|-------|-------------|---------------|
| `OneCompIVBolus` | One-compartment with IV bolus dosing | Simple IV drugs, initial PK characterization |
| `OneCompOralFirstOrder` | One-compartment with first-order oral absorption | Immediate-release oral formulations |
| `TwoCompIVBolus` | Two-compartment with IV bolus dosing | Distribution phase modeling, antibodies |
| `TwoCompOral` | Two-compartment with first-order oral absorption | Oral drugs with tissue distribution |
| `ThreeCompIVBolus` | Three-compartment with IV bolus dosing | Deep tissue distribution, long half-life drugs |
| `TransitAbsorption` | Transit compartment absorption model | Delayed/complex oral absorption, GI transit |
| `MichaelisMentenElimination` | Saturable (nonlinear) elimination | High-dose drugs, biologics, enzyme saturation |

#### Pharmacodynamic Models

| Model | Description | Key Use Cases |
|-------|-------------|---------------|
| `DirectEmax` | Direct effect Emax model | Immediate drug effects, receptor binding |
| `SigmoidEmax` | Sigmoid Emax with Hill coefficient | Steep dose-response, cooperative binding |
| `BiophaseEquilibration` | Effect compartment model | CNS effects, delayed equilibration |
| `IndirectResponseTurnover` | Indirect response turnover model | Enzyme/receptor modulation, biomarkers |

### Non-Compartmental Analysis (NCA)

FDA/EMA-compliant NCA for exposure assessment:

```python
from openpkpd.nca import run_nca

result = run_nca(times, conc, dose)
print(f"Cmax: {result.cmax}, t½: {result.t_half}")
```

### Clinical Trial Simulation

Complete trial simulation capabilities:

```python
from openpkpd import trial

design = trial.parallel_design(2)
regimen = trial.dosing_qd(100.0, 28)
pop = trial.generate_virtual_population(100)
```

### Visualization

Dual-backend (Matplotlib/Plotly) visualization:

```python
from openpkpd import viz

viz.set_backend("matplotlib")
viz.plot_conc_time(result)
viz.plot_vpc(pop_result, observed_data)
```

### Population Variability

Model variability at multiple levels:

- **IIV (Inter-Individual Variability)**: Log-normal distribution of parameters across individuals
- **IOV (Inter-Occasion Variability)**: Parameter variation between dosing occasions
- **Covariates**: Linear, power, or exponential effects of patient characteristics

### Artifacts & Reproducibility

Every simulation can be serialized to a JSON artifact containing:

- Complete model specification
- Solver settings and grid
- Results with full precision
- Semantic version fingerprint

Artifacts can be replayed to verify reproducibility across versions.

## Documentation Contents

| Section | Description |
|---------|-------------|
| [Models](models.md) | Complete PK and PD model reference |
| [NCA](nca.md) | Non-compartmental analysis (FDA/EMA compliant) |
| [Trial Simulation](trial.md) | Clinical trial design and simulation |
| [Visualization](visualization.md) | Matplotlib and Plotly plotting |
| [Population Simulation](population.md) | IIV, IOV, and covariate modeling |
| [Sensitivity Analysis](sensitivity.md) | Parameter sensitivity methods |
| [Architecture](architecture.md) | System design and boundaries |
| [Semantics](semantics.md) | Versioning and numerical semantics |
| [Reproducibility](reproducibility.md) | Artifacts and validation |
| [CLI Reference](cli.md) | Command-line interface |
| [Python Bindings](python.md) | Python API reference |

## Version Information

```
OpenPKPD Version: 0.1.0
Event Semantics: 1.0.0
Solver Semantics: 1.0.0
Artifact Schema: 1.0.0
```

## License

OpenPKPD is open source software. See the repository for license details.
