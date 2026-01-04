# OpenPKPD Documentation

**OpenPKPD** is a transparent, validated pharmacokinetics and pharmacodynamics (PK/PD) modeling infrastructure built for reproducibility and scientific rigor.

## Key Features

- **Comprehensive PK/PD Models**: One/Two/Three-compartment IV & oral, transit absorption, Michaelis-Menten, direct Emax, sigmoid Emax, biophase equilibration, indirect response
- **IV Infusion Support**: Zero-order infusion with duration specification, overlapping infusions
- **Population Simulation**: Inter-individual variability (IIV), inter-occasion variability (IOV), and covariate effects
- **Parameter Estimation (NLME)**: FOCE-I, SAEM, and Laplacian estimation methods with standard error computation
- **Non-Compartmental Analysis (NCA)**: FDA/EMA-compliant exposure metrics with bioequivalence analysis
- **Visual Predictive Checks (VPC)**: Standard VPC, prediction-corrected VPC (pcVPC), stratification, bootstrap CIs
- **Clinical Trial Simulation**: Parallel, crossover, dose-escalation, bioequivalence designs with power analysis
- **Model Import**: NONMEM (.ctl) and Monolix (.mlxtran) model file parsing
- **Data Import**: CDISC/SDTM format support (PC, EX, DM domains) in CSV and XPT formats
- **Residual Error Models**: Additive, proportional, combined, and exponential error
- **Sensitivity Analysis**: Single-subject and population-level parameter sensitivity
- **Reproducible Artifacts**: JSON-serialized execution artifacts with semantic versioning
- **Multi-Language Support**: Julia core with Python bindings and CLI interface
- **Professional Visualization**: Matplotlib and Plotly backends with estimation diagnostics

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/openpkpd/openpkpd.git
cd openpkpd

# Install Julia dependencies
julia --project=packages/core -e 'using Pkg; Pkg.instantiate()'

# (Optional) Install Python bindings
cd packages/python
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[all]"
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
./packages/cli/bin/openpkpd version

# Run simulation from spec
./packages/cli/bin/openpkpd simulate --spec simulation.json --out result.json

# Run NCA analysis
./packages/cli/bin/openpkpd nca --spec nca_spec.json --out nca_result.json

# Run parameter estimation
./packages/cli/bin/openpkpd estimate --spec estimate_spec.json --out fit.json

# Compute VPC
./packages/cli/bin/openpkpd vpc --spec vpc_spec.json --out vpc_result.json

# Run clinical trial simulation
./packages/cli/bin/openpkpd trial --spec trial_spec.json --out trial_result.json

# Import NONMEM model
./packages/cli/bin/openpkpd import --input run001.ctl --format nonmem --out model.json

# Replay an artifact
./packages/cli/bin/openpkpd replay --artifact validation/golden/pk_iv_bolus.json

# Validate golden artifacts
./packages/cli/bin/openpkpd validate-golden
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
| [Parameter Estimation](estimation.md) | NLME estimation (FOCE-I, SAEM, Laplacian) |
| [Population Simulation](population.md) | IIV, IOV, and covariate modeling |
| [NCA](nca.md) | Non-compartmental analysis (FDA/EMA compliant) |
| [VPC](vpc.md) | Visual Predictive Checks |
| [Trial Simulation](trial.md) | Clinical trial design and simulation |
| [Data Import](data.md) | CDISC/SDTM data format support |
| [Model Import](import.md) | NONMEM and Monolix model parsing |
| [Visualization](visualization.md) | Matplotlib and Plotly plotting |
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
