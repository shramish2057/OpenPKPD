<p align="center">
  <h1 align="center">OpenPKPD</h1>
  <p align="center">
    <strong>Transparent, validated pharmacokinetics and pharmacodynamics modeling infrastructure</strong>
  </p>
  <p align="center">
    <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License: MIT"></a>
    <a href="https://julialang.org/"><img src="https://img.shields.io/badge/Julia-1.9+-purple.svg" alt="Julia"></a>
    <a href="https://python.org/"><img src="https://img.shields.io/badge/Python-3.9+-blue.svg" alt="Python"></a>
  </p>
  <p align="center">
    <a href="#features">Features</a> •
    <a href="#installation">Installation</a> •
    <a href="#quick-start">Quick Start</a> •
    <a href="#documentation">Documentation</a> •
    <a href="#contributing">Contributing</a>
  </p>
</p>

---

## Overview

**OpenPKPD** is a reference-grade PK/PD simulation platform designed for research, method development, and reproducible scientific computation. It emphasizes deterministic execution, transparent numerical semantics, and complete artifact serialization.

## Features

| Category | Features |
|----------|----------|
| **PK Models** | One-compartment IV bolus, oral first-order absorption |
| **PD Models** | Direct Emax, indirect response turnover |
| **Population** | IIV, IOV, static & time-varying covariates |
| **Sensitivity** | Single-subject and population-level analysis |
| **Interfaces** | Julia API, Python bindings, CLI |
| **Reproducibility** | Versioned artifacts with deterministic replay |

## Installation

### Julia (Core)

```bash
git clone https://github.com/openpkpd/openpkpd.git
cd openpkpd

# Install dependencies
julia --project=core/OpenPKPDCore -e 'using Pkg; Pkg.instantiate()'

# Verify installation
julia --project=core/OpenPKPDCore -e 'using OpenPKPDCore; println("v", OPENPKPD_VERSION)'
```

### Python (Optional)

```bash
cd python
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

### CLI

```bash
./bin/openpkpd version
```

## Quick Start

### Julia

```julia
using OpenPKPDCore

# Define model
spec = ModelSpec(
    OneCompIVBolus(),
    "example",
    OneCompIVBolusParams(5.0, 50.0),  # CL=5 L/h, V=50 L
    [DoseEvent(0.0, 100.0)]            # 100 mg at t=0
)

# Configure simulation
grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run and access results
result = simulate(spec, grid, solver)
println(result.observations[:conc])
```

### Python

```python
import openpkpd

openpkpd.init_julia()
result = openpkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)]
)
print(result["observations"]["conc"])
```

### CLI

```bash
# Check version
./bin/openpkpd version

# Replay an artifact
./bin/openpkpd replay --artifact validation/golden/pk_iv_bolus.json

# Validate all golden artifacts
./bin/openpkpd validate-golden
```

## Repository Structure

```
openpkpd/
├── core/OpenPKPDCore/    # Core Julia package
│   ├── src/              # Source code
│   └── test/             # Test suite
├── cli/OpenPKPDCLI/      # Command-line interface
├── python/               # Python bindings
├── docs/                 # Documentation (MkDocs)
│   └── examples/         # Executable examples
├── validation/           # Golden artifacts
│   ├── golden/           # Reference outputs
│   └── scripts/          # Validation runners
├── scripts/              # Development tools
└── bin/                  # CLI entry point
```

## Testing

```bash
# Julia unit tests
julia --project=core/OpenPKPDCore -e 'using Pkg; Pkg.test()'

# Golden artifact validation
./bin/openpkpd validate-golden

# Python tests
cd python && source .venv/bin/activate && pytest tests/

# Documentation build
mkdocs build --strict
```

## Semantic Versioning

OpenPKPD uses independent version numbers for numerical behavior:

| Version | Current | Scope |
|---------|---------|-------|
| Event Semantics | 1.0.0 | Dose handling |
| Solver Semantics | 1.0.0 | ODE solver behavior |
| Artifact Schema | 1.0.0 | JSON format |

Any change to numerical output requires a version bump.

## Documentation

Full documentation: [openpkpd.github.io/openpkpd](https://openpkpd.github.io/openpkpd/)

Build locally:
```bash
pip install -r docs/requirements.txt
mkdocs serve
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Citation

```bibtex
@software{openpkpd,
  title = {OpenPKPD: Transparent PK/PD Modeling Infrastructure},
  url = {https://github.com/openpkpd/openpkpd},
  version = {0.1.0}
}
```

---

<p align="center">
  <sub>Built for reproducible pharmacometric research</sub>
</p>
