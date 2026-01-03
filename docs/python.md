# Python Bindings

OpenPKPD provides Python bindings through the `openpkpd` package, enabling seamless integration with Python data science workflows.

## Installation

```bash
cd python
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

**Requirements**:

- Python 3.9+
- Julia 1.9+ (installed and in PATH)
- juliacall package (installed automatically)

---

## Quick Start

```python
import openpkpd

# Initialize Julia (required once per session)
openpkpd.init_julia()

# Check version
print(openpkpd.version())

# Run a simulation
result = openpkpd.simulate_pk_iv_bolus(
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=24.0,
    saveat=[float(t) for t in range(25)]
)

print("Times:", result["t"])
print("Concentrations:", result["observations"]["conc"])
```

---

## API Reference

### Initialization

#### `init_julia(repo_root=None)`

Initialize the Julia runtime and load OpenPKPDCore.

```python
import openpkpd

# Auto-detect repository root
openpkpd.init_julia()

# Or specify explicitly
openpkpd.init_julia("/path/to/openpkpd")
```

**Notes**:

- Safe to call multiple times (no-op after first call)
- Must be called before any simulation functions
- Automatically activates the OpenPKPDCore Julia project

#### `version()`

Get the OpenPKPD version string.

```python
print(openpkpd.version())  # "0.1.0"
```

---

### Single Simulation

#### `simulate_pk_iv_bolus(...)`

Run a one-compartment IV bolus simulation.

```python
result = openpkpd.simulate_pk_iv_bolus(
    cl=5.0,              # Clearance (L/h)
    v=50.0,              # Volume of distribution (L)
    doses=[              # List of dose events
        {"time": 0.0, "amount": 100.0},
        {"time": 12.0, "amount": 50.0}
    ],
    t0=0.0,              # Start time
    t1=24.0,             # End time
    saveat=[0.0, 1.0, 2.0, ...],  # Output times
    alg="Tsit5",         # ODE solver (optional)
    reltol=1e-10,        # Relative tolerance (optional)
    abstol=1e-12,        # Absolute tolerance (optional)
    maxiters=10_000_000  # Max iterations (optional)
)
```

**Parameters**:

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `cl` | float | Yes | Clearance |
| `v` | float | Yes | Volume of distribution |
| `doses` | list[dict] | Yes | Dose events |
| `t0` | float | Yes | Start time |
| `t1` | float | Yes | End time |
| `saveat` | list[float] | Yes | Output time points |
| `alg` | str | No | ODE solver (default: "Tsit5") |
| `reltol` | float | No | Relative tolerance |
| `abstol` | float | No | Absolute tolerance |
| `maxiters` | int | No | Maximum iterations |

**Returns**: `dict` with keys:

- `t`: List of time points
- `states`: Dict of state trajectories (`A_central`)
- `observations`: Dict of observables (`conc`)
- `metadata`: Dict of simulation metadata

---

#### `simulate_pk_oral_first_order(...)`

Run a one-compartment oral first-order absorption simulation.

```python
result = openpkpd.simulate_pk_oral_first_order(
    ka=1.5,              # Absorption rate constant (1/h)
    cl=5.0,              # Clearance (L/h)
    v=50.0,              # Volume of distribution (L)
    doses=[{"time": 0.0, "amount": 200.0}],
    t0=0.0,
    t1=24.0,
    saveat=[float(t) for t in range(25)]
)
```

**Parameters**:

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `ka` | float | Yes | Absorption rate constant |
| `cl` | float | Yes | Clearance |
| `v` | float | Yes | Volume of distribution |
| `doses` | list[dict] | Yes | Dose events |
| `t0` | float | Yes | Start time |
| `t1` | float | Yes | End time |
| `saveat` | list[float] | Yes | Output time points |

**Returns**: Same structure as `simulate_pk_iv_bolus`

---

### Population Simulation

#### `simulate_population_iv_bolus(...)`

Run a population IV bolus simulation with inter-individual variability.

```python
result = openpkpd.simulate_population_iv_bolus(
    cl=5.0,              # Typical CL
    v=50.0,              # Typical V
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=24.0,
    saveat=[float(t) for t in range(25)],
    n=100,               # Number of individuals
    seed=12345,          # RNG seed for reproducibility
    omegas={             # IIV standard deviations
        "CL": 0.3,
        "V": 0.2
    }
)
```

**Parameters**:

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `cl` | float | Yes | Typical clearance |
| `v` | float | Yes | Typical volume |
| `doses` | list[dict] | Yes | Dose events |
| `t0` | float | Yes | Start time |
| `t1` | float | Yes | End time |
| `saveat` | list[float] | Yes | Output time points |
| `n` | int | Yes | Population size |
| `seed` | int | Yes | Random seed |
| `omegas` | dict | Yes | IIV omegas (parameter: SD) |

**Returns**: `dict` with keys:

- `individuals`: List of individual SimResult dicts
- `params`: List of realized parameter dicts
- `summaries`: Dict of PopulationSummary dicts
- `metadata`: Dict of population metadata

**Example - Accessing Population Results**:

```python
result = openpkpd.simulate_population_iv_bolus(...)

# Individual results
for i, individual in enumerate(result["individuals"]):
    print(f"Individual {i}: Cmax = {max(individual['observations']['conc'])}")

# Realized parameters
for i, params in enumerate(result["params"]):
    print(f"Individual {i}: CL={params['CL']:.2f}, V={params['V']:.2f}")

# Population summary
summary = result["summaries"]["conc"]
print("Mean concentrations:", summary["mean"])
print("5th percentile:", summary["quantiles"]["0.05"])
print("95th percentile:", summary["quantiles"]["0.95"])
```

---

### Artifact Operations

#### `write_single_artifact(path, *, model, grid, solver=None)`

Write a single simulation to a JSON artifact file.

```python
openpkpd.write_single_artifact(
    "my_simulation.json",
    model={
        "kind": "OneCompIVBolus",
        "params": {"CL": 5.0, "V": 50.0},
        "doses": [{"time": 0.0, "amount": 100.0}]
    },
    grid={
        "t0": 0.0,
        "t1": 24.0,
        "saveat": [float(t) for t in range(25)]
    },
    solver={  # Optional
        "alg": "Tsit5",
        "reltol": 1e-10,
        "abstol": 1e-12,
        "maxiters": 10_000_000
    }
)
```

**Model Kinds**:

- `"OneCompIVBolus"`: Requires `params: {"CL": float, "V": float}`
- `"OneCompOralFirstOrder"`: Requires `params: {"KA": float, "CL": float, "V": float}`

---

#### `replay_artifact(path)`

Replay any artifact and return Python dict.

```python
# Replay single artifact
result = openpkpd.replay_artifact("validation/golden/pk_iv_bolus.json")
print("Concentrations:", result["observations"]["conc"])

# Replay population artifact
result = openpkpd.replay_artifact("validation/golden/population_iv_bolus.json")
print("Number of individuals:", len(result["individuals"]))

# Replay sensitivity artifact
result = openpkpd.replay_artifact("validation/golden/sensitivity_single.json")
print("Metrics:", result["metrics"])
```

**Supported Artifact Types**:

| Type | Return Structure |
|------|------------------|
| Single | `{"t", "states", "observations", "metadata"}` |
| Population | `{"individuals", "params", "summaries", "metadata"}` |
| Sensitivity (single) | `{"plan", "observation", "base_series", "pert_series", "metrics", "metadata"}` |
| Sensitivity (population) | `{"plan", "observation", "probs", "base_mean", "pert_mean", "metrics_mean", "metadata"}` |

---

## Integration Examples

### With NumPy

```python
import numpy as np
import openpkpd

openpkpd.init_julia()

result = openpkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=list(np.linspace(0, 24, 49))
)

t = np.array(result["t"])
conc = np.array(result["observations"]["conc"])

# Calculate metrics
print(f"Cmax: {np.max(conc):.2f}")
print(f"AUC: {np.trapz(conc, t):.2f}")
```

### With Pandas

```python
import pandas as pd
import openpkpd

openpkpd.init_julia()

result = openpkpd.simulate_population_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)],
    n=100, seed=12345,
    omegas={"CL": 0.3, "V": 0.2}
)

# Create parameter DataFrame
params_df = pd.DataFrame(result["params"])
print(params_df.describe())

# Create concentration DataFrame
data = []
for i, ind in enumerate(result["individuals"]):
    for t, c in zip(ind["t"], ind["observations"]["conc"]):
        data.append({"id": i, "time": t, "conc": c})

conc_df = pd.DataFrame(data)
print(conc_df.groupby("time")["conc"].describe())
```

### With Matplotlib

```python
import matplotlib.pyplot as plt
import openpkpd

openpkpd.init_julia()

result = openpkpd.simulate_population_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) * 0.5 for t in range(49)],
    n=100, seed=12345,
    omegas={"CL": 0.3, "V": 0.2}
)

# Plot individual trajectories
plt.figure(figsize=(10, 6))
for ind in result["individuals"][:20]:  # First 20 individuals
    plt.plot(ind["t"], ind["observations"]["conc"], alpha=0.3, color="blue")

# Plot population summary
summary = result["summaries"]["conc"]
t = result["individuals"][0]["t"]
plt.plot(t, summary["mean"], "k-", linewidth=2, label="Mean")
plt.fill_between(t, summary["quantiles"]["0.05"], summary["quantiles"]["0.95"],
                 alpha=0.2, label="90% CI")

plt.xlabel("Time (h)")
plt.ylabel("Concentration (mg/L)")
plt.legend()
plt.title("Population PK Simulation")
plt.savefig("population_pk.png")
```

---

## Error Handling

```python
import openpkpd

try:
    openpkpd.init_julia()
    result = openpkpd.simulate_pk_iv_bolus(
        cl=-5.0,  # Invalid: negative clearance
        v=50.0,
        doses=[{"time": 0.0, "amount": 100.0}],
        t0=0.0, t1=24.0,
        saveat=[0.0, 1.0]
    )
except Exception as e:
    print(f"Simulation error: {e}")
```

---

## Performance Tips

1. **Initialize Once**: Call `init_julia()` once at the start, not per simulation
2. **Batch Simulations**: Julia JIT compilation makes subsequent calls faster
3. **Use Appropriate Tolerances**: Higher tolerance = faster, less accurate
4. **Parallel Python**: For many independent simulations, consider `multiprocessing`

```python
# Good: Initialize once
import openpkpd
openpkpd.init_julia()

for params in parameter_sets:
    result = openpkpd.simulate_pk_iv_bolus(...)

# Bad: Don't call init_julia repeatedly (though it's safe, it's unnecessary)
```

---

## Troubleshooting

### Julia Not Found

```
RuntimeError: Julia not found
```

**Solution**: Install Julia and ensure it's in your PATH:

```bash
# macOS
brew install julia

# Ubuntu
sudo apt install julia

# Verify
julia --version
```

### First Call Slow

The first simulation call may take 30-60 seconds due to Julia compilation.

**Solution**: This is normal. Subsequent calls will be fast.

### Memory Issues

```
OutOfMemoryError
```

**Solution**: Reduce population size or output grid density:

```python
# Use sparser output grid
saveat = [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]  # Instead of every 0.1h
```

---

## Module Contents

```python
# All exported functions
from openpkpd import (
    init_julia,
    version,
    simulate_pk_iv_bolus,
    simulate_pk_oral_first_order,
    simulate_population_iv_bolus,
    write_single_artifact,
    replay_artifact,
)
```
