# Visual Predictive Checks (VPC)

Visual Predictive Checks (VPC) are a model evaluation tool that compares observed data to model predictions. OpenPKPD supports standard VPC, prediction-corrected VPC (pcVPC), stratification, and bootstrap confidence intervals.

## Overview

A VPC compares:
- **Observed data percentiles** (e.g., 5th, 50th, 95th) computed from actual observations
- **Simulated percentiles** from model-based population simulations

If the model adequately describes the data, observed percentiles should fall within the confidence intervals of simulated percentiles.

## Basic Usage

### Julia

```julia
using OpenPKPDCore

# Define VPC configuration
vpc_config = VPCConfig(
    pi_levels=[0.05, 0.50, 0.95],   # Prediction interval levels
    ci_level=0.95,                   # Confidence interval for percentiles
    binning=QuantileBinning(n_bins=10),
    prediction_corrected=false,
    stratify_by=Symbol[],
    lloq=nothing,                    # Lower limit of quantification
    n_bootstrap=500,
    seed=UInt64(12345)
)

# Compute VPC
vpc_result = compute_vpc(
    observed_data,    # ObservedData struct
    pop_spec,         # PopulationSpec with model
    grid,             # SimGrid
    solver;           # SolverSpec
    config=vpc_config,
    error_spec=residual_error_spec  # Optional residual error
)
```

### Python

```python
from openpkpd.analysis import compute_vpc, VPCConfig

config = VPCConfig(
    pi_levels=[0.05, 0.50, 0.95],
    ci_level=0.95,
    n_bins=10,
    binning_strategy="quantile",
    prediction_corrected=False,
    n_bootstrap=500,
    seed=12345
)

vpc_result = compute_vpc(
    observed_data,
    pop_spec,
    grid,
    solver,
    config=config
)

# Visualize
from openpkpd import viz
viz.plot_vpc(vpc_result)
```

## Configuration Options

### Prediction Intervals

```python
config = VPCConfig(
    pi_levels=[0.05, 0.50, 0.95]  # 5th, 50th, 95th percentiles
)
```

Common choices:
- `[0.05, 0.50, 0.95]` - Standard 90% prediction interval
- `[0.10, 0.50, 0.90]` - 80% prediction interval
- `[0.025, 0.50, 0.975]` - 95% prediction interval

### Confidence Intervals

```python
config = VPCConfig(
    ci_level=0.95,      # 95% CI for percentiles
    n_bootstrap=500     # Number of bootstrap samples
)
```

Higher `n_bootstrap` gives more stable CIs but increases computation time.

### Binning Strategies

#### Quantile Binning (Default)

Bins have equal numbers of observations.

```julia
binning = QuantileBinning(n_bins=10)
```

#### Equal Width Binning

Bins have equal time width.

```julia
binning = EqualWidthBinning(n_bins=10)
```

#### Custom Bins

```python
config = VPCConfig(
    binning_strategy="custom",
    bin_edges=[0, 1, 2, 4, 8, 12, 24]
)
```

### Lower Limit of Quantification (LLOQ)

Handle BLQ (below limit of quantification) data:

```python
config = VPCConfig(
    lloq=0.1  # Values below 0.1 are BLQ
)
```

## Prediction-Corrected VPC (pcVPC)

pcVPC normalizes observations and predictions to account for design differences (e.g., different doses, times). This is especially useful for sparse data or studies with varying protocols.

### How pcVPC Works

1. For each bin, compute the typical prediction (median of simulated data)
2. Normalize observations: `obs_corrected = obs * median_pred / pred`
3. Compute percentiles on corrected data

### Usage

```julia
vpc_config = VPCConfig(
    pi_levels=[0.05, 0.50, 0.95],
    ci_level=0.95,
    binning=QuantileBinning(n_bins=10),
    prediction_corrected=true,  # Enable pcVPC
    n_bootstrap=500,
    seed=UInt64(12345)
)
```

```python
config = VPCConfig(
    prediction_corrected=True
)
```

## Stratification

Stratify VPC by categorical covariates (e.g., study, dose group, formulation):

```julia
vpc_config = VPCConfig(
    stratify_by=[:STUDY, :DOSE_GRP],
    # ... other settings
)
```

```python
config = VPCConfig(
    stratify_by=["STUDY", "DOSE_GRP"]
)
```

This produces separate VPC panels for each stratum.

## VPC Result Structure

The `VPCResult` contains:

| Field | Description |
|-------|-------------|
| `config` | VPC configuration used |
| `bins` | List of VPCBin objects |
| `observed_data` | Original observed data |
| `n_simulations` | Number of simulated datasets |

Each `VPCBin` contains:

| Field | Description |
|-------|-------------|
| `time_lower` | Bin lower time bound |
| `time_upper` | Bin upper time bound |
| `time_median` | Median time in bin |
| `obs_percentiles` | Observed data percentiles |
| `sim_percentiles` | Simulated percentiles (median) |
| `sim_ci_lower` | Lower CI for simulated percentiles |
| `sim_ci_upper` | Upper CI for simulated percentiles |
| `n_obs` | Number of observations in bin |

## Visualization

### Standard VPC Plot

```python
from openpkpd import viz

fig = viz.plot_vpc(
    vpc_result,
    log_scale=False,
    title="Visual Predictive Check",
    xlabel="Time (h)",
    ylabel="Concentration (mg/L)",
    figsize=(10, 6)
)
```

### VPC from Dict (Legacy)

```python
# For VPC results returned as dicts
fig = viz.plot_vpc(vpc_dict)
```

### Customization

```python
fig = viz.plot_vpc(
    vpc_result,
    show_observed=True,      # Show observed points
    show_pi=True,            # Show prediction intervals
    show_ci=True,            # Show confidence intervals
    obs_alpha=0.3,           # Observed point transparency
    ci_alpha=0.2,            # CI ribbon transparency
    colors={
        "obs": "black",
        "median": "blue",
        "pi": "lightblue",
        "ci": "lightgray"
    }
)
```

## CLI Usage

```bash
./packages/cli/bin/openpkpd vpc --spec vpc_spec.json --out vpc_result.json
```

**Specification Format:**
```json
{
  "observed_data": "observed.csv",
  "population_spec": {
    "base_model": {
      "kind": "OneCompIVBolus",
      "params": {"CL": 5.0, "V": 50.0},
      "doses": [{"time": 0.0, "amount": 100.0}]
    },
    "iiv": {
      "kind": "LogNormalIIV",
      "omegas": {"CL": 0.3, "V": 0.2},
      "n": 1000
    }
  },
  "config": {
    "pi_levels": [0.05, 0.50, 0.95],
    "ci_level": 0.95,
    "binning": {
      "strategy": "quantile",
      "n_bins": 10
    },
    "prediction_corrected": false,
    "lloq": null,
    "n_bootstrap": 500,
    "seed": 12345
  },
  "error_spec": {
    "kind": "proportional",
    "sigma": 0.1
  },
  "grid": {
    "t0": 0.0,
    "t1": 24.0,
    "saveat": [0, 0.5, 1, 2, 4, 8, 12, 24]
  },
  "solver": {
    "alg": "Tsit5",
    "reltol": 1e-8,
    "abstol": 1e-10
  }
}
```

## Interpreting VPC Results

### Good Model Fit

- Observed percentiles (lines) fall within simulated CIs (ribbons)
- Median observed data follows median simulated data
- No systematic trends in deviations

### Signs of Misspecification

| Pattern | Possible Issue |
|---------|----------------|
| Observed median below simulated | Model overpredicts |
| Observed PI outside simulated PI | Variability underestimated |
| Trend in residuals | Structural model issue |
| Time-dependent bias | Absorption/elimination misspecified |

## Tips and Best Practices

1. **Sample Size**: Use at least 500-1000 simulations for stable percentiles

2. **Binning**: Use quantile binning when observation times vary; use equal-width when times are similar

3. **pcVPC**: Use for studies with multiple doses or varying protocols

4. **Stratification**: Stratify when you suspect different subpopulations behave differently

5. **LLOQ**: Always specify LLOQ if your data has BLQ observations

6. **Residual Error**: Include residual error for realistic variability

## See Also

- [Parameter Estimation](estimation.md) - Fitting models to data
- [Population Simulation](population.md) - Running population simulations
- [Visualization](visualization.md) - Plotting functions
