# Python Bindings

OpenPKPD provides comprehensive Python bindings through the `openpkpd` package, enabling seamless integration with Python data science workflows.

## Installation

```bash
cd python
python3 -m venv .venv
source .venv/bin/activate
pip install -e .

# With visualization support
pip install -e ".[viz]"

# With all optional dependencies
pip install -e ".[all]"
```

**Requirements**:

- Python 3.10+
- Julia 1.10+ (installed and in PATH)
- juliacall package (installed automatically)

---

## Quick Start

```python
import openpkpd

# Initialize Julia (required once per session)
openpkpd.init_julia()

# Check version
print(openpkpd.version())

# Run a simple simulation
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

## Module Overview

| Module | Description |
|--------|-------------|
| `openpkpd` | Core simulation functions |
| `openpkpd.nca` | Non-compartmental analysis |
| `openpkpd.trial` | Clinical trial simulation |
| `openpkpd.viz` | Visualization (matplotlib/plotly) |

---

## Core Simulation Functions

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

### One-Compartment PK Models

#### `simulate_pk_iv_bolus(...)`

One-compartment IV bolus simulation.

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
    saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
)
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `cl` | float | Yes | Clearance (L/h) |
| `v` | float | Yes | Volume of distribution (L) |
| `doses` | list[dict] | Yes | Dose events with `time` and `amount` |
| `t0` | float | Yes | Start time |
| `t1` | float | Yes | End time |
| `saveat` | list[float] | Yes | Output time points |

#### `simulate_pk_oral_first_order(...)`

One-compartment oral first-order absorption.

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

---

### Two-Compartment PK Models

#### `simulate_pk_twocomp_iv_bolus(...)`

Two-compartment IV bolus simulation.

```python
result = openpkpd.simulate_pk_twocomp_iv_bolus(
    cl=5.0,              # Clearance (L/h)
    v1=50.0,             # Central volume (L)
    q=2.0,               # Inter-compartmental clearance (L/h)
    v2=100.0,            # Peripheral volume (L)
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=24.0,
    saveat=[0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
)
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `cl` | float | Yes | Clearance (L/h) |
| `v1` | float | Yes | Central volume (L) |
| `q` | float | Yes | Inter-compartmental clearance (L/h) |
| `v2` | float | Yes | Peripheral volume (L) |

#### `simulate_pk_twocomp_oral(...)`

Two-compartment oral absorption.

```python
result = openpkpd.simulate_pk_twocomp_oral(
    ka=1.5,              # Absorption rate constant (1/h)
    cl=5.0,              # Clearance (L/h)
    v1=50.0,             # Central volume (L)
    q=2.0,               # Inter-compartmental clearance (L/h)
    v2=100.0,            # Peripheral volume (L)
    doses=[{"time": 0.0, "amount": 200.0}],
    t0=0.0,
    t1=24.0,
    saveat=[float(t) for t in range(25)]
)
```

---

### Three-Compartment PK Model

#### `simulate_pk_threecomp_iv_bolus(...)`

Three-compartment IV bolus simulation.

```python
result = openpkpd.simulate_pk_threecomp_iv_bolus(
    cl=5.0,              # Clearance (L/h)
    v1=10.0,             # Central volume (L)
    q2=20.0,             # Clearance to shallow peripheral (L/h)
    v2=50.0,             # Shallow peripheral volume (L)
    q3=2.0,              # Clearance to deep peripheral (L/h)
    v3=200.0,            # Deep peripheral volume (L)
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=72.0,
    saveat=[0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0, 72.0]
)
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `cl` | float | Yes | Clearance (L/h) |
| `v1` | float | Yes | Central volume (L) |
| `q2` | float | Yes | Clearance to shallow peripheral (L/h) |
| `v2` | float | Yes | Shallow peripheral volume (L) |
| `q3` | float | Yes | Clearance to deep peripheral (L/h) |
| `v3` | float | Yes | Deep peripheral volume (L) |

---

### Advanced PK Models

#### `simulate_pk_transit_absorption(...)`

Transit compartment absorption model for complex oral absorption.

```python
result = openpkpd.simulate_pk_transit_absorption(
    ktr=2.0,             # Transit rate constant (1/h)
    n_transit=3,         # Number of transit compartments
    cl=5.0,              # Clearance (L/h)
    v=50.0,              # Volume of distribution (L)
    doses=[{"time": 0.0, "amount": 200.0}],
    t0=0.0,
    t1=24.0,
    saveat=[float(t) for t in range(25)]
)
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `ktr` | float | Yes | Transit rate constant (1/h) |
| `n_transit` | int | Yes | Number of transit compartments |
| `cl` | float | Yes | Clearance (L/h) |
| `v` | float | Yes | Volume of distribution (L) |

#### `simulate_pk_michaelis_menten(...)`

Michaelis-Menten (saturable) elimination model.

```python
result = openpkpd.simulate_pk_michaelis_menten(
    vmax=10.0,           # Maximum elimination rate (mg/h)
    km=2.0,              # Michaelis constant (mg/L)
    v=50.0,              # Volume of distribution (L)
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0,
    t1=48.0,
    saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 36.0, 48.0]
)
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `vmax` | float | Yes | Maximum elimination rate (mg/h) |
| `km` | float | Yes | Michaelis constant (mg/L) |
| `v` | float | Yes | Volume of distribution (L) |

---

### PK-PD Models

#### `simulate_pkpd_direct_emax(...)`

Direct Emax PD model coupled with one-compartment PK.

```python
result = openpkpd.simulate_pkpd_direct_emax(
    # PK parameters
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    # PD parameters
    e0=0.0,              # Baseline effect
    emax=100.0,          # Maximum effect
    ec50=2.0,            # Concentration at 50% effect
    # Simulation grid
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)]
)

# Access PD output
print("Effect:", result["observations"]["effect"])
```

#### `simulate_pkpd_sigmoid_emax(...)`

Sigmoid Emax (Hill equation) PD model.

```python
result = openpkpd.simulate_pkpd_sigmoid_emax(
    # PK parameters
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    # PD parameters
    e0=10.0,             # Baseline effect
    emax=40.0,           # Maximum effect
    ec50=0.8,            # EC50 (mg/L)
    gamma=2.0,           # Hill coefficient
    # Simulation grid
    t0=0.0, t1=24.0,
    saveat=[0.0, 1.0, 4.0, 8.0, 12.0, 24.0]
)
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `e0` | float | Yes | Baseline effect |
| `emax` | float | Yes | Maximum effect |
| `ec50` | float | Yes | Concentration at 50% effect |
| `gamma` | float | Yes | Hill coefficient (steepness) |

#### `simulate_pkpd_biophase_equilibration(...)`

Effect compartment (biophase) model for delayed effects.

```python
result = openpkpd.simulate_pkpd_biophase_equilibration(
    # PK parameters
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    # PD parameters
    ke0=0.5,             # Equilibration rate constant (1/h)
    e0=0.0,              # Baseline effect
    emax=100.0,          # Maximum effect
    ec50=1.0,            # EC50 (mg/L)
    # Simulation grid
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)]
)

# Access effect compartment concentration
print("Ce:", result["states"]["A_effect"])
```

#### `simulate_pkpd_indirect_response(...)`

Indirect response turnover model.

```python
result = openpkpd.simulate_pkpd_indirect_response(
    # PK parameters
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    # PD parameters
    kin=10.0,            # Zero-order production rate
    kout=0.5,            # First-order elimination rate
    ic50=2.0,            # IC50 for inhibition
    imax=0.9,            # Maximum inhibition fraction
    baseline=None,       # Baseline (defaults to kin/kout)
    # Simulation grid
    t0=0.0, t1=72.0,
    saveat=[float(t) for t in range(73)]
)

# Access response
print("Response:", result["observations"]["response"])
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `kin` | float | Yes | Zero-order production rate |
| `kout` | float | Yes | First-order elimination rate |
| `ic50` | float | Yes | Concentration for 50% inhibition |
| `imax` | float | Yes | Maximum inhibition (0-1) |
| `baseline` | float | No | Baseline response (default: kin/kout) |

---

### Population Simulation

#### `simulate_population_iv_bolus(...)`

Population IV bolus simulation with inter-individual variability.

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

**Returns**: `dict` with keys:

- `individuals`: List of individual SimResult dicts
- `params`: List of realized parameter dicts
- `summaries`: Dict of PopulationSummary dicts
- `metadata`: Dict of population metadata

**Example - Accessing Population Results**:

```python
# Individual results
for i, individual in enumerate(result["individuals"][:5]):
    print(f"Individual {i}: Cmax = {max(individual['observations']['conc']):.2f}")

# Realized parameters
for i, params in enumerate(result["params"][:5]):
    print(f"Individual {i}: CL={params['CL']:.2f}, V={params['V']:.2f}")

# Population summary
summary = result["summaries"]["conc"]
print("Mean concentrations:", summary["mean"])
print("5th percentile:", summary["quantiles"]["0.05"])
print("95th percentile:", summary["quantiles"]["0.95"])
```

#### `simulate_population_oral(...)`

Population oral absorption simulation.

```python
result = openpkpd.simulate_population_oral(
    ka=1.5,              # Typical Ka
    cl=5.0,              # Typical CL
    v=50.0,              # Typical V
    doses=[{"time": 0.0, "amount": 200.0}],
    t0=0.0,
    t1=24.0,
    saveat=[float(t) for t in range(25)],
    n=100,
    seed=12345,
    omegas={"KA": 0.4, "CL": 0.3, "V": 0.2}
)
```

---

### Sensitivity Analysis

#### `run_sensitivity(...)`

Parameter sensitivity analysis for single or population simulations.

```python
# Single-subject sensitivity
result = openpkpd.run_sensitivity(
    model_kind="OneCompIVBolus",
    params={"CL": 5.0, "V": 50.0},
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=24.0,
    saveat=[float(t) for t in range(25)],
    sensitivity_param="CL",
    perturbation=0.1,  # 10% perturbation
    observation="conc"
)

print("AUC ratio:", result.metrics.auc_ratio)
print("Cmax ratio:", result.metrics.cmax_ratio)
```

---

### PK/PD Metrics

```python
import openpkpd

result = openpkpd.simulate_pk_iv_bolus(...)

# Exposure metrics
cmax = openpkpd.cmax(result)
t_max = openpkpd.tmax(result)
auc = openpkpd.auc_trapezoid(result)

# Half-life calculation
t_half = openpkpd.half_life(cl=5.0, v=50.0)

# PD metrics (for PKPD results)
e_min = openpkpd.emin(pkpd_result)
time_below = openpkpd.time_below(result, threshold=0.5)
auc_above = openpkpd.auc_above_baseline(result, baseline=0.0)
```

---

## NCA Module

FDA/EMA-compliant non-compartmental analysis.

### Basic NCA

```python
from openpkpd.nca import run_nca, NCAConfig

# Run NCA on concentration-time data
result = run_nca(
    times=[0, 0.5, 1, 2, 4, 8, 12, 24],
    conc=[0, 1.8, 2.0, 1.5, 1.0, 0.5, 0.25, 0.06],
    dose=100.0,
    config=NCAConfig(method="log_linear")
)

print(f"Cmax: {result.cmax:.2f}")
print(f"Tmax: {result.tmax:.2f}")
print(f"AUC0-t: {result.auc_0_t:.2f}")
print(f"AUC0-inf: {result.auc_0_inf:.2f}")
print(f"t½: {result.t_half:.2f}")
print(f"CL/F: {result.cl_f:.2f}")
```

### NCAConfig Options

```python
from openpkpd.nca import NCAConfig

config = NCAConfig(
    method="lin_log_mixed",          # "linear", "log_linear", "lin_log_mixed"
    lambda_z_min_points=3,           # Minimum points for lambda_z
    lambda_z_r2_threshold=0.9,       # R² threshold for quality
    extrapolation_max_pct=20.0,      # Max % extrapolation warning
    blq_handling="zero"              # "zero", "missing", "lloq_half"
)
```

### Individual NCA Functions

```python
from openpkpd import nca

# Exposure metrics
cmax = nca.nca_cmax(times, conc)
tmax = nca.nca_tmax(times, conc)

# AUC calculations
auc_0_t = nca.auc_0_t(times, conc, method="log_linear")
auc_0_inf, extra_pct = nca.auc_0_inf(times, conc, lambda_z)

# Terminal phase
lambda_z, t_half, r_squared = nca.estimate_lambda_z(times, conc)
```

### Population NCA

```python
from openpkpd.nca import run_population_nca, summarize_population_nca

# Run NCA for population
pop_results = run_population_nca(
    population_result,
    dose=100.0,
    config=NCAConfig()
)

# Summarize across population
summary = summarize_population_nca(pop_results)
print(f"Cmax: {summary['cmax']['mean']:.2f} (CV: {summary['cmax']['cv']:.1f}%)")
print(f"AUC: {summary['auc_0_inf']['mean']:.2f} (CV: {summary['auc_0_inf']['cv']:.1f}%)")
```

### Bioequivalence Analysis

```python
from openpkpd.nca import bioequivalence_90ci, tost_analysis, be_conclusion

# 90% confidence interval for ratio
lower, upper = bioequivalence_90ci(test_values, reference_values)
print(f"90% CI: ({lower:.2f}, {upper:.2f})")

# TOST analysis
result = tost_analysis(
    test_values,
    reference_values,
    theta_lower=0.80,
    theta_upper=1.25
)

# BE conclusion
conclusion = be_conclusion(lower, upper, theta_lower=0.80, theta_upper=1.25)
print(f"Bioequivalent: {conclusion}")
```

---

## Trial Module

Clinical trial simulation for pharmacometrics.

### Study Designs

```python
from openpkpd import trial

# Parallel design (2 arms: treatment vs placebo)
design = trial.parallel_design(n_arms=2, randomization_ratio=[1, 1])

# 2x2 crossover
design = trial.crossover_2x2(washout_duration=14.0)

# 3-period Williams design
design = trial.williams_design(washout_duration=7.0)

# 3+3 dose escalation
design = trial.dose_escalation_3plus3(
    starting_dose=10.0,
    dose_levels=[10, 25, 50, 100, 200]
)

# Bioequivalence design
design = trial.bioequivalence_design(
    n_periods=2,
    washout_duration=14.0,
    reference_formulation="tablet",
    test_formulation="capsule"
)

# Adaptive design with interim analyses
design = trial.adaptive_design(
    n_stages=3,
    interim_looks=[0.33, 0.67, 1.0],
    adaptation_rules=["sample_size", "dose_selection"]
)
```

### Dosing Regimens

```python
from openpkpd import trial

# Once daily
regimen = trial.dosing_qd(dose=100.0, duration_days=28)

# Twice daily
regimen = trial.dosing_bid(dose=50.0, duration_days=14)

# Three times daily
regimen = trial.dosing_tid(dose=25.0, duration_days=7)

# Four times daily
regimen = trial.dosing_qid(dose=25.0, duration_days=7)

# Custom schedule
regimen = trial.dosing_custom(
    dose=100.0,
    times_per_day=[0.0, 8.0, 16.0],  # Specific times
    duration_days=14
)

# Titration regimen
regimen = trial.titration_regimen(
    start_dose=25.0,
    target_dose=100.0,
    steps=[25, 50, 75, 100],
    days_per_step=7
)
```

### Virtual Population Generation

```python
from openpkpd import trial

# Default healthy volunteer population
pop = trial.generate_virtual_population(
    n=100,
    spec=trial.healthy_volunteer_spec(),
    seed=42
)

# Custom demographics
spec = trial.DemographicSpec(
    age_mean=55.0,
    age_sd=12.0,
    age_min=18.0,
    age_max=80.0,
    weight_mean=85.0,
    weight_sd=18.0,
    weight_min=50.0,
    weight_max=150.0,
    female_fraction=0.45,
    race_distribution={"white": 0.6, "black": 0.2, "asian": 0.15, "other": 0.05}
)

pop = trial.generate_virtual_population(n=200, spec=spec, seed=42)

# Summarize population
summary = trial.summarize_population(pop)
print(f"Age: {summary['age']['mean']:.1f} ({summary['age']['min']}-{summary['age']['max']})")
```

### Trial Simulation

```python
from openpkpd import trial

# Define trial specification
spec = trial.TrialSpec(
    name="Phase 2 Dose Finding",
    design=trial.parallel_design(3),
    arms=[
        trial.TreatmentArm("Placebo", regimen=trial.dosing_qd(0.0, 28)),
        trial.TreatmentArm("Low Dose", regimen=trial.dosing_qd(50.0, 28)),
        trial.TreatmentArm("High Dose", regimen=trial.dosing_qd(100.0, 28)),
    ],
    population=trial.generate_virtual_population(150, seed=42),
    dropout=trial.DropoutSpec(rate=0.05, pattern="exponential"),
    compliance=trial.ComplianceSpec(mean=0.90, sd=0.10)
)

# Run trial simulation
result = trial.simulate_trial(spec, seed=12345)

# Access results
for arm_name, arm_result in result.arms.items():
    print(f"{arm_name}: n={arm_result.n_completed}")
```

### Power Analysis

```python
from openpkpd import trial

# Analytical power calculation
power = trial.estimate_power_analytical(
    n_per_arm=50,
    effect_size=0.5,      # Cohen's d
    sd=1.0,
    alpha=0.05
)
print(f"Power: {power.power:.1%}")

# Sample size estimation
result = trial.estimate_sample_size(
    target_power=0.80,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05
)
print(f"Required n per arm: {result.n_per_arm}")

# Alpha spending for interim analyses
alpha = trial.alpha_spending_function(
    information_fraction=0.5,
    total_alpha=0.05,
    method="obrien_fleming"
)
```

### Arm Comparison

```python
from openpkpd import trial

# Compare treatment arms
comparison = trial.compare_arms(
    treatment_values=[1.2, 1.5, 1.1, 1.8, 1.4],
    control_values=[0.9, 1.0, 0.8, 1.1, 0.95],
    test="ttest"
)

print(f"Treatment effect: {comparison.difference:.2f}")
print(f"95% CI: ({comparison.ci_lower:.2f}, {comparison.ci_upper:.2f})")
print(f"p-value: {comparison.p_value:.4f}")

# Responder analysis
result = trial.responder_analysis(
    values=[1.2, 0.8, 1.5, 0.6, 1.1, 2.0],
    threshold=1.0
)
print(f"Responder rate: {result.rate:.1%}")
```

---

## Visualization Module

Professional visualization with matplotlib and plotly backends.

### Backend Selection

```python
from openpkpd import viz

# Set backend
viz.set_backend("matplotlib")  # or "plotly"

# Check available backends
print(viz.available_backends())

# Get current backend
print(viz.get_backend())
```

### PK Plots

```python
from openpkpd import viz
import openpkpd

result = openpkpd.simulate_pk_iv_bolus(...)

# Concentration-time plot
fig = viz.plot_conc_time(result, log_scale=False)
fig.show()

# Multiple simulations
results = [sim1, sim2, sim3]
fig = viz.plot_multi_conc_time(results, labels=["Low", "Medium", "High"])

# Population spaghetti plot
pop_result = openpkpd.simulate_population_iv_bolus(...)
fig = viz.plot_spaghetti(pop_result, n_subjects=50, alpha=0.3)

# Mean with confidence ribbon
fig = viz.plot_mean_ribbon(pop_result, ci_levels=[0.05, 0.95], show_median=True)

# Individual fits (grid of subjects)
fig = viz.plot_individual_fits(pop_result, n_cols=4)
```

### NCA Plots

```python
from openpkpd import viz
from openpkpd.nca import run_nca

nca_result = run_nca(times, conc, dose=100.0)

# Lambda-z fit visualization
fig = viz.plot_lambda_z_fit(nca_result, times, conc, show_excluded=True)

# AUC visualization
fig = viz.plot_auc_visualization(times, conc, nca_result, show_extrapolation=True)

# Dose proportionality
results = [nca_25mg, nca_50mg, nca_100mg]
doses = [25, 50, 100]
fig = viz.plot_dose_proportionality(results, doses, metric="auc")
```

### PKPD Plots

```python
from openpkpd import viz

pkpd_result = openpkpd.simulate_pkpd_direct_emax(...)

# Effect-concentration plot
fig = viz.plot_effect_conc(pkpd_result, color_by_time=True)

# Hysteresis plot (counterclockwise loop)
fig = viz.plot_hysteresis(pkpd_result, arrow_frequency=5)

# Dose-response curve
dose_results = [result_10mg, result_25mg, result_50mg, result_100mg]
doses = [10, 25, 50, 100]
fig = viz.plot_dose_response(dose_results, doses, fit_model="emax")
```

### Population Plots

```python
from openpkpd import viz

# Visual Predictive Check (VPC)
fig = viz.plot_vpc(
    pop_result,
    observed_data=observed_df,
    prediction_intervals=[0.05, 0.50, 0.95]
)

# Parameter distributions
fig = viz.plot_parameter_distributions(
    pop_result,
    parameters=["CL", "V"],
    plot_type="histogram"  # or "kde", "boxplot"
)

# Forest plot for subgroup analysis
effects = [
    {"label": "Overall", "estimate": 0.5, "ci_lower": 0.2, "ci_upper": 0.8},
    {"label": "Age < 65", "estimate": 0.6, "ci_lower": 0.25, "ci_upper": 0.95},
    {"label": "Age >= 65", "estimate": 0.4, "ci_lower": 0.1, "ci_upper": 0.7},
]
fig = viz.plot_forest(effects, reference_line=0.0)

# Boxplot comparison
fig = viz.plot_boxplot(pop_result, groups=["Treatment", "Control"], metric="cmax")
```

### Trial Plots

```python
from openpkpd import viz
from openpkpd import trial

# Power curve
power_results = [
    {"n": 20, "power": 0.45},
    {"n": 40, "power": 0.70},
    {"n": 60, "power": 0.85},
    {"n": 80, "power": 0.92},
]
fig = viz.plot_power_curve(power_results, target_power=0.80)

# Tornado plot (sensitivity)
sensitivity = [
    {"parameter": "CL", "low": -0.3, "high": 0.25},
    {"parameter": "V", "low": -0.15, "high": 0.18},
    {"parameter": "Ka", "low": -0.4, "high": 0.35},
]
fig = viz.plot_tornado(sensitivity, baseline_value=0.0)

# Kaplan-Meier survival
fig = viz.plot_kaplan_meier(
    time_to_event=[10, 15, 20, 25, 30],
    event_occurred=[1, 1, 0, 1, 0],
    groups=["Treatment", "Treatment", "Control", "Control", "Treatment"]
)

# Endpoint distribution by arm
fig = viz.plot_endpoint_distribution(trial_result, endpoint="response", by_arm=True)
```

### Themes

```python
from openpkpd import viz

# Set theme
viz.set_theme("publication")  # or "presentation", "default"

# Available themes
print(viz.available_themes())

# Access color palette
colors = viz.OPENPKPD_COLORS
print(colors["primary"])
```

---

## Artifact Operations

### Writing Artifacts

```python
openpkpd.write_single_artifact(
    "simulation.json",
    model={
        "kind": "TwoCompIVBolus",
        "params": {"CL": 5.0, "V1": 50.0, "Q": 2.0, "V2": 100.0},
        "doses": [{"time": 0.0, "amount": 100.0}]
    },
    grid={"t0": 0.0, "t1": 24.0, "saveat": list(range(25))},
    solver={"alg": "Tsit5", "reltol": 1e-10, "abstol": 1e-12}
)
```

**Supported Model Kinds**:

- `"OneCompIVBolus"`, `"OneCompOralFirstOrder"`
- `"TwoCompIVBolus"`, `"TwoCompOral"`
- `"ThreeCompIVBolus"`
- `"TransitAbsorption"`, `"MichaelisMentenElimination"`

### Replaying Artifacts

```python
# Replay any artifact type
result = openpkpd.replay_artifact("validation/golden/pk_iv_bolus.json")
print("Concentrations:", result["observations"]["conc"])

# Population artifact
result = openpkpd.replay_artifact("validation/golden/population_iv_bolus.json")
print("Number of individuals:", len(result["individuals"]))

# Sensitivity artifact
result = openpkpd.replay_artifact("validation/golden/sensitivity_single.json")
print("Metrics:", result["metrics"])
```

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
4. **Sparse Output**: Use fewer `saveat` points for large populations

```python
# Good: Initialize once
import openpkpd
openpkpd.init_julia()

for params in parameter_sets:
    result = openpkpd.simulate_pk_iv_bolus(...)
```

---

## Complete API Reference

### Core Functions

| Function | Description |
|----------|-------------|
| `init_julia()` | Initialize Julia runtime |
| `version()` | Get OpenPKPD version |

### Simulation Functions

| Function | Model Type |
|----------|------------|
| `simulate_pk_iv_bolus` | One-compartment IV bolus |
| `simulate_pk_oral_first_order` | One-compartment oral |
| `simulate_pk_twocomp_iv_bolus` | Two-compartment IV bolus |
| `simulate_pk_twocomp_oral` | Two-compartment oral |
| `simulate_pk_threecomp_iv_bolus` | Three-compartment IV bolus |
| `simulate_pk_transit_absorption` | Transit compartment |
| `simulate_pk_michaelis_menten` | Michaelis-Menten elimination |
| `simulate_pkpd_direct_emax` | Direct Emax PD |
| `simulate_pkpd_sigmoid_emax` | Sigmoid Emax (Hill) PD |
| `simulate_pkpd_biophase_equilibration` | Effect compartment PD |
| `simulate_pkpd_indirect_response` | Indirect response PD |
| `simulate_population_iv_bolus` | Population IV bolus |
| `simulate_population_oral` | Population oral |
| `run_sensitivity` | Sensitivity analysis |

### NCA Functions

| Function | Description |
|----------|-------------|
| `run_nca` | Full NCA analysis |
| `run_population_nca` | Population NCA |
| `summarize_population_nca` | NCA summary statistics |
| `nca_cmax`, `nca_tmax` | Peak metrics |
| `auc_0_t`, `auc_0_inf` | AUC calculations |
| `estimate_lambda_z` | Terminal slope |
| `bioequivalence_90ci` | BE confidence interval |
| `tost_analysis` | Two one-sided t-tests |

### Trial Functions

| Function | Description |
|----------|-------------|
| `parallel_design`, `crossover_2x2` | Study designs |
| `dosing_qd`, `dosing_bid` | Dosing regimens |
| `generate_virtual_population` | Virtual subjects |
| `simulate_trial` | Trial simulation |
| `estimate_power_analytical` | Power calculation |
| `estimate_sample_size` | Sample size |

### Visualization Functions

| Function | Category |
|----------|----------|
| `plot_conc_time`, `plot_spaghetti` | PK plots |
| `plot_lambda_z_fit`, `plot_auc_visualization` | NCA plots |
| `plot_effect_conc`, `plot_hysteresis` | PKPD plots |
| `plot_vpc`, `plot_forest` | Population plots |
| `plot_power_curve`, `plot_tornado` | Trial plots |
