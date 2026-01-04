# Parameter Estimation (NLME)

OpenPKPD provides nonlinear mixed-effects (NLME) parameter estimation for population pharmacokinetic/pharmacodynamic modeling. Three estimation methods are available: FOCE-I, SAEM, and Laplacian approximation.

## Overview

Parameter estimation fits model parameters to observed data while accounting for:
- **Fixed effects (theta)**: Population typical values (e.g., CL, V)
- **Random effects (eta)**: Individual deviations from population means
- **Variance components (omega)**: Between-subject variability
- **Residual error (sigma)**: Unexplained variability

## Estimation Methods

### FOCE-I (First-Order Conditional Estimation with Interaction)

The gold standard for NLME estimation, used by NONMEM and other tools.

**Algorithm:**
1. For each individual, find the mode of the conditional distribution (eta)
2. Linearize the model around the conditional mode
3. Optimize population parameters using the linearized likelihood
4. Iterate until convergence

**Advantages:**
- Well-established with decades of validation
- Handles interaction between eta and epsilon
- Good balance of speed and accuracy

**Julia:**
```julia
config = EstimationConfig(
    FOCEIMethod(
        max_inner_iter=100,  # Max iterations for eta optimization
        inner_tol=1e-6,       # Tolerance for eta convergence
        centered=false        # Use uncentered parameterization
    ),
    theta_init=[5.0, 50.0],
    theta_lower=[0.1, 1.0],
    theta_upper=[100.0, 500.0],
    omega_init=diagm([0.09, 0.04]),
    omega_structure=:diagonal,
    sigma_init=ResidualErrorSpec(ProportionalError(), (sigma=0.1,), :conc, UInt64(0)),
    max_iter=500,
    compute_se=true
)

result = estimate(observed_data, model_spec, config; grid=grid, solver=solver)
```

**Python:**
```python
from openpkpd.estimation import estimate, EstimationConfig, FOCEIMethod

config = EstimationConfig(
    method=FOCEIMethod(max_inner_iter=100, inner_tol=1e-6),
    theta_init=[5.0, 50.0],
    theta_lower=[0.1, 1.0],
    theta_upper=[100.0, 500.0],
    omega_init=[[0.09, 0], [0, 0.04]],
    sigma_init={"kind": "proportional", "value": 0.1},
    max_iter=500,
    compute_se=True
)

result = estimate(observed_data, model_spec, config)
```

### SAEM (Stochastic Approximation EM)

A Monte Carlo-based method that can handle complex models.

**Algorithm:**
1. **E-step**: Sample from the conditional distribution of eta using MCMC
2. **M-step**: Update population parameters using sufficient statistics
3. Use stochastic approximation to smooth updates
4. Continue until convergence

**Advantages:**
- Robust for complex models
- Handles multiple modes in the likelihood
- Used by Monolix

**Julia:**
```julia
config = EstimationConfig(
    SAEMMethod(
        n_burn=200,    # Burn-in iterations
        n_iter=300,    # Main iterations
        n_chains=3     # MCMC chains per subject
    ),
    theta_init=[5.0, 50.0],
    # ... same as FOCE-I
)
```

**Python:**
```python
from openpkpd.estimation import SAEMMethod

config = EstimationConfig(
    method=SAEMMethod(n_burn=200, n_iter=300, n_chains=3),
    # ... same as FOCE-I
)
```

### Laplacian Approximation

A simpler approximation that works well for sparse data.

**Algorithm:**
1. Approximate the marginal likelihood using Laplace's method
2. Find the mode of the joint density
3. Use the Hessian at the mode for the approximation

**Advantages:**
- Fast computation
- Good baseline for comparison
- Works well with sparse data

**Julia:**
```julia
config = EstimationConfig(
    LaplacianMethod(),
    theta_init=[5.0, 50.0],
    # ... same as FOCE-I
)
```

## Residual Error Models

OpenPKPD supports four residual error models:

### Additive Error

$$Y = F + \epsilon, \quad \epsilon \sim N(0, \sigma^2)$$

```julia
error_spec = ResidualErrorSpec(AdditiveError(), (sigma=0.5,), :conc, seed)
```

### Proportional Error

$$Y = F \cdot (1 + \epsilon), \quad \epsilon \sim N(0, \sigma^2)$$

```julia
error_spec = ResidualErrorSpec(ProportionalError(), (sigma=0.1,), :conc, seed)
```

### Combined Error

$$Y = F + F \cdot \epsilon_1 + \epsilon_2$$

```julia
error_spec = ResidualErrorSpec(CombinedError(), (sigma_add=0.5, sigma_prop=0.1), :conc, seed)
```

### Exponential Error

$$Y = F \cdot e^{\epsilon}, \quad \epsilon \sim N(0, \sigma^2)$$

```julia
error_spec = ResidualErrorSpec(ExponentialError(), (sigma=0.1,), :conc, seed)
```

## Estimation Results

The `EstimationResult` contains:

| Field | Description |
|-------|-------------|
| `theta` | Final population parameter estimates |
| `theta_se` | Standard errors for theta |
| `omega` | Final variance-covariance matrix |
| `omega_se` | Standard errors for omega elements |
| `sigma` | Final residual error parameters |
| `sigma_se` | Standard errors for sigma |
| `ofv` | Objective function value (-2LL) |
| `aic` | Akaike Information Criterion |
| `bic` | Bayesian Information Criterion |
| `individuals` | Per-subject estimates and diagnostics |
| `convergence` | Whether estimation converged |
| `n_iterations` | Number of iterations |

### Individual Estimates

Each individual result contains:

| Field | Description |
|-------|-------------|
| `subject_id` | Subject identifier |
| `eta` | Individual random effects |
| `ipred` | Individual predictions |
| `cwres` | Conditional weighted residuals |
| `ofv_contribution` | Individual contribution to OFV |

## Diagnostics

### Goodness-of-Fit Plots

**Python:**
```python
from openpkpd import viz

# Observed vs Population Predicted
viz.plot_goodness_of_fit(result, plot_type="obs_vs_pred")

# Observed vs Individual Predicted
viz.plot_goodness_of_fit(result, plot_type="obs_vs_ipred")

# CWRES vs Time
viz.plot_goodness_of_fit(result, plot_type="cwres_vs_time")

# CWRES vs PRED
viz.plot_goodness_of_fit(result, plot_type="cwres_vs_pred")

# CWRES QQ Plot
viz.plot_goodness_of_fit(result, plot_type="cwres_qq")

# ETA Distribution
viz.plot_goodness_of_fit(result, plot_type="eta_dist")

# 4-Panel Summary
viz.plot_estimation_summary(result)
```

### Model Comparison

```python
# Compare models using AIC/BIC
print(f"Model 1 AIC: {result1.aic:.2f}, BIC: {result1.bic:.2f}")
print(f"Model 2 AIC: {result2.aic:.2f}, BIC: {result2.bic:.2f}")

# Likelihood ratio test (nested models)
delta_ofv = result1.ofv - result2.ofv
df = n_params2 - n_params1
from scipy.stats import chi2
p_value = 1 - chi2.cdf(delta_ofv, df)
```

## Configuration Options

### Theta (Fixed Effects)

| Parameter | Description |
|-----------|-------------|
| `theta_init` | Initial parameter values |
| `theta_lower` | Lower bounds |
| `theta_upper` | Upper bounds |

### Omega (Between-Subject Variability)

| Parameter | Description |
|-----------|-------------|
| `omega_init` | Initial variance-covariance matrix |
| `omega_structure` | `:diagonal` or `:block` |

### Optimization

| Parameter | Description |
|-----------|-------------|
| `max_iter` | Maximum outer iterations |
| `compute_se` | Whether to compute standard errors |

## CLI Usage

```bash
./packages/cli/bin/openpkpd estimate --spec estimate_spec.json --out fit.json
```

**Specification Format:**
```json
{
  "data": "observed.csv",
  "model": {
    "kind": "OneCompOralFirstOrder",
    "params": {"Ka": 1.0, "CL": 5.0, "V": 50.0}
  },
  "estimation": {
    "method": "FOCE-I",
    "theta_init": [1.0, 5.0, 50.0],
    "theta_lower": [0.1, 0.1, 1.0],
    "theta_upper": [10.0, 100.0, 500.0],
    "omega_init": [[0.09, 0, 0], [0, 0.09, 0], [0, 0, 0.04]],
    "omega_structure": "diagonal",
    "sigma": {"kind": "proportional", "value": 0.1},
    "max_iter": 500,
    "compute_se": true
  },
  "grid": {"t0": 0.0, "t1": 24.0},
  "solver": {"alg": "Tsit5", "reltol": 1e-8, "abstol": 1e-10}
}
```

## Tips and Best Practices

1. **Initial Values**: Start with reasonable initial values based on prior knowledge or NCA results

2. **Bounds**: Set physiologically plausible bounds to prevent estimation from exploring unrealistic parameter space

3. **Convergence**: Check that OFV is stable and standard errors are reasonable

4. **Residual Error**: Choose the error model that matches your data (proportional for concentrations that span orders of magnitude)

5. **Omega Structure**: Start with diagonal omega; add correlations only if supported by data

6. **Validation**: Use VPC to validate the final model

## See Also

- [Models](models.md) - Available PK/PD models
- [Population](population.md) - Population simulation
- [VPC](vpc.md) - Visual Predictive Checks
- [Data Import](data.md) - Loading observed data
