"""
OpenPKPD Parameter Estimation Module

This module provides Python bindings for NLME (Nonlinear Mixed Effects)
parameter estimation methods including FOCE-I, SAEM, and Laplacian.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Union
from pathlib import Path

from .._core import _require_julia


@dataclass
class IndividualEstimate:
    """Individual-level estimation results."""
    subject_id: str
    eta: List[float]
    eta_se: Optional[List[float]]
    ipred: List[float]
    pred: List[float]
    cwres: List[float]
    iwres: List[float]
    wres: List[float]
    ofv_contribution: float


@dataclass
class EstimationResult:
    """Result from NLME parameter estimation."""
    method: str
    theta: List[float]
    theta_se: Optional[List[float]]
    theta_rse: Optional[List[float]]
    theta_ci_lower: Optional[List[float]]
    theta_ci_upper: Optional[List[float]]
    omega: List[List[float]]
    omega_se: Optional[List[List[float]]]
    omega_corr: List[List[float]]
    sigma: Dict[str, Any]
    sigma_se: Optional[Dict[str, Any]]
    ofv: float
    aic: float
    bic: float
    convergence: bool
    n_iterations: int
    gradient_norm: float
    condition_number: float
    eigenvalue_ratio: float
    covariance_successful: bool
    individuals: List[IndividualEstimate]
    runtime_seconds: float
    warnings: List[str]


@dataclass
class EstimationConfig:
    """Configuration for parameter estimation."""
    method: str = "foce"  # foce, saem, or laplacian
    theta_init: Optional[List[float]] = None
    theta_lower: Optional[List[float]] = None
    theta_upper: Optional[List[float]] = None
    theta_names: Optional[List[str]] = None
    omega_init: Optional[List[List[float]]] = None
    omega_names: Optional[List[str]] = None
    omega_structure: str = "diagonal"  # diagonal or block
    sigma_type: str = "proportional"  # additive, proportional, combined
    sigma_init: float = 0.1
    max_iter: int = 500
    tol: float = 1e-6
    compute_se: bool = True
    compute_ci: bool = True
    ci_level: float = 0.95
    verbose: bool = False
    seed: int = 12345


def estimate(
    observed_data: Dict[str, Any],
    model_kind: str,
    config: EstimationConfig,
    grid: Dict[str, Any],
    solver: Optional[Dict[str, Any]] = None,
) -> EstimationResult:
    """
    Run NLME parameter estimation.

    Supports three estimation methods:
    - FOCE-I: First-Order Conditional Estimation with Interaction (gold standard)
    - SAEM: Stochastic Approximation Expectation Maximization
    - Laplacian: Simple Laplacian approximation

    Args:
        observed_data: Dict with keys:
            - subjects: List of subject data dicts
            - Each subject has: subject_id, times, observations, doses
        model_kind: PK model type (e.g., "OneCompIVBolus", "TwoCompOral")
        config: EstimationConfig with initial values and method settings
        grid: Simulation grid dict with t0, t1, saveat
        solver: Optional solver settings

    Returns:
        EstimationResult containing estimated parameters, SEs, and diagnostics

    Example:
        >>> config = EstimationConfig(
        ...     method="foce",
        ...     theta_init=[10.0, 50.0],
        ...     theta_names=["CL", "V"],
        ...     omega_init=[[0.09, 0.0], [0.0, 0.04]],
        ...     sigma_type="proportional",
        ...     sigma_init=0.1
        ... )
        >>> result = estimate(observed_data, "OneCompIVBolus", config, grid)
        >>> print(f"CL = {result.theta[0]} +/- {result.theta_se[0]}")
    """
    jl = _require_julia()

    # Build observed data structure
    SubjectData = jl.OpenPKPDCore.SubjectData
    ObservedData = jl.OpenPKPDCore.ObservedData
    DoseEvent = jl.OpenPKPDCore.DoseEvent

    subjects = []
    for subj in observed_data["subjects"]:
        doses = [DoseEvent(float(d["time"]), float(d["amount"])) for d in subj.get("doses", [])]
        subjects.append(SubjectData(
            subj["subject_id"],
            [float(t) for t in subj["times"]],
            [float(o) for o in subj["observations"]],
            doses
        ))

    obs = ObservedData(subjects)

    # Build model spec (just for type inference - estimation creates individual specs)
    model_spec = _build_base_model_spec(jl, model_kind, config.theta_init)

    # Build estimation method
    if config.method == "foce":
        method = jl.OpenPKPDCore.FOCEIMethod()
    elif config.method == "saem":
        method = jl.OpenPKPDCore.SAEMMethod(n_burn=100, n_iter=200)
    elif config.method == "laplacian":
        method = jl.OpenPKPDCore.LaplacianMethod()
    else:
        raise ValueError(f"Unknown estimation method: {config.method}")

    # Build sigma specification
    sigma_spec = _build_sigma_spec(jl, config.sigma_type, config.sigma_init)

    # Build omega matrix
    import numpy as np
    omega_init = np.array(config.omega_init) if config.omega_init else np.eye(len(config.theta_init)) * 0.09

    # Build theta bounds
    n_theta = len(config.theta_init)
    theta_lower = config.theta_lower if config.theta_lower else [1e-6] * n_theta
    theta_upper = config.theta_upper if config.theta_upper else [1e6] * n_theta

    # Build omega structure
    omega_structure = jl.OpenPKPDCore.DiagonalOmega() if config.omega_structure == "diagonal" else jl.OpenPKPDCore.BlockOmega()

    # Build estimation config
    est_config = jl.OpenPKPDCore.EstimationConfig(
        method,
        theta_init=[float(x) for x in config.theta_init],
        theta_lower=[float(x) for x in theta_lower],
        theta_upper=[float(x) for x in theta_upper],
        theta_names=[jl.Symbol(n) for n in (config.theta_names or [f"theta_{i}" for i in range(n_theta)])],
        omega_init=omega_init.tolist(),
        omega_names=[jl.Symbol(n) for n in (config.omega_names or [f"eta_{i}" for i in range(omega_init.shape[0])])],
        omega_structure=omega_structure,
        sigma_init=sigma_spec,
        max_iter=config.max_iter,
        tol=config.tol,
        compute_se=config.compute_se,
        compute_ci=config.compute_ci,
        ci_level=config.ci_level,
        verbose=config.verbose,
        seed=jl.UInt64(config.seed),
    )

    # Build grid
    grid_jl = jl.OpenPKPDCore.SimGrid(
        float(grid["t0"]),
        float(grid["t1"]),
        [float(x) for x in grid["saveat"]],
    )

    # Build solver
    if solver is None:
        solver_jl = jl.OpenPKPDCore.SolverSpec(jl.Symbol("Tsit5"), 1e-10, 1e-12, 10**7)
    else:
        solver_jl = jl.OpenPKPDCore.SolverSpec(
            jl.Symbol(solver.get("alg", "Tsit5")),
            float(solver.get("reltol", 1e-10)),
            float(solver.get("abstol", 1e-12)),
            int(solver.get("maxiters", 10**7)),
        )

    # Run estimation
    result = jl.OpenPKPDCore.estimate(obs, model_spec, est_config, grid=grid_jl, solver=solver_jl)

    # Convert result to Python
    return _convert_estimation_result(result, config.method)


def _build_base_model_spec(jl, model_kind: str, theta_init: List[float]):
    """Build a model spec for the specified model kind."""
    ModelSpec = jl.OpenPKPDCore.ModelSpec
    DoseEvent = jl.OpenPKPDCore.DoseEvent

    doses = [DoseEvent(0.0, 100.0)]  # Placeholder

    if model_kind == "OneCompIVBolus":
        Kind = jl.OpenPKPDCore.OneCompIVBolus
        Params = jl.OpenPKPDCore.OneCompIVBolusParams
        params = Params(float(theta_init[0]), float(theta_init[1]))
    elif model_kind == "OneCompOralFirstOrder":
        Kind = jl.OpenPKPDCore.OneCompOralFirstOrder
        Params = jl.OpenPKPDCore.OneCompOralFirstOrderParams
        params = Params(float(theta_init[0]), float(theta_init[1]), float(theta_init[2]))
    elif model_kind == "TwoCompIVBolus":
        Kind = jl.OpenPKPDCore.TwoCompIVBolus
        Params = jl.OpenPKPDCore.TwoCompIVBolusParams
        params = Params(float(theta_init[0]), float(theta_init[1]), float(theta_init[2]), float(theta_init[3]))
    elif model_kind == "TwoCompOral":
        Kind = jl.OpenPKPDCore.TwoCompOral
        Params = jl.OpenPKPDCore.TwoCompOralParams
        params = Params(*[float(x) for x in theta_init[:5]])
    else:
        raise ValueError(f"Unsupported model kind: {model_kind}")

    return ModelSpec(Kind(), "estimation", params, doses)


def _build_sigma_spec(jl, sigma_type: str, sigma_init: float):
    """Build residual error specification."""
    if sigma_type == "additive":
        kind = jl.OpenPKPDCore.AdditiveError()
        params = jl.OpenPKPDCore.AdditiveErrorParams(float(sigma_init))
    elif sigma_type == "proportional":
        kind = jl.OpenPKPDCore.ProportionalError()
        params = jl.OpenPKPDCore.ProportionalErrorParams(float(sigma_init))
    elif sigma_type == "combined":
        kind = jl.OpenPKPDCore.CombinedError()
        params = jl.OpenPKPDCore.CombinedErrorParams(float(sigma_init), float(sigma_init))
    elif sigma_type == "exponential":
        kind = jl.OpenPKPDCore.ExponentialError()
        params = jl.OpenPKPDCore.ExponentialErrorParams(float(sigma_init))
    else:
        raise ValueError(f"Unknown sigma type: {sigma_type}")

    return jl.OpenPKPDCore.ResidualErrorSpec(kind, params, jl.Symbol("conc"), jl.UInt64(1))


def _convert_estimation_result(result, method: str) -> EstimationResult:
    """Convert Julia EstimationResult to Python dataclass."""
    # Convert individuals
    individuals = []
    for ind in result.individuals:
        individuals.append(IndividualEstimate(
            subject_id=str(ind.subject_id),
            eta=list(ind.eta),
            eta_se=list(ind.eta_se) if ind.eta_se is not None else None,
            ipred=list(ind.ipred),
            pred=list(ind.pred),
            cwres=list(ind.cwres),
            iwres=list(ind.iwres),
            wres=list(ind.wres),
            ofv_contribution=float(ind.ofv_contribution),
        ))

    # Convert omega to list of lists
    omega = [[float(result.omega[i, j]) for j in range(result.omega.shape[1])]
             for i in range(result.omega.shape[0])]
    omega_corr = [[float(result.omega_corr[i, j]) for j in range(result.omega_corr.shape[1])]
                  for i in range(result.omega_corr.shape[0])]

    omega_se = None
    if result.omega_se is not None:
        omega_se = [[float(result.omega_se[i, j]) for j in range(result.omega_se.shape[1])]
                    for i in range(result.omega_se.shape[0])]

    return EstimationResult(
        method=method,
        theta=list(result.theta),
        theta_se=list(result.theta_se) if result.theta_se is not None else None,
        theta_rse=list(result.theta_rse) if result.theta_rse is not None else None,
        theta_ci_lower=list(result.theta_ci_lower) if result.theta_ci_lower is not None else None,
        theta_ci_upper=list(result.theta_ci_upper) if result.theta_ci_upper is not None else None,
        omega=omega,
        omega_se=omega_se,
        omega_corr=omega_corr,
        sigma={"type": "estimated"},  # Simplified
        sigma_se=None,
        ofv=float(result.ofv),
        aic=float(result.aic),
        bic=float(result.bic),
        convergence=bool(result.convergence),
        n_iterations=int(result.n_iterations),
        gradient_norm=float(result.gradient_norm),
        condition_number=float(result.condition_number) if hasattr(result, 'condition_number') else float('nan'),
        eigenvalue_ratio=float(result.eigenvalue_ratio) if hasattr(result, 'eigenvalue_ratio') else float('nan'),
        covariance_successful=bool(result.covariance_successful),
        individuals=individuals,
        runtime_seconds=float(result.runtime_seconds),
        warnings=list(result.warnings) if hasattr(result, 'warnings') else [],
    )


__all__ = [
    "estimate",
    "EstimationConfig",
    "EstimationResult",
    "IndividualEstimate",
]
