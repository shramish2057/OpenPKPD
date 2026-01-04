"""
OpenPKPD Visual Predictive Check (VPC) Module

This module provides Python bindings for VPC analysis including
standard VPC and prediction-corrected VPC (pcVPC).
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from .._core import _require_julia


@dataclass
class VPCPercentile:
    """Percentile data for a VPC bin."""
    level: float
    observed: float
    simulated_median: float
    simulated_lower: float
    simulated_upper: float


@dataclass
class VPCBin:
    """A single bin in VPC analysis."""
    bin_id: int
    time_min: float
    time_max: float
    time_midpoint: float
    n_observed: int
    percentiles: List[VPCPercentile]


@dataclass
class VPCResult:
    """Result from VPC analysis."""
    n_subjects_observed: int
    n_observations_observed: int
    n_simulations: int
    prediction_corrected: bool
    pi_levels: List[float]
    ci_level: float
    bins: List[VPCBin]


@dataclass
class VPCConfig:
    """Configuration for VPC analysis."""
    pi_levels: List[float] = None
    ci_level: float = 0.95
    n_bins: int = 8
    binning: str = "quantile"  # quantile, equal_width
    n_simulations: int = 200
    n_bootstrap: int = 500
    prediction_corrected: bool = False
    stratify_by: Optional[List[str]] = None
    lloq: Optional[float] = None
    seed: int = 12345

    def __post_init__(self):
        if self.pi_levels is None:
            self.pi_levels = [0.05, 0.50, 0.95]


def compute_vpc(
    observed_data: Dict[str, Any],
    population_spec: Dict[str, Any],
    grid: Dict[str, Any],
    config: Optional[VPCConfig] = None,
    solver: Optional[Dict[str, Any]] = None,
    error_spec: Optional[Dict[str, Any]] = None,
) -> VPCResult:
    """
    Compute Visual Predictive Check.

    Args:
        observed_data: Dict with subject data containing:
            - subjects: List of dicts with subject_id, times, observations, doses
        population_spec: Dict with population specification:
            - model: Model specification with kind, params
            - iiv: IIV specification with omegas, n, seed
        grid: Simulation grid with t0, t1, saveat
        config: VPCConfig with analysis settings (optional)
        solver: Solver settings (optional)
        error_spec: Residual error specification (optional)

    Returns:
        VPCResult with binned percentile comparisons

    Example:
        >>> config = VPCConfig(
        ...     n_simulations=200,
        ...     pi_levels=[0.05, 0.50, 0.95],
        ...     n_bins=8
        ... )
        >>> result = compute_vpc(observed, pop_spec, grid, config)
        >>> for bin in result.bins:
        ...     print(f"Bin {bin.bin_id}: {bin.n_observed} observations")
    """
    jl = _require_julia()

    if config is None:
        config = VPCConfig()

    # Build observed data
    obs = _build_observed_data(jl, observed_data)

    # Build population spec
    pop_spec = _build_population_spec(jl, population_spec)

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

    # Build VPC config
    binning = jl.OpenPKPDCore.QuantileBinning(config.n_bins) if config.binning == "quantile" else jl.OpenPKPDCore.EqualWidthBinning(config.n_bins)

    vpc_config = jl.OpenPKPDCore.VPCConfig(
        pi_levels=[float(x) for x in config.pi_levels],
        ci_level=float(config.ci_level),
        binning=binning,
        n_simulations=config.n_simulations,
        n_bootstrap=config.n_bootstrap,
        prediction_corrected=config.prediction_corrected,
        stratify_by=[jl.Symbol(s) for s in config.stratify_by] if config.stratify_by else [],
        lloq=float(config.lloq) if config.lloq else None,
        seed=jl.UInt64(config.seed),
    )

    # Build error spec if provided
    if error_spec is not None:
        error_jl = _build_error_spec(jl, error_spec)
    else:
        error_jl = None

    # Run VPC
    if config.prediction_corrected:
        result = jl.OpenPKPDCore.compute_pcvpc(
            obs, pop_spec, grid_jl, solver_jl,
            config=vpc_config, error_spec=error_jl
        )
    else:
        result = jl.OpenPKPDCore.compute_vpc(
            obs, pop_spec, grid_jl, solver_jl,
            config=vpc_config, error_spec=error_jl
        )

    return _convert_vpc_result(result, config)


def compute_pcvpc(
    observed_data: Dict[str, Any],
    population_spec: Dict[str, Any],
    grid: Dict[str, Any],
    config: Optional[VPCConfig] = None,
    solver: Optional[Dict[str, Any]] = None,
    error_spec: Optional[Dict[str, Any]] = None,
) -> VPCResult:
    """
    Compute Prediction-Corrected Visual Predictive Check (pcVPC).

    pcVPC normalizes both observed and simulated data to remove the effect
    of model-predicted trends, making it useful when there's high variability
    in dosing or covariates.

    Args:
        Same as compute_vpc

    Returns:
        VPCResult with prediction-corrected percentiles
    """
    if config is None:
        config = VPCConfig(prediction_corrected=True)
    else:
        config.prediction_corrected = True

    return compute_vpc(observed_data, population_spec, grid, config, solver, error_spec)


def _build_observed_data(jl, data: Dict[str, Any]):
    """Build ObservedData from Python dict."""
    SubjectData = jl.OpenPKPDCore.SubjectData
    ObservedData = jl.OpenPKPDCore.ObservedData
    DoseEvent = jl.OpenPKPDCore.DoseEvent

    subjects = []
    for subj in data["subjects"]:
        doses = [DoseEvent(float(d["time"]), float(d["amount"])) for d in subj.get("doses", [])]
        subjects.append(SubjectData(
            subj["subject_id"],
            [float(t) for t in subj["times"]],
            [float(o) for o in subj["observations"]],
            doses
        ))

    return ObservedData(subjects)


def _build_population_spec(jl, spec: Dict[str, Any]):
    """Build PopulationSpec from Python dict."""
    from ..bridge import _to_julia_vector

    ModelSpec = jl.OpenPKPDCore.ModelSpec
    DoseEvent = jl.OpenPKPDCore.DoseEvent
    IIVSpec = jl.OpenPKPDCore.IIVSpec
    LogNormalIIV = jl.OpenPKPDCore.LogNormalIIV
    PopulationSpec = jl.OpenPKPDCore.PopulationSpec

    model = spec["model"]
    doses_raw = model.get("doses", [{"time": 0.0, "amount": 100.0}])
    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses_raw]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    # Build model based on kind
    kind_str = model["kind"]
    params = model["params"]

    if kind_str == "OneCompIVBolus":
        Kind = jl.OpenPKPDCore.OneCompIVBolus
        Params = jl.OpenPKPDCore.OneCompIVBolusParams
        model_params = Params(float(params["CL"]), float(params["V"]))
    elif kind_str == "OneCompOralFirstOrder":
        Kind = jl.OpenPKPDCore.OneCompOralFirstOrder
        Params = jl.OpenPKPDCore.OneCompOralFirstOrderParams
        model_params = Params(float(params["Ka"]), float(params["CL"]), float(params["V"]))
    elif kind_str == "TwoCompIVBolus":
        Kind = jl.OpenPKPDCore.TwoCompIVBolus
        Params = jl.OpenPKPDCore.TwoCompIVBolusParams
        model_params = Params(float(params["CL"]), float(params["V1"]), float(params["Q"]), float(params["V2"]))
    else:
        raise ValueError(f"Unsupported model kind: {kind_str}")

    base = ModelSpec(Kind(), "vpc_model", model_params, doses_vec)

    # Build IIV spec
    iiv_spec = spec.get("iiv", {})
    omegas = {jl.Symbol(k): float(v) for k, v in iiv_spec.get("omegas", {"CL": 0.3}).items()}
    n = iiv_spec.get("n", 100)
    seed = iiv_spec.get("seed", 12345)

    iiv = IIVSpec(LogNormalIIV(), omegas, jl.UInt64(seed), n)

    return PopulationSpec(base, iiv, None, None, [])


def _build_error_spec(jl, spec: Dict[str, Any]):
    """Build ResidualErrorSpec from Python dict."""
    error_type = spec.get("type", "proportional")
    sigma = spec.get("sigma", 0.1)

    if error_type == "additive":
        kind = jl.OpenPKPDCore.AdditiveError()
        params = jl.OpenPKPDCore.AdditiveErrorParams(float(sigma))
    elif error_type == "proportional":
        kind = jl.OpenPKPDCore.ProportionalError()
        params = jl.OpenPKPDCore.ProportionalErrorParams(float(sigma))
    elif error_type == "combined":
        kind = jl.OpenPKPDCore.CombinedError()
        params = jl.OpenPKPDCore.CombinedErrorParams(float(sigma), float(spec.get("sigma_prop", sigma)))
    else:
        kind = jl.OpenPKPDCore.ProportionalError()
        params = jl.OpenPKPDCore.ProportionalErrorParams(float(sigma))

    return jl.OpenPKPDCore.ResidualErrorSpec(kind, params, jl.Symbol("conc"), jl.UInt64(1))


def _convert_vpc_result(result, config: VPCConfig) -> VPCResult:
    """Convert Julia VPCResult to Python dataclass."""
    bins = []
    for bin in result.bins:
        percentiles = []
        for p in bin.percentiles:
            percentiles.append(VPCPercentile(
                level=float(p.level),
                observed=float(p.observed),
                simulated_median=float(p.simulated_median),
                simulated_lower=float(p.simulated_lower),
                simulated_upper=float(p.simulated_upper),
            ))

        bins.append(VPCBin(
            bin_id=int(bin.bin_id),
            time_min=float(bin.time_min),
            time_max=float(bin.time_max),
            time_midpoint=float(bin.time_midpoint),
            n_observed=int(bin.n_observed),
            percentiles=percentiles,
        ))

    return VPCResult(
        n_subjects_observed=int(result.n_subjects_observed),
        n_observations_observed=int(result.n_observations_observed),
        n_simulations=int(result.n_simulations),
        prediction_corrected=config.prediction_corrected,
        pi_levels=config.pi_levels,
        ci_level=config.ci_level,
        bins=bins,
    )


__all__ = [
    "compute_vpc",
    "compute_pcvpc",
    "VPCConfig",
    "VPCResult",
    "VPCBin",
    "VPCPercentile",
]
