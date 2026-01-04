"""
OpenPKPD Core - Internal module for Julia bridge and utilities.

This module provides the core infrastructure for interfacing with the
OpenPKPD Julia core, including initialization, type conversion, and
data classes.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

_JL = None


# ============================================================================
# Data Classes
# ============================================================================

@dataclass(frozen=True)
class RepoPaths:
    """Repository path configuration."""
    repo_root: Path
    core_project: Path


@dataclass(frozen=True)
class SensitivityMetrics:
    """Results from sensitivity analysis."""
    max_abs_delta: float
    max_rel_delta: float
    l2_norm_delta: float


@dataclass(frozen=True)
class SensitivityResult:
    """Full sensitivity analysis result."""
    plan_name: str
    observation: str
    base_series: List[float]
    pert_series: List[float]
    metrics: SensitivityMetrics
    metadata: Dict[str, Any]


@dataclass(frozen=True)
class PopulationSensitivityResult:
    """Population sensitivity analysis result."""
    plan_name: str
    observation: str
    probs: List[float]
    base_mean: List[float]
    pert_mean: List[float]
    metrics_mean: SensitivityMetrics
    metadata: Dict[str, Any]


# ============================================================================
# Initialization
# ============================================================================

def version() -> str:
    """
    Get the OpenPKPD version string.

    Returns:
        str: The version string (e.g., "0.1.0")

    Example:
        >>> openpkpd.version()
        '0.1.0'
    """
    jl = _require_julia()
    return str(jl.OpenPKPDCore.OPENPKPD_VERSION)


def _detect_repo_root(start: Optional[Path] = None) -> Path:
    """Detect the OpenPKPD repository root directory."""
    here = (start or Path(__file__)).resolve()
    for p in [here] + list(here.parents):
        if (p / "packages" / "core" / "Project.toml").exists():
            return p
    raise RuntimeError("Could not locate repo root (packages/core/Project.toml not found).")


def init_julia(repo_root: Optional[Union[str, Path]] = None) -> None:
    """
    Initialize the Julia runtime and load OpenPKPDCore.

    This function activates the Julia project and loads the core simulation
    engine. It's safe to call multiple times - subsequent calls are no-ops.

    Args:
        repo_root: Optional path to the OpenPKPD repository root.
                   If not provided, auto-detection is used.

    Example:
        >>> import openpkpd
        >>> openpkpd.init_julia()
        >>> # Now ready to run simulations
    """
    global _JL
    if _JL is not None:
        return

    root = Path(repo_root).resolve() if repo_root else _detect_repo_root()
    core_project = root / "packages" / "core"

    from juliacall import Main as jl  # type: ignore

    # Activate and instantiate the exact Julia project in this repo
    jl.seval("import Pkg")
    jl.Pkg.activate(str(core_project))
    jl.Pkg.instantiate()

    jl.seval("using OpenPKPDCore")

    _JL = jl


def _require_julia() -> Any:
    """Ensure Julia is initialized and return the Julia main module."""
    if _JL is None:
        init_julia()
    return _JL


# ============================================================================
# Result Conversion Utilities
# ============================================================================

def _simresult_to_py(res: Any) -> Dict[str, Any]:
    """
    Convert OpenPKPDCore.SimResult to a Python dictionary.

    Args:
        res: Julia SimResult object

    Returns:
        Dict with keys: t, states, observations, metadata
    """
    t = list(res.t)
    states = {str(k): list(v) for (k, v) in res.states.items()}
    obs = {str(k): list(v) for (k, v) in res.observations.items()}
    meta = dict(res.metadata)
    return {"t": t, "states": states, "observations": obs, "metadata": meta}


def _popresult_to_py(popres: Any) -> Dict[str, Any]:
    """
    Convert OpenPKPDCore.PopulationResult to a Python dictionary.

    Args:
        popres: Julia PopulationResult object

    Returns:
        Dict with keys: individuals, params, summaries, metadata
    """
    individuals = [_simresult_to_py(r) for r in popres.individuals]
    params = [{str(k): float(v) for (k, v) in d.items()} for d in popres.params]
    summaries = {}

    for (k, s) in popres.summaries.items():
        summaries[str(k)] = {
            "observation": str(s.observation),
            "probs": list(s.probs),
            "mean": list(s.mean),
            "median": list(s.median),
            "quantiles": {str(p): list(v) for (p, v) in s.quantiles.items()},
        }

    meta = dict(popres.metadata)
    return {"individuals": individuals, "params": params, "summaries": summaries, "metadata": meta}


def _to_julia_vector(jl: Any, items: list, item_type: Any) -> Any:
    """Convert a Python list to a Julia Vector of the specified type."""
    vec = jl.Vector[item_type](jl.undef, len(items))
    for i, item in enumerate(items):
        vec[i] = item  # PythonCall uses 0-based indexing from Python
    return vec


def _to_julia_float_vector(jl: Any, items: list) -> Any:
    """Convert a Python list to a Julia Vector{Float64}."""
    vec = jl.Vector[jl.Float64](jl.undef, len(items))
    for i, item in enumerate(items):
        vec[i] = float(item)
    return vec


def _create_dose_event(jl: Any, dose_dict: Dict[str, Any]) -> Any:
    """
    Create a Julia DoseEvent from a Python dict.

    Supports both bolus and infusion administration:
    - Bolus: {"time": 0.0, "amount": 100.0}  (duration defaults to 0.0)
    - Infusion: {"time": 0.0, "amount": 100.0, "duration": 1.0}  (100 mg over 1 hour)

    Args:
        jl: Julia module
        dose_dict: Dict with 'time', 'amount', and optional 'duration'

    Returns:
        Julia DoseEvent object
    """
    DoseEvent = jl.OpenPKPDCore.DoseEvent
    return DoseEvent(
        float(dose_dict["time"]),
        float(dose_dict["amount"]),
        float(dose_dict.get("duration", 0.0))
    )


def _create_dose_events(jl: Any, doses: List[Dict[str, Any]]) -> Any:
    """
    Create a Julia Vector of DoseEvents from a list of Python dicts.

    Args:
        jl: Julia module
        doses: List of dose dicts with 'time', 'amount', and optional 'duration'

    Returns:
        Julia Vector{DoseEvent}
    """
    DoseEvent = jl.OpenPKPDCore.DoseEvent
    dose_objs = [_create_dose_event(jl, d) for d in doses]
    return _to_julia_vector(jl, dose_objs, DoseEvent)
