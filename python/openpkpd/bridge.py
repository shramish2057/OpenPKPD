from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

_JL = None

@dataclass(frozen=True)
class RepoPaths:
    repo_root: Path
    core_project: Path

def version() -> str:
    jl = _require_julia()
    return str(jl.OpenPKPDCore.OPENPKPD_VERSION)


def _detect_repo_root(start: Optional[Path] = None) -> Path:
    here = (start or Path(__file__)).resolve()
    for p in [here] + list(here.parents):
        if (p / "core" / "OpenPKPDCore" / "Project.toml").exists():
            return p
    raise RuntimeError("Could not locate repo root (core/OpenPKPDCore/Project.toml not found).")


def init_julia(repo_root: Optional[Union[str, Path]] = None) -> None:
    """
    Initialize Julia and activate core/OpenPKPDCore deterministically.
    Safe to call multiple times.
    """
    global _JL
    if _JL is not None:
        return

    root = Path(repo_root).resolve() if repo_root else _detect_repo_root()
    core_project = root / "core" / "OpenPKPDCore"

    from juliacall import Main as jl  # type: ignore

    # Activate and instantiate the exact Julia project in this repo
    jl.seval("import Pkg")
    jl.Pkg.activate(str(core_project))
    jl.Pkg.instantiate()

    jl.seval("using OpenPKPDCore")

    _JL = jl


def _require_julia() -> Any:
    if _JL is None:
        init_julia()
    return _JL


def _simresult_to_py(res: Any) -> Dict[str, Any]:
    """
    Convert OpenPKPDCore.SimResult to Python dicts with lists of floats.
    """
    t = list(res.t)
    states = {str(k): list(v) for (k, v) in res.states.items()}
    obs = {str(k): list(v) for (k, v) in res.observations.items()}
    meta = dict(res.metadata)
    return {"t": t, "states": states, "observations": obs, "metadata": meta}


def _popresult_to_py(popres: Any) -> Dict[str, Any]:
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


def replay_artifact(path: Union[str, Path]) -> Dict[str, Any]:
    """
    Replay any artifact JSON (single, population, sensitivity) and return a Python dict.

    For single artifacts: returns SimResult dict.
    For population artifacts: returns PopulationResult dict.
    For sensitivity artifacts: returns dict with base and pert series and metrics.
    """
    jl = _require_julia()
    artifact = jl.OpenPKPDCore.read_execution_json(str(Path(path).resolve()))

    atype = "single"
    if "artifact_type" in artifact:
        atype = str(artifact["artifact_type"])

    if atype == "population":
        res = jl.OpenPKPDCore.replay_population_execution(artifact)
        return _popresult_to_py(res)

    if atype == "sensitivity_single":
        res = jl.OpenPKPDCore.replay_sensitivity_execution(artifact)
        return {
            "plan": {"name": str(res.plan.name)},
            "observation": str(res.observation),
            "base_series": list(res.base_metric_series),
            "pert_series": list(res.pert_metric_series),
            "metrics": {
                "max_abs_delta": float(res.metrics.max_abs_delta),
                "max_rel_delta": float(res.metrics.max_rel_delta),
                "l2_norm_delta": float(res.metrics.l2_norm_delta),
            },
            "metadata": dict(res.metadata),
        }

    if atype == "sensitivity_population":
        res = jl.OpenPKPDCore.replay_population_sensitivity_execution(artifact)
        return {
            "plan": {"name": str(res.plan.name)},
            "observation": str(res.observation),
            "probs": list(res.probs),
            "base_mean": list(res.base_summary_mean),
            "pert_mean": list(res.pert_summary_mean),
            "metrics_mean": {
                "max_abs_delta": float(res.metrics_mean.max_abs_delta),
                "max_rel_delta": float(res.metrics_mean.max_rel_delta),
                "l2_norm_delta": float(res.metrics_mean.l2_norm_delta),
            },
            "metadata": dict(res.metadata),
        }

    res = jl.OpenPKPDCore.replay_execution(artifact)
    return _simresult_to_py(res)

def _to_julia_vector(jl: Any, items: list, item_type: Any) -> Any:
    """Convert a Python list to a Julia Vector of the specified type."""
    vec = jl.Vector[item_type](jl.undef, len(items))
    for i, item in enumerate(items):
        vec[i] = item  # PythonCall uses 0-based indexing from Python
    return vec


def write_single_artifact(
    path: Union[str, Path],
    *,
    model: Dict[str, Any],
    grid: Dict[str, Any],
    solver: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Run a single PK simulation and write a Julia-native execution artifact.

    model:
      {
        "kind": "OneCompIVBolus" | "OneCompOralFirstOrder",
        "params": {...},
        "doses": [{"time": float, "amount": float}]
      }

    grid:
      {"t0": float, "t1": float, "saveat": [float]}

    solver:
      optional, defaults match core golden settings
    """
    jl = _require_julia()

    DoseEvent = jl.OpenPKPDCore.DoseEvent
    ModelSpec = jl.OpenPKPDCore.ModelSpec
    SimGrid = jl.OpenPKPDCore.SimGrid
    SolverSpec = jl.OpenPKPDCore.SolverSpec

    kind = model["kind"]
    dose_list = [DoseEvent(float(d["time"]), float(d["amount"])) for d in model["doses"]]
    doses = _to_julia_vector(jl, dose_list, DoseEvent)

    if kind == "OneCompIVBolus":
        OneCompIVBolus = jl.OpenPKPDCore.OneCompIVBolus
        Params = jl.OpenPKPDCore.OneCompIVBolusParams
        params = Params(float(model["params"]["CL"]), float(model["params"]["V"]))
        spec = ModelSpec(OneCompIVBolus(), "py_artifact", params, doses)
    elif kind == "OneCompOralFirstOrder":
        OneCompOralFirstOrder = jl.OpenPKPDCore.OneCompOralFirstOrder
        Params = jl.OpenPKPDCore.OneCompOralFirstOrderParams
        params = Params(
            float(model["params"]["KA"]),
            float(model["params"]["CL"]),
            float(model["params"]["V"]),
        )
        spec = ModelSpec(OneCompOralFirstOrder(), "py_artifact", params, doses)
    else:
        raise ValueError(f"Unsupported model kind: {kind}")

    grid_jl = SimGrid(
        float(grid["t0"]),
        float(grid["t1"]),
        [float(x) for x in grid["saveat"]],
    )

    if solver is None:
        solver_jl = SolverSpec(jl.Symbol("Tsit5"), 1e-10, 1e-12, 10**7)
    else:
        solver_jl = SolverSpec(
            jl.Symbol(solver.get("alg", "Tsit5")),
            float(solver.get("reltol", 1e-10)),
            float(solver.get("abstol", 1e-12)),
            int(solver.get("maxiters", 10**7)),
        )

    res = jl.OpenPKPDCore.simulate(spec, grid_jl, solver_jl)

    jl.OpenPKPDCore.write_execution_json(
        str(Path(path).resolve()),
        model_spec=spec,
        grid=grid_jl,
        solver=solver_jl,
        result=res,
    )


def simulate_pk_iv_bolus(
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    jl = _require_julia()

    DoseEvent = jl.OpenPKPDCore.DoseEvent
    ModelSpec = jl.OpenPKPDCore.ModelSpec
    OneCompIVBolus = jl.OpenPKPDCore.OneCompIVBolus
    OneCompIVBolusParams = jl.OpenPKPDCore.OneCompIVBolusParams
    SimGrid = jl.OpenPKPDCore.SimGrid
    SolverSpec = jl.OpenPKPDCore.SolverSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    spec = ModelSpec(OneCompIVBolus(), "py_iv_bolus", OneCompIVBolusParams(float(cl), float(v)), doses_vec)
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.OpenPKPDCore.simulate(spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pk_oral_first_order(
    ka: float,
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    jl = _require_julia()

    DoseEvent = jl.OpenPKPDCore.DoseEvent
    ModelSpec = jl.OpenPKPDCore.ModelSpec
    OneCompOralFirstOrder = jl.OpenPKPDCore.OneCompOralFirstOrder
    OneCompOralFirstOrderParams = jl.OpenPKPDCore.OneCompOralFirstOrderParams
    SimGrid = jl.OpenPKPDCore.SimGrid
    SolverSpec = jl.OpenPKPDCore.SolverSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "py_oral_first_order",
        OneCompOralFirstOrderParams(float(ka), float(cl), float(v)),
        doses_vec,
    )
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.OpenPKPDCore.simulate(spec, grid, solver)
    return _simresult_to_py(res)


def simulate_population_iv_bolus(
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    n: int,
    seed: int,
    omegas: Dict[str, float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    jl = _require_julia()

    DoseEvent = jl.OpenPKPDCore.DoseEvent
    ModelSpec = jl.OpenPKPDCore.ModelSpec
    OneCompIVBolus = jl.OpenPKPDCore.OneCompIVBolus
    OneCompIVBolusParams = jl.OpenPKPDCore.OneCompIVBolusParams
    SimGrid = jl.OpenPKPDCore.SimGrid
    SolverSpec = jl.OpenPKPDCore.SolverSpec
    IIVSpec = jl.OpenPKPDCore.IIVSpec
    LogNormalIIV = jl.OpenPKPDCore.LogNormalIIV
    PopulationSpec = jl.OpenPKPDCore.PopulationSpec
    IndividualCovariates = jl.OpenPKPDCore.IndividualCovariates

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    base = ModelSpec(OneCompIVBolus(), "py_pop_iv", OneCompIVBolusParams(float(cl), float(v)), doses_vec)

    omega_j = {jl.Symbol(k): float(val) for (k, val) in omegas.items()}
    iiv = IIVSpec(LogNormalIIV(), omega_j, jl.UInt64(int(seed)), int(n))

    # No IOV, no covariate model, no covariates
    pop = PopulationSpec(base, iiv, None, None, [])

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.OpenPKPDCore.simulate_population(pop, grid, solver)
    return _popresult_to_py(res)
