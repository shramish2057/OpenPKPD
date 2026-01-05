"""
OpenPKPD Model Import Module

This module provides Python bindings for importing models from
NONMEM control files (.ctl) and Monolix project files (.mlxtran).
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Union
from pathlib import Path

from .._core import _require_julia


@dataclass
class ImportedModel:
    """Result from importing an external model file."""
    source_format: str  # "nonmem" or "monolix"
    source_file: str
    model_kind: str
    params: Dict[str, float]
    theta_init: List[float]
    theta_names: List[str]
    omega_init: List[List[float]]
    omega_names: List[str]
    sigma_type: str
    sigma_init: float
    warnings: List[str]
    metadata: Dict[str, Any]


@dataclass
class NONMEMControlFile:
    """Parsed NONMEM control file."""
    problem: str
    data_file: Optional[str]
    input_columns: List[str]
    advan: int
    trans: int
    theta_specs: List[Dict[str, Any]]
    omega_block: List[List[float]]
    sigma_block: List[List[float]]
    pk_code: List[str]
    error_code: List[str]


@dataclass
class MonolixProject:
    """Parsed Monolix project."""
    model_type: str
    structural_model: str
    parameters: Dict[str, Any]
    data_file: Optional[str]


def import_nonmem(path: Union[str, Path], doses: Optional[List[Dict[str, float]]] = None) -> ImportedModel:
    """
    Import a model from a NONMEM control file.

    Supports ADVAN1-4, ADVAN10, ADVAN11 with common TRANS options. Extracts:
    - Model structure (ADVAN/TRANS)
    - THETA initial estimates and bounds
    - OMEGA structure and values
    - SIGMA values

    Args:
        path: Path to NONMEM control file (.ctl or .mod)
        doses: Optional list of dose events, e.g., [{"time": 0.0, "amount": 100.0}]

    Returns:
        ImportedModel with OpenPKPD-compatible model specification

    Supported ADVAN/TRANS combinations:
        - ADVAN1/TRANS2: One-compartment IV bolus (CL, V)
        - ADVAN2/TRANS2: One-compartment oral (Ka, CL, V)
        - ADVAN3/TRANS4: Two-compartment IV bolus (CL, V1, Q, V2)
        - ADVAN4/TRANS4: Two-compartment oral (Ka, CL, V1, Q, V2)
        - ADVAN10: Michaelis-Menten elimination
        - ADVAN11: Three-compartment model

    Example:
        >>> model = import_nonmem("run001.ctl", doses=[{"time": 0.0, "amount": 100.0}])
        >>> print(f"Model type: {model.model_kind}")
        >>> print(f"Parameters: {model.params}")
    """
    jl = _require_julia()

    # Read the control file text
    path = Path(path).resolve()
    with open(path, 'r') as f:
        ctl_text = f.read()

    # Parse the control file
    ctl = jl.OpenPKPDCore.parse_nonmem_control(ctl_text)

    # Create dose events if provided
    jl_doses = jl.OpenPKPDCore.DoseEvent[]
    if doses:
        for d in doses:
            dose = jl.OpenPKPDCore.DoseEvent(float(d.get("time", 0.0)), float(d.get("amount", 100.0)))
            jl_doses = jl.seval("push!")([jl_doses, dose])

    # Convert to OpenPKPD format
    result = jl.OpenPKPDCore.convert_nonmem_to_openpkpd(ctl, doses=jl_doses)

    return _convert_julia_import_result(result, "nonmem", str(path))


def import_monolix(path: Union[str, Path], doses: Optional[List[Dict[str, float]]] = None) -> ImportedModel:
    """
    Import a model from a Monolix project file.

    Parses .mlxtran files and extracts:
    - Structural model definition
    - Parameter initial values
    - Random effect structure

    Args:
        path: Path to Monolix project file (.mlxtran)
        doses: Optional list of dose events, e.g., [{"time": 0.0, "amount": 100.0}]

    Returns:
        ImportedModel with OpenPKPD-compatible model specification

    Example:
        >>> model = import_monolix("project.mlxtran", doses=[{"time": 0.0, "amount": 100.0}])
        >>> print(f"Model type: {model.model_kind}")
    """
    jl = _require_julia()

    # Read the project file text
    path = Path(path).resolve()
    with open(path, 'r') as f:
        mlx_text = f.read()

    # Parse the project file
    project = jl.OpenPKPDCore.parse_monolix_project(mlx_text)

    # Create dose events if provided
    jl_doses = jl.OpenPKPDCore.DoseEvent[]
    if doses:
        for d in doses:
            dose = jl.OpenPKPDCore.DoseEvent(float(d.get("time", 0.0)), float(d.get("amount", 100.0)))
            jl_doses = jl.seval("push!")([jl_doses, dose])

    # Convert to OpenPKPD format
    result = jl.OpenPKPDCore.convert_monolix_to_openpkpd(project, doses=jl_doses)

    return _convert_julia_import_result(result, "monolix", str(path))


def import_model(path: Union[str, Path], format: Optional[str] = None) -> ImportedModel:
    """
    Import a model from NONMEM or Monolix format.

    Auto-detects format based on file extension if not specified.

    Args:
        path: Path to model file
        format: Optional format override ("nonmem" or "monolix")

    Returns:
        ImportedModel with OpenPKPD-compatible model specification

    Example:
        >>> model = import_model("run001.ctl")
        >>> model = import_model("project.mlxtran")
    """
    path = Path(path)

    if format is None:
        if path.suffix in [".ctl", ".mod"]:
            format = "nonmem"
        elif path.suffix == ".mlxtran":
            format = "monolix"
        else:
            raise ValueError(f"Cannot auto-detect format for {path.suffix}. Specify format='nonmem' or format='monolix'")

    if format == "nonmem":
        return import_nonmem(path)
    elif format == "monolix":
        return import_monolix(path)
    else:
        raise ValueError(f"Unknown format: {format}")


def parse_nonmem_control(path: Union[str, Path]) -> NONMEMControlFile:
    """
    Parse a NONMEM control file without converting to OpenPKPD format.

    Useful for inspecting the control file structure.

    Args:
        path: Path to NONMEM control file

    Returns:
        NONMEMControlFile with parsed sections

    Example:
        >>> ctl = parse_nonmem_control("run001.ctl")
        >>> print(f"Problem: {ctl.problem}")
        >>> print(f"ADVAN{ctl.advan} TRANS{ctl.trans}")
    """
    jl = _require_julia()
    result = jl.OpenPKPDCore.parse_nonmem_control(str(Path(path).resolve()))
    return _convert_nonmem_control(result)


def parse_monolix_project(path: Union[str, Path]) -> MonolixProject:
    """
    Parse a Monolix project file without converting to OpenPKPD format.

    Args:
        path: Path to Monolix project file

    Returns:
        MonolixProject with parsed structure

    Example:
        >>> mlx = parse_monolix_project("project.mlxtran")
        >>> print(f"Model: {mlx.structural_model}")
    """
    jl = _require_julia()
    result = jl.OpenPKPDCore.parse_monolix_project(str(Path(path).resolve()))
    return _convert_monolix_project(result)


def _convert_julia_import_result(result, format: str, path: str) -> ImportedModel:
    """Convert Julia NONMEMConversionResult or MonolixConversionResult to Python ImportedModel."""
    # Check for conversion errors
    if result.errors and len(result.errors) > 0:
        raise ValueError(f"Import failed: {'; '.join(str(e) for e in result.errors)}")

    model_spec = result.model_spec
    if model_spec is None:
        raise ValueError("Import failed: no model spec generated")

    # Extract model kind from model type
    model = model_spec.model
    model_kind = type(model).__name__

    # Extract params as dict
    params = {}
    model_params = model_spec.params
    # Get field names from the Julia struct
    try:
        for field in model_params._jl_field_names:
            val = getattr(model_params, str(field))
            if isinstance(val, (int, float)):
                params[str(field)] = float(val)
    except AttributeError:
        # Fallback: try introspection
        for attr in dir(model_params):
            if not attr.startswith('_'):
                try:
                    val = getattr(model_params, attr)
                    if isinstance(val, (int, float)):
                        params[attr] = float(val)
                except Exception:
                    pass

    # Extract theta info
    theta_init = list(params.values())
    theta_names = list(params.keys())

    # Extract omega from IIV spec
    omega_init = []
    omega_names = []
    iiv_spec = result.iiv_spec
    if iiv_spec is not None:
        try:
            omegas = iiv_spec.omegas
            for k, v in omegas.items():
                omega_names.append(str(k))
                omega_init.append([float(v)])
        except Exception:
            pass

    # Extract sigma from error spec
    sigma_type = "proportional"
    sigma_init = 0.1
    error_spec = result.error_spec
    if error_spec is not None:
        try:
            sigma_type = str(error_spec.kind) if hasattr(error_spec, 'kind') else "proportional"
            sigma_init = float(error_spec.sigma) if hasattr(error_spec, 'sigma') else 0.1
        except Exception:
            pass

    # Get warnings
    warnings = [str(w) for w in result.warnings] if result.warnings else []

    return ImportedModel(
        source_format=format,
        source_file=path,
        model_kind=model_kind,
        params=params,
        theta_init=theta_init,
        theta_names=theta_names,
        omega_init=omega_init if omega_init else [[0.09]],
        omega_names=omega_names if omega_names else ["eta_1"],
        sigma_type=sigma_type,
        sigma_init=sigma_init,
        warnings=warnings,
        metadata={"parameter_mapping": dict(result.parameter_mapping) if hasattr(result, 'parameter_mapping') else {}},
    )


def _convert_imported_model(raw_result, model_spec, pop_spec, metadata, format: str, path: str) -> ImportedModel:
    """Convert Julia import result to Python ImportedModel (legacy fallback)."""
    # Extract model kind
    model_kind = str(type(model_spec.kind).__name__) if hasattr(model_spec, 'kind') else "unknown"

    # Extract params as dict
    params = {}
    for field in dir(model_spec.params):
        if not field.startswith("_"):
            val = getattr(model_spec.params, field)
            if isinstance(val, (int, float)):
                params[field] = float(val)

    # Extract theta info
    theta_init = list(params.values())
    theta_names = list(params.keys())

    # Extract omega
    omega_init = []
    omega_names = []
    if pop_spec is not None and pop_spec.iiv is not None:
        omegas = pop_spec.iiv.omegas
        for k, v in omegas.items():
            omega_names.append(str(k))
            omega_init.append([float(v)])

    # Determine sigma type (default to proportional)
    sigma_type = "proportional"
    sigma_init = 0.1

    warnings = list(metadata.get("warnings", [])) if metadata else []

    return ImportedModel(
        source_format=format,
        source_file=path,
        model_kind=model_kind,
        params=params,
        theta_init=theta_init,
        theta_names=theta_names,
        omega_init=omega_init if omega_init else [[0.09]],
        omega_names=omega_names if omega_names else ["eta_1"],
        sigma_type=sigma_type,
        sigma_init=sigma_init,
        warnings=warnings,
        metadata=dict(metadata) if metadata else {},
    )


def _convert_nonmem_control(result) -> NONMEMControlFile:
    """Convert Julia NONMEM control file to Python dataclass."""
    theta_specs = []
    for spec in result.theta_specs:
        theta_specs.append({
            "init": float(spec.init),
            "lower": float(spec.lower) if spec.lower is not None else None,
            "upper": float(spec.upper) if spec.upper is not None else None,
            "fixed": bool(spec.fixed),
        })

    omega_block = [[float(x) for x in row] for row in result.omega_block]
    sigma_block = [[float(x) for x in row] for row in result.sigma_block]

    return NONMEMControlFile(
        problem=str(result.problem),
        data_file=str(result.data_file) if result.data_file else None,
        input_columns=[str(c) for c in result.input_columns],
        advan=int(result.advan),
        trans=int(result.trans),
        theta_specs=theta_specs,
        omega_block=omega_block,
        sigma_block=sigma_block,
        pk_code=[str(line) for line in result.pk_code],
        error_code=[str(line) for line in result.error_code],
    )


def _convert_monolix_project(result) -> MonolixProject:
    """Convert Julia Monolix project to Python dataclass."""
    return MonolixProject(
        model_type=str(result.model_type) if hasattr(result, "model_type") else "unknown",
        structural_model=str(result.structural_model) if hasattr(result, "structural_model") else "",
        parameters=dict(result.parameters) if hasattr(result, "parameters") else {},
        data_file=str(result.data_file) if hasattr(result, "data_file") and result.data_file else None,
    )


__all__ = [
    "import_nonmem",
    "import_monolix",
    "import_model",
    "parse_nonmem_control",
    "parse_monolix_project",
    "ImportedModel",
    "NONMEMControlFile",
    "MonolixProject",
]
