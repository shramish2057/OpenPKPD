# Model Import

OpenPKPD can import model definitions from NONMEM control files (.ctl) and Monolix project files (.mlxtran), enabling migration of existing models to OpenPKPD.

## NONMEM Import

### Supported Features

| Feature | Support |
|---------|---------|
| ADVAN 1-4, 11 | Full |
| TRANS 1, 2, 4 | Full |
| $THETA | Full |
| $OMEGA (DIAGONAL) | Full |
| $OMEGA (BLOCK) | Full |
| $SIGMA | Full |
| IIV (exponential) | Full |
| Simple covariates | Partial |
| $PK user code | Limited |
| $ERROR user code | Limited |

### Basic Usage

```julia
using OpenPKPDCore

# Parse NONMEM control file
nmctl = parse_nonmem_control("run001.ctl")

# Convert to OpenPKPD format
model_spec, pop_spec, mapping = convert_nonmem_to_openpkpd(nmctl)

# Use for simulation
result = simulate(model_spec, grid, solver)
```

```python
from openpkpd.import_ import parse_nonmem, convert_to_openpkpd

# Parse and convert
nmctl = parse_nonmem("run001.ctl")
model_spec, pop_spec, mapping = convert_to_openpkpd(nmctl)
```

### ADVAN/TRANS Mapping

| ADVAN | TRANS | OpenPKPD Model | Parameters |
|-------|-------|----------------|------------|
| 1 | 2 | OneCompIVBolus | CL, V |
| 2 | 2 | OneCompOralFirstOrder | Ka, CL, V |
| 3 | 4 | TwoCompIVBolus | CL, V1, Q, V2 |
| 4 | 4 | TwoCompOral | Ka, CL, V1, Q, V2 |
| 11 | 4 | ThreeCompIVBolus | CL, V1, Q2, V2, Q3, V3 |

### Example NONMEM Control File

```
$PROBLEM One-Compartment PK Model

$DATA data.csv IGNORE=@
$INPUT ID TIME DV AMT EVID MDV WT

$SUBROUTINES ADVAN2 TRANS2

$PK
TVCL = THETA(1)
TVV  = THETA(2)
TVKA = THETA(3)

CL = TVCL * EXP(ETA(1))
V  = TVV  * EXP(ETA(2))
KA = TVKA * EXP(ETA(3))

$ERROR
IPRED = F
Y = F * (1 + ERR(1))

$THETA
(0, 5, 100)    ; CL
(0, 50, 500)   ; V
(0, 1.5, 10)   ; KA

$OMEGA DIAGONAL
0.09  ; IIV CL
0.04  ; IIV V
0.16  ; IIV KA

$SIGMA
0.01  ; Proportional error

$EST METHOD=1 INTER MAXEVAL=9999
```

### Parsed Structure

```julia
nmctl = parse_nonmem_control("run001.ctl")

# Access parsed components
println(nmctl.problem)        # "One-Compartment PK Model"
println(nmctl.advan)          # 2
println(nmctl.trans)          # 2
println(nmctl.theta)          # Vector of THETASpec
println(nmctl.omega)          # OMEGABlock
println(nmctl.sigma)          # SIGMABlock
```

### Parameter Mapping

The conversion returns a mapping dictionary:

```julia
model_spec, pop_spec, mapping = convert_nonmem_to_openpkpd(nmctl)

# mapping contains:
# {
#   :THETA1 => :CL,
#   :THETA2 => :V,
#   :THETA3 => :Ka,
#   :ETA1 => :CL,
#   :ETA2 => :V,
#   :ETA3 => :Ka
# }
```

## Monolix Import

### Supported Features

| Feature | Support |
|---------|---------|
| pk= macros | Full |
| compartmental_model | Full |
| parameter blocks | Full |
| [INDIVIDUAL] | Full |
| omega (diagonal) | Full |
| omega (block) | Partial |
| Covariates | Partial |
| Custom models | Limited |

### Basic Usage

```julia
using OpenPKPDCore

# Parse Monolix project
mlx = parse_monolix_project("project.mlxtran")

# Convert to OpenPKPD format
model_spec, pop_spec, mapping = convert_monolix_to_openpkpd(mlx)
```

```python
from openpkpd.import_ import parse_monolix, convert_to_openpkpd

mlx = parse_monolix("project.mlxtran")
model_spec, pop_spec, mapping = convert_to_openpkpd(mlx)
```

### Example Monolix Project

```
<DATAFILE>
[FILEINFO]
file = 'data.csv'
delimiter = comma
header = {ID, TIME, DV, AMT, EVID, MDV, WT}

[CONTENT]
ID = {use=identifier}
TIME = {use=time}
DV = {use=observation, name=Y, type=continuous}
AMT = {use=amount}
EVID = {use=eventidentifier}

<MODEL>
[LONGITUDINAL]
input = {Cl, V, ka}

pk=compartmental_model(
  model_type = 1-cpt,
  admin = oral,
  distribution = one-compartment,
  elimination = linear,
  p_ka = ka,
  p_Cl = Cl,
  p_V = V
)

DEFINITION:
Y = {distribution=normal, prediction=Cc, errorModel=proportional1}

<PARAMETER>
Cl = {value=5, method=MLE}
V = {value=50, method=MLE}
ka = {value=1.5, method=MLE}

<MONOLIX>
[TASKS]
globalSettings={}
estimatePopulationParameters( method = SAEM )

[SETTINGS]
GLOBAL:
seed = 12345
```

### pk= Macro Mapping

| Monolix Macro | OpenPKPD Model |
|---------------|----------------|
| `model_type=1-cpt, admin=iv` | OneCompIVBolus |
| `model_type=1-cpt, admin=oral` | OneCompOralFirstOrder |
| `model_type=2-cpt, admin=iv` | TwoCompIVBolus |
| `model_type=2-cpt, admin=oral` | TwoCompOral |
| `elimination=Michaelis-Menten` | MichaelisMentenElimination |

## CLI Usage

### Import NONMEM Model

```bash
./packages/cli/bin/openpkpd import \
    --input run001.ctl \
    --format nonmem \
    --out model.json
```

### Import Monolix Model

```bash
./packages/cli/bin/openpkpd import \
    --input project.mlxtran \
    --format monolix \
    --out model.json
```

### Output Format

The output JSON contains the converted model specification:

```json
{
  "model_spec": {
    "kind": "OneCompOralFirstOrder",
    "params": {
      "Ka": 1.5,
      "CL": 5.0,
      "V": 50.0
    },
    "doses": []
  },
  "population_spec": {
    "iiv": {
      "kind": "LogNormalIIV",
      "omegas": {
        "Ka": 0.4,
        "CL": 0.3,
        "V": 0.2
      }
    },
    "error_spec": {
      "kind": "proportional",
      "sigma": 0.1
    }
  },
  "mapping": {
    "THETA1": "CL",
    "THETA2": "V",
    "THETA3": "Ka"
  },
  "source": {
    "format": "nonmem",
    "file": "run001.ctl"
  }
}
```

## Workflow Examples

### NONMEM to OpenPKPD Estimation

```python
from openpkpd.import_ import parse_nonmem
from openpkpd.data import read_cdisc_csv, cdisc_to_population
from openpkpd.estimation import estimate

# 1. Import NONMEM model
nmctl = parse_nonmem("run001.ctl")
model_spec, _, mapping = convert_to_openpkpd(nmctl)

# 2. Read data
dataset = read_cdisc_csv("pc.csv", "ex.csv", "dm.csv")
pop_spec, observed = cdisc_to_population(dataset, model_spec)

# 3. Re-estimate in OpenPKPD
result = estimate(observed, model_spec, config)

# 4. Compare parameters
print("NONMEM THETA:", nmctl.theta)
print("OpenPKPD Theta:", result.theta)
```

### Model Translation Validation

```python
# Run simulation with both tools and compare
import numpy as np

# OpenPKPD simulation
openpkpd_result = simulate(model_spec, grid, solver)

# Compare to NONMEM output
nonmem_pred = np.loadtxt("nonmem_pred.csv")
openpkpd_pred = np.array(openpkpd_result["observations"]["conc"])

# Check agreement
max_diff = np.max(np.abs(nonmem_pred - openpkpd_pred))
print(f"Max difference: {max_diff:.6f}")
```

## Limitations

### NONMEM

| Feature | Status |
|---------|--------|
| Complex $PK code | Manual translation needed |
| TIME-VARYING covariates in $PK | Limited support |
| $DES (custom ODEs) | Not supported |
| PREDPP subroutines | ADVAN 1-4, 11 only |
| $MIX (mixture models) | Not supported |

### Monolix

| Feature | Status |
|---------|--------|
| Custom structural models | Manual translation needed |
| Complex covariate models | Limited support |
| Occasion effects (IOV) | Partial support |
| Mixture models | Not supported |

## Troubleshooting

### Unsupported ADVAN

```
Error: ADVAN 6 not supported
```

ADVAN 6 uses custom ODEs ($DES). You'll need to implement the model manually in OpenPKPD.

### Parameter Mismatch

```
Warning: Could not map THETA(4)
```

Some parameters in complex $PK blocks may not map directly. Check the mapping dictionary and manually adjust if needed.

### Covariate Model

```
Warning: Covariate model not fully translated
```

Complex covariate relationships (e.g., piecewise, categorical) may need manual implementation.

## See Also

- [Data Import](data.md) - CDISC data format
- [Parameter Estimation](estimation.md) - Re-estimating imported models
- [Models](models.md) - Available OpenPKPD models
