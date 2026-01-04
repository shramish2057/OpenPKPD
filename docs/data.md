# Data Import

OpenPKPD supports importing clinical pharmacokinetic data from CDISC/SDTM format, enabling analysis of real-world study data.

## CDISC/SDTM Format

CDISC (Clinical Data Interchange Standards Consortium) defines standard formats for clinical trial data. OpenPKPD reads the following domains:

| Domain | Description | Key Variables |
|--------|-------------|---------------|
| **PC** | Pharmacokinetic Concentrations | USUBJID, PCTESTCD, PCSTRESN, PCTPTNUM |
| **EX** | Exposure (Dosing) | USUBJID, EXDOSE, EXROUTE, EXSTDTC, EXDUR |
| **DM** | Demographics | USUBJID, AGE, SEX, RACE, WEIGHT |
| **PP** | PK Parameters (optional) | USUBJID, PPTESTCD, PPSTRESN |

## Reading CDISC Data

### CSV Format

```julia
using OpenPKPDCore

# Read CDISC domains from CSV files
dataset = read_cdisc_csv(
    "pc.csv",   # Concentrations
    "ex.csv",   # Dosing
    "dm.csv"    # Demographics (optional)
)
```

```python
from openpkpd.data import read_cdisc_csv

dataset = read_cdisc_csv("pc.csv", "ex.csv", "dm.csv")
```

### SAS Transport (XPT) Format

```julia
# Read SAS transport files
dataset = read_cdisc_xpt(
    "pc.xpt",
    "ex.xpt",
    "dm.xpt"
)
```

```python
from openpkpd.data import read_cdisc_xpt

dataset = read_cdisc_xpt("pc.xpt", "ex.xpt", "dm.xpt")
```

## Data Validation

Validate your CDISC data against standards:

```julia
warnings = validate_cdisc_dataset(dataset)

for w in warnings
    println("Warning: ", w)
end
```

```python
from openpkpd.data import validate_cdisc_dataset

warnings = validate_cdisc_dataset(dataset)
for w in warnings:
    print(f"Warning: {w}")
```

### Validation Checks

| Check | Description |
|-------|-------------|
| Required variables | Ensures mandatory fields are present |
| Data types | Validates numeric/string field types |
| Subject consistency | Cross-checks USUBJID across domains |
| Time ordering | Verifies chronological order |
| Dose-observation matching | Ensures subjects with PK have dosing |

## Converting to OpenPKPD Format

Convert CDISC data to OpenPKPD's internal format:

```julia
# Convert to PopulationSpec and ObservedData
pop_spec, observed = cdisc_to_population(
    dataset,
    model_spec;
    covariate_mapping=Dict(
        :WT => :WEIGHT,
        :CRCL => :CREATININE_CL
    )
)
```

```python
from openpkpd.data import cdisc_to_population

pop_spec, observed = cdisc_to_population(
    dataset,
    model_spec,
    covariate_mapping={"WT": "WEIGHT", "CRCL": "CREATININE_CL"}
)
```

### ObservedData Structure

| Field | Description |
|-------|-------------|
| `subject_ids` | Unique subject identifiers |
| `times` | Observation times |
| `observations` | Observed concentration values |
| `dvid` | Dependent variable ID (for multiple endpoints) |
| `doses` | Dose events extracted from EX domain |
| `covariates` | Subject-level covariates from DM |

## Handling Infusions

The EX domain supports infusion duration via the EXDUR field (ISO 8601 duration):

```csv
USUBJID,EXDOSE,EXROUTE,EXSTDTC,EXDUR
001,100,IV,2024-01-01T08:00,PT1H
002,100,IV,2024-01-01T08:00,PT30M
```

Duration codes:
- `PT1H` = 1 hour
- `PT30M` = 30 minutes
- `PT2H30M` = 2.5 hours

OpenPKPD automatically parses these into DoseEvent duration fields.

## PC Domain Format

Example PC domain CSV:

```csv
USUBJID,PCTESTCD,PCTEST,PCSTRESN,PCSTRESU,PCTPTNUM,PCTPT,PCSPEC
001,DRUG,Drug Concentration,2.5,ng/mL,1,1 hour post-dose,PLASMA
001,DRUG,Drug Concentration,1.8,ng/mL,2,2 hours post-dose,PLASMA
001,DRUG,Drug Concentration,1.2,ng/mL,4,4 hours post-dose,PLASMA
```

Key variables:
- `PCSTRESN`: Numeric result (concentration)
- `PCTPTNUM`: Numeric time point
- `PCSPEC`: Specimen type (PLASMA, SERUM, etc.)

## EX Domain Format

Example EX domain CSV:

```csv
USUBJID,EXDOSE,EXDOSU,EXROUTE,EXSTDTC,EXENDTC,EXDUR
001,100,mg,IV,2024-01-01T08:00,2024-01-01T09:00,PT1H
002,100,mg,ORAL,2024-01-01T08:00,,
```

Key variables:
- `EXDOSE`: Dose amount
- `EXROUTE`: Route of administration
- `EXDUR`: Duration for infusions (ISO 8601)

## DM Domain Format

Example DM domain CSV:

```csv
USUBJID,AGE,AGEU,SEX,RACE,WEIGHT,WEIGHTU,HEIGHT,HEIGHTU
001,45,YEARS,M,WHITE,80,kg,175,cm
002,52,YEARS,F,ASIAN,65,kg,162,cm
```

## BLQ Handling

Handle below-limit-of-quantification (BLQ) data:

```python
from openpkpd.data import read_cdisc_csv

dataset = read_cdisc_csv(
    "pc.csv", "ex.csv", "dm.csv",
    blq_handling="lloq_half"  # Replace BLQ with LLOQ/2
)
```

BLQ handling options:
- `"zero"`: Set BLQ values to 0
- `"missing"`: Exclude BLQ observations
- `"lloq_half"`: Set BLQ to LLOQ/2
- `"m3"`: Use M3 method (likelihood-based)

## Multiple Analytes

Handle studies with multiple analytes (parent + metabolite):

```python
# Filter by analyte
parent_data = dataset.filter_analyte("DRUG")
metabolite_data = dataset.filter_analyte("METAB")

# Or use DVID for multiple endpoints
pop_spec, observed = cdisc_to_population(
    dataset,
    model_spec,
    dvid_mapping={1: "conc_parent", 2: "conc_metabolite"}
)
```

## CLI Usage

```bash
# Read CDISC data and convert to OpenPKPD format
./packages/cli/bin/openpkpd data --format cdisc \
    --pc pc.csv --ex ex.csv --dm dm.csv \
    --out dataset.json
```

## Example Workflow

### From CDISC to Estimation

```python
from openpkpd.data import read_cdisc_csv, cdisc_to_population
from openpkpd.estimation import estimate, EstimationConfig, FOCEIMethod

# 1. Read CDISC data
dataset = read_cdisc_csv("pc.csv", "ex.csv", "dm.csv")

# 2. Validate
warnings = validate_cdisc_dataset(dataset)
print(f"Validation warnings: {len(warnings)}")

# 3. Convert to OpenPKPD format
model_spec = create_model_spec("OneCompOralFirstOrder", initial_params)
pop_spec, observed = cdisc_to_population(dataset, model_spec)

# 4. Run estimation
config = EstimationConfig(
    method=FOCEIMethod(),
    theta_init=[1.0, 5.0, 50.0],
    # ... other settings
)
result = estimate(observed, model_spec, config)
```

### From CDISC to VPC

```python
from openpkpd.data import read_cdisc_csv, cdisc_to_population
from openpkpd.analysis import compute_vpc

# 1. Read and convert data
dataset = read_cdisc_csv("pc.csv", "ex.csv", "dm.csv")
pop_spec, observed = cdisc_to_population(dataset, model_spec)

# 2. Compute VPC
vpc_result = compute_vpc(observed, pop_spec, grid, solver)
```

## See Also

- [Model Import](import.md) - Importing NONMEM/Monolix models
- [Parameter Estimation](estimation.md) - Fitting models to data
- [VPC](vpc.md) - Visual Predictive Checks
