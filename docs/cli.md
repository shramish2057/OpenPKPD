# CLI Reference

OpenPKPD provides a command-line interface for common operations.

## Installation

The CLI is located in `bin/openpkpd` and requires Julia with the OpenPKPDCore and OpenPKPDCLI packages.

```bash
# Make executable
chmod +x bin/openpkpd

# Run CLI
./bin/openpkpd <command> [options]
```

---

## Commands

### version

Display version information for all OpenPKPD components.

```bash
./bin/openpkpd version
```

**Output**:

```
OpenPKPD 0.1.0
Event semantics: 1.0.0
Solver semantics: 1.0.0
Artifact schema: 1.0.0
```

---

### replay

Replay an execution artifact to verify reproducibility.

```bash
./bin/openpkpd replay --artifact <path> [--out <output_path>]
```

**Options**:

| Option | Required | Description |
|--------|----------|-------------|
| `--artifact` | Yes | Path to artifact JSON file |
| `--out` | No | Path to write replayed artifact |

**Supported Artifact Types**:

- Single execution (`artifact_type: "single"` or missing)
- Population execution (`artifact_type: "population"`)
- Single sensitivity (`artifact_type: "sensitivity_single"`)
- Population sensitivity (`artifact_type: "sensitivity_population"`)

**Examples**:

```bash
# Replay and verify
./bin/openpkpd replay --artifact validation/golden/pk_iv_bolus.json

# Replay and save output
./bin/openpkpd replay --artifact my_simulation.json --out replayed.json

# Replay population artifact
./bin/openpkpd replay --artifact validation/golden/population_iv_bolus.json
```

---

### validate-golden

Run the full golden artifact validation suite.

```bash
./bin/openpkpd validate-golden
```

**What It Does**:

1. Finds all golden artifacts in `validation/golden/`
2. Replays each artifact
3. Compares replayed results to stored results
4. Reports any discrepancies

**Output**:

```
Validating golden artifacts...
  pk_iv_bolus.json: PASS
  pk_oral_first_order.json: PASS
  population_iv_bolus.json: PASS
  pkpd_direct_emax.json: PASS
  sensitivity_single.json: PASS
  ...
All 10 golden artifacts validated successfully.
```

**Exit Codes**:

| Code | Meaning |
|------|---------|
| 0 | All validations passed |
| 1 | One or more validations failed |

---

## Artifact Format

All OpenPKPD artifacts are JSON files with a consistent structure.

### Common Fields

```json
{
  "artifact_schema_version": "1.0.0",
  "semantics_fingerprint": {
    "artifact_schema_version": "1.0.0",
    "event_semantics_version": "1.0.0",
    "solver_semantics_version": "1.0.0"
  }
}
```

### Single Execution Artifact

```json
{
  "artifact_schema_version": "1.0.0",
  "execution_mode": "pk",
  "model_spec": {
    "kind": "OneCompIVBolus",
    "name": "example",
    "params": {"CL": 5.0, "V": 50.0},
    "doses": [{"time": 0.0, "amount": 100.0}]
  },
  "grid": {
    "t0": 0.0,
    "t1": 24.0,
    "saveat": [0.0, 1.0, 2.0, ...]
  },
  "solver": {
    "alg": "Tsit5",
    "reltol": 1e-10,
    "abstol": 1e-12,
    "maxiters": 10000000
  },
  "result": {
    "t": [0.0, 1.0, 2.0, ...],
    "states": {"A_central": [...]},
    "observations": {"conc": [...]},
    "metadata": {...}
  }
}
```

### Population Artifact

```json
{
  "artifact_type": "population",
  "population_spec": {
    "base_model_spec": {...},
    "iiv": {
      "kind": "LogNormalIIV",
      "omegas": {"CL": 0.3, "V": 0.2},
      "seed": 12345,
      "n": 100
    },
    "iov": null,
    "covariate_model": null,
    "covariates": []
  },
  "result": {
    "individuals": [...],
    "params": [...],
    "summaries": {
      "conc": {
        "observation": "conc",
        "probs": [0.05, 0.95],
        "mean": [...],
        "median": [...],
        "quantiles": {...}
      }
    },
    "metadata": {...}
  }
}
```

---

## Integration Examples

### CI/CD Pipeline

```yaml
# .github/workflows/validate.yml
name: Golden Validation
on: [push, pull_request]
jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: julia-actions/setup-julia@v1
      - run: julia --project=core/OpenPKPDCore -e 'using Pkg; Pkg.instantiate()'
      - run: ./bin/openpkpd validate-golden
```

### Batch Replay

```bash
#!/bin/bash
# replay_all.sh - Replay all artifacts in a directory

for artifact in artifacts/*.json; do
    echo "Replaying: $artifact"
    ./bin/openpkpd replay --artifact "$artifact" --out "replayed/$(basename $artifact)"
done
```

### Compare Artifacts

```bash
# Generate new artifact
julia --project=core/OpenPKPDCore my_simulation.jl

# Compare with golden
diff <(jq -S . validation/golden/pk_iv_bolus.json) <(jq -S . my_artifact.json)
```

---

## Troubleshooting

### Julia Not Found

```
Error: julia: command not found
```

**Solution**: Ensure Julia is installed and in your PATH.

### Package Not Installed

```
Error: ArgumentError: Package OpenPKPDCore not found
```

**Solution**: Install dependencies:

```bash
julia --project=core/OpenPKPDCore -e 'using Pkg; Pkg.instantiate()'
julia --project=cli/OpenPKPDCLI -e 'using Pkg; Pkg.instantiate()'
```

### Artifact Schema Mismatch

```
Warning: Artifact schema version mismatch
```

**Cause**: Artifact was created with a different version of OpenPKPD.

**Solution**: Re-generate the artifact or update your OpenPKPD installation.

---

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `JULIA_PROJECT` | Julia project path | Auto-detected |
| `OPENPKPD_ROOT` | Repository root | Auto-detected |

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | General error |
| 2 | Invalid arguments |
| 3 | File not found |
| 4 | Validation failure |
