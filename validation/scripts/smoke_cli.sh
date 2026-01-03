#!/usr/bin/env bash
set -euo pipefail

# Instantiate CLI dependencies first
julia -e 'using Pkg; Pkg.activate("cli/OpenPKPDCLI"); Pkg.instantiate()'

./bin/openpkpd version
./bin/openpkpd replay --artifact validation/golden/pk_iv_bolus.json
./bin/openpkpd validate-golden
