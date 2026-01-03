#!/usr/bin/env bash
set -euo pipefail

python -m pip install -e python
python -m pip install pytest

pytest -q python/tests

python python/examples/write_pk_iv_bolus.py
