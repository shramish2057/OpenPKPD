# Contributing

## Principles

- Deterministic execution is mandatory.
- Numerical behavior must be explicit and versioned.
- No silent changes to solver semantics.

## Pull requests

A PR must include:

- Tests covering new behavior
- Documentation for any new model or solver option
- Validation impact notes if numerical outputs change

## Style

- Julia formatting uses JuliaFormatter.
- Avoid hidden defaults. Prefer explicit configuration.
