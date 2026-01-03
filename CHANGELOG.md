# Changelog

All notable changes to OpenPKPD will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Comprehensive documentation covering all features
- Python bindings with full simulation support
- CLI for artifact replay and validation
- Time-varying covariates (step and linear interpolation)
- Inter-occasion variability (IOV) support
- Coupled PK-PD simulation with indirect response models
- Sensitivity analysis (single and population)
- Golden artifact validation system

### Changed
- Repository restructured for professional organization
- Documentation moved to MkDocs with Material theme

## [0.1.0] - 2024-01-03

### Added
- Initial release
- One-compartment IV bolus PK model
- One-compartment oral first-order absorption PK model
- Direct Emax PD model
- Indirect response turnover PD model
- Population simulation with IIV
- Covariate effects (linear, power, exponential)
- JSON artifact serialization with semantic versioning
- Deterministic replay system
- Golden artifact validation

### Semantics Versions
- Event Semantics: 1.0.0
- Solver Semantics: 1.0.0
- Artifact Schema: 1.0.0

[Unreleased]: https://github.com/openpkpd/openpkpd/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/openpkpd/openpkpd/releases/tag/v0.1.0
