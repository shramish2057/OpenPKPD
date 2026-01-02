# Design

## Non negotiable invariants

- Deterministic execution
- Explicit solver configuration
- Model specification separated from numerical execution
- Serializable model definitions
- Versioned numerical semantics
- No hidden global state

## Separation of concerns

- Model specification: states, parameters, equations, observation mappings
- Solver layer: integration method, tolerances, event handling
- Engine: validation, binding, execution, metadata emission
