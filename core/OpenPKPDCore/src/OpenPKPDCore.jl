module OpenPKPDCore

# ------------------------------------------------------------------
# External dependencies
# ------------------------------------------------------------------
using SciMLBase
using DifferentialEquations

# ------------------------------------------------------------------
# Core specs and shared types 
# ------------------------------------------------------------------
include("specs/specs.jl")
include("specs/sensitivity.jl")

# ------------------------------------------------------------------
# PK model definitions
# ------------------------------------------------------------------
include("models/onecomp_iv_bolus.jl")
include("models/onecomp_oral_first_order.jl")
include("models/pk_interface.jl")

# ------------------------------------------------------------------
# PD model definitions
# ------------------------------------------------------------------
include("pd/direct_emax.jl")
include("pd/indirect_response_turnover.jl")

# ------------------------------------------------------------------
# Numerical semantics 
# These define versioned scientific meaning
# ------------------------------------------------------------------
include("engine/semantics.jl")
include("engine/solver_semantics.jl")
include("engine/semantics_fingerprint.jl")

# ------------------------------------------------------------------
# Perturbation + sensitivity core
# ------------------------------------------------------------------
include("engine/perturb.jl")
include("engine/sensitivity_metrics.jl")
include("engine/sensitivity.jl")
include("engine/sensitivity_population.jl")

# ------------------------------------------------------------------
# Core simulation engine
# ------------------------------------------------------------------
include("engine/events.jl")
include("engine/solve.jl")

# ------------------------------------------------------------------
# PK–PD execution layers
# ------------------------------------------------------------------
include("engine/pkpd.jl")
include("engine/pkpd_coupled.jl")

# ------------------------------------------------------------------
# Population engine
# Defines PopulationSpec, PopulationResult, IIV, etc.
# ------------------------------------------------------------------
include("engine/population.jl")

# ------------------------------------------------------------------
# Serialization
# ------------------------------------------------------------------
include("serialization/schema.jl")
include("serialization/serialize.jl")
include("serialization/deserialize.jl")
include("serialization/serialize_population.jl")
include("serialization/deserialize_population.jl")

end
