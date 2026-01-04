# Parameter Estimation Module for NLME Models
# Entry point for all estimation functionality

# Types and configuration
include("estimation_types.jl")

# Core estimation algorithm
include("estimate.jl")

# Individual method implementations
include("laplacian.jl")
include("foce.jl")
include("saem.jl")

# Gradients and optimization utilities
include("gradients.jl")

# Standard errors and covariance matrix
include("standard_errors.jl")

# Diagnostics (CWRES, IWRES, etc.)
include("diagnostics.jl")
