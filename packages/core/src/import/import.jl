# Model Import Module
# Import NONMEM and Monolix models to OpenPKPD format
#
# This module provides parsers for external pharmacometrics tools:
# - NONMEM: Parse .ctl control files
# - Monolix: Parse .mlxtran project files

# NONMEM import
include("nonmem_types.jl")
include("nonmem_parser.jl")
include("nonmem_converter.jl")

# Monolix import (coming soon)
include("monolix_types.jl")
include("monolix_parser.jl")
include("monolix_converter.jl")
