# NONMEM Control File Types
# Structures for representing parsed NONMEM .ctl files

export NONMEMControlFile, THETASpec, OMEGABlock, SIGMABlock
export SubroutineSpec, DataSpec, InputColumn

"""
Specification for a single THETA parameter.

Fields:
- lower: Lower bound (can be -Inf)
- init: Initial estimate
- upper: Upper bound (can be Inf)
- fixed: Whether parameter is fixed
- name: Optional parameter name/label
"""
struct THETASpec
    lower::Float64
    init::Float64
    upper::Float64
    fixed::Bool
    name::String

    function THETASpec(lower::Float64, init::Float64, upper::Float64, fixed::Bool=false, name::String="")
        if !fixed && lower > init
            error("THETA lower bound ($lower) > init ($init)")
        end
        if !fixed && init > upper
            error("THETA init ($init) > upper bound ($upper)")
        end
        new(lower, init, upper, fixed, name)
    end
end

# Convenience constructor for simple (INIT) format
THETASpec(init::Float64) = THETASpec(-Inf, init, Inf, false, "")

"""
OMEGA block specification.

Represents variance-covariance matrix for inter-individual variability.

Fields:
- values: Variance/covariance values (diagonal or block)
- structure: :diagonal, :block, or :same
- dimension: Size of the block
- fixed: Whether block is fixed
"""
struct OMEGABlock
    values::Vector{Float64}
    structure::Symbol
    dimension::Int
    fixed::Bool

    function OMEGABlock(values::Vector{Float64}, structure::Symbol=:diagonal, fixed::Bool=false)
        dim = if structure == :diagonal
            length(values)
        else
            # For block, n*(n+1)/2 = length, solve for n
            n = Int((-1 + sqrt(1 + 8*length(values))) / 2)
            n
        end
        new(values, structure, dim, fixed)
    end
end

"""
SIGMA block specification.

Represents residual error variance.

Fields:
- values: Variance values
- structure: :diagonal or :block
- fixed: Whether fixed
"""
struct SIGMABlock
    values::Vector{Float64}
    structure::Symbol
    fixed::Bool

    function SIGMABlock(values::Vector{Float64}, structure::Symbol=:diagonal, fixed::Bool=false)
        new(values, structure, fixed)
    end
end

"""
\$SUBROUTINES specification (ADVAN/TRANS).

Fields:
- advan: ADVAN number (1-13)
- trans: TRANS number (1-6)
- other: Other subroutines specified
"""
struct SubroutineSpec
    advan::Int
    trans::Int
    other::Vector{String}

    function SubroutineSpec(advan::Int, trans::Int=1, other::Vector{String}=String[])
        advan in 1:13 || error("ADVAN must be 1-13, got $advan")
        trans in 1:6 || error("TRANS must be 1-6, got $trans")
        new(advan, trans, other)
    end
end

"""
\$DATA specification.

Fields:
- filename: Path to data file
- ignore: Characters/conditions to ignore
- accept: Conditions to accept
"""
struct DataSpec
    filename::String
    ignore::Vector{String}
    accept::Vector{String}

    DataSpec(filename::String) = new(filename, String[], String[])
    DataSpec(filename::String, ignore::Vector{String}, accept::Vector{String}) = new(filename, ignore, accept)
end

"""
\$INPUT column specification.

Fields:
- name: Column name
- drop: Whether to drop this column
- alias: Alias for standard names (ID, TIME, DV, etc.)
"""
struct InputColumn
    name::String
    drop::Bool
    alias::String

    InputColumn(name::String) = new(name, false, "")
    InputColumn(name::String, drop::Bool) = new(name, drop, "")
    InputColumn(name::String, drop::Bool, alias::String) = new(name, drop, alias)
end

"""
Complete NONMEM control file representation.

This structure captures the essential elements needed to convert
a NONMEM model to OpenPKPD format.

Fields:
- problem: Problem description from \$PROBLEM
- data: Data file specification
- input: Input column definitions
- subroutines: ADVAN/TRANS specification
- thetas: Vector of THETA specifications
- omegas: Vector of OMEGA blocks
- sigmas: Vector of SIGMA blocks
- pk_code: Lines from \$PK block
- error_code: Lines from \$ERROR block
- estimation: Estimation method settings
- tables: Table output specifications
- raw_text: Original control file text
"""
struct NONMEMControlFile
    problem::String
    data::Union{Nothing,DataSpec}
    input::Vector{InputColumn}
    subroutines::Union{Nothing,SubroutineSpec}
    thetas::Vector{THETASpec}
    omegas::Vector{OMEGABlock}
    sigmas::Vector{SIGMABlock}
    pk_code::Vector{String}
    error_code::Vector{String}
    estimation::Dict{String,Any}
    tables::Vector{Dict{String,Any}}
    raw_text::String
end

# ADVAN/TRANS to Model Kind mapping
const ADVAN_TRANS_MAP = Dict{Tuple{Int,Int},Tuple{Symbol,Vector{Symbol}}}(
    # ADVAN1: One-compartment IV bolus
    (1, 1) => (:OneCompIVBolus, [:K, :V]),      # K parameterization
    (1, 2) => (:OneCompIVBolus, [:CL, :V]),     # CL/V parameterization

    # ADVAN2: One-compartment oral first-order
    (2, 1) => (:OneCompOralFirstOrder, [:KA, :K, :V]),
    (2, 2) => (:OneCompOralFirstOrder, [:KA, :CL, :V]),

    # ADVAN3: Two-compartment IV bolus
    (3, 1) => (:TwoCompIVBolus, [:K, :K12, :K21, :V]),     # Micro-constants
    (3, 3) => (:TwoCompIVBolus, [:CL, :V, :Q, :VSS]),      # CL, V, Q, Vss
    (3, 4) => (:TwoCompIVBolus, [:CL, :V1, :Q, :V2]),      # CL, V1, Q, V2

    # ADVAN4: Two-compartment oral first-order
    (4, 1) => (:TwoCompOral, [:KA, :K, :K23, :K32, :V]),
    (4, 3) => (:TwoCompOral, [:KA, :CL, :V, :Q, :VSS]),
    (4, 4) => (:TwoCompOral, [:KA, :CL, :V1, :Q, :V2]),

    # ADVAN11: Three-compartment IV bolus
    (11, 1) => (:ThreeCompIVBolus, [:K, :K12, :K21, :K13, :K31, :V]),
    (11, 4) => (:ThreeCompIVBolus, [:CL, :V1, :Q2, :V2, :Q3, :V3]),

    # ADVAN10: Michaelis-Menten elimination
    (10, 1) => (:MichaelisMentenElimination, [:VM, :KM, :V]),
)

"""
Get the OpenPKPD model kind and expected parameters for an ADVAN/TRANS combination.

Returns (model_kind_symbol, parameter_symbols) or nothing if not supported.
"""
function get_model_mapping(advan::Int, trans::Int)
    key = (advan, trans)
    if haskey(ADVAN_TRANS_MAP, key)
        return ADVAN_TRANS_MAP[key]
    end
    return nothing
end

export ADVAN_TRANS_MAP, get_model_mapping
