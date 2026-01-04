# NONMEM to OpenPKPD Converter
# Converts parsed NONMEM control files to OpenPKPD model specifications

export convert_nonmem_to_openpkpd, NONMEMConversionResult

"""
Result of converting a NONMEM control file to OpenPKPD format.

Fields:
- model_spec: The converted ModelSpec (or nothing if conversion failed)
- iiv_spec: Inter-individual variability specification
- error_spec: Residual error specification
- warnings: Any warnings generated during conversion
- errors: Any errors that prevented conversion
- parameter_mapping: Mapping from NONMEM THETAs to OpenPKPD parameters
"""
struct NONMEMConversionResult
    model_spec::Union{Nothing,ModelSpec}
    iiv_spec::Union{Nothing,IIVSpec}
    error_spec::Union{Nothing,ResidualErrorSpec}
    warnings::Vector{String}
    errors::Vector{String}
    parameter_mapping::Dict{Int,Symbol}
end

"""
Convert a NONMEM control file to OpenPKPD format.

Arguments:
- ctl: Parsed NONMEMControlFile
- doses: Vector of DoseEvent (required, as NONMEM data is external)
- name: Model name (optional, defaults to problem description)

Returns:
- NONMEMConversionResult containing converted specs and diagnostics
"""
function convert_nonmem_to_openpkpd(
    ctl::NONMEMControlFile;
    doses::Vector{DoseEvent}=DoseEvent[],
    name::String=""
)::NONMEMConversionResult
    warnings = String[]
    errors = String[]
    parameter_mapping = Dict{Int,Symbol}()

    # Use problem description as name if not provided
    if isempty(name)
        name = isempty(ctl.problem) ? "imported_model" : ctl.problem
    end

    # Check for subroutines
    if ctl.subroutines === nothing
        push!(errors, "No \$SUBROUTINES found - cannot determine model type")
        return NONMEMConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping)
    end

    # Get model mapping
    mapping = get_model_mapping(ctl.subroutines.advan, ctl.subroutines.trans)
    if mapping === nothing
        push!(errors, "Unsupported ADVAN$(ctl.subroutines.advan) TRANS$(ctl.subroutines.trans) combination")
        return NONMEMConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping)
    end

    model_kind_sym, param_symbols = mapping

    # Check we have enough THETAs
    n_params = length(param_symbols)
    if length(ctl.thetas) < n_params
        push!(errors, "Expected at least $n_params THETAs for $(model_kind_sym), got $(length(ctl.thetas))")
        return NONMEMConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping)
    end

    # Extract parameter values
    theta_values = [ctl.thetas[i].init for i in 1:n_params]

    # Build parameter mapping
    for (i, sym) in enumerate(param_symbols)
        parameter_mapping[i] = sym
    end

    # Create model spec based on kind
    model_spec = _create_model_spec(model_kind_sym, theta_values, param_symbols, name, doses, warnings)

    if model_spec === nothing
        push!(errors, "Failed to create model spec for $(model_kind_sym)")
        return NONMEMConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping)
    end

    # Convert OMEGA to IIV spec
    iiv_spec = _convert_omega_to_iiv(ctl.omegas, param_symbols, warnings)

    # Convert SIGMA to residual error spec
    error_spec = _convert_sigma_to_error(ctl.sigmas, warnings)

    # Add warning if there are extra THETAs
    if length(ctl.thetas) > n_params
        push!(warnings, "Extra THETAs detected ($(length(ctl.thetas)) vs $n_params expected) - may be covariate effects")
    end

    return NONMEMConversionResult(model_spec, iiv_spec, error_spec, warnings, errors, parameter_mapping)
end

"""
Create a ModelSpec from model kind and parameters.
"""
function _create_model_spec(
    kind_sym::Symbol,
    values::Vector{Float64},
    param_symbols::Vector{Symbol},
    name::String,
    doses::Vector{DoseEvent},
    warnings::Vector{String}
)::Union{Nothing,ModelSpec}

    # Map NONMEM parameter names to OpenPKPD
    param_map = Dict{Symbol,Float64}()
    for (i, sym) in enumerate(param_symbols)
        param_map[sym] = values[i]
    end

    if kind_sym == :OneCompIVBolus
        CL = get(param_map, :CL, get(param_map, :K, 0.0) * get(param_map, :V, 1.0))
        V = get(param_map, :V, 1.0)
        params = OneCompIVBolusParams(CL, V)
        return ModelSpec(OneCompIVBolus(), name, params, doses)

    elseif kind_sym == :OneCompOralFirstOrder
        Ka = get(param_map, :KA, get(param_map, :Ka, 1.0))
        CL = get(param_map, :CL, get(param_map, :K, 0.0) * get(param_map, :V, 1.0))
        V = get(param_map, :V, 1.0)
        params = OneCompOralFirstOrderParams(Ka, CL, V)
        return ModelSpec(OneCompOralFirstOrder(), name, params, doses)

    elseif kind_sym == :TwoCompIVBolus
        CL = get(param_map, :CL, 1.0)
        V1 = get(param_map, :V1, get(param_map, :V, 1.0))
        Q = get(param_map, :Q, 1.0)
        V2 = get(param_map, :V2, get(param_map, :VSS, V1) - V1)
        params = TwoCompIVBolusParams(CL, V1, Q, V2)
        return ModelSpec(TwoCompIVBolus(), name, params, doses)

    elseif kind_sym == :TwoCompOral
        Ka = get(param_map, :KA, get(param_map, :Ka, 1.0))
        CL = get(param_map, :CL, 1.0)
        V1 = get(param_map, :V1, get(param_map, :V, 1.0))
        Q = get(param_map, :Q, 1.0)
        V2 = get(param_map, :V2, get(param_map, :VSS, V1) - V1)
        params = TwoCompOralParams(Ka, CL, V1, Q, V2)
        return ModelSpec(TwoCompOral(), name, params, doses)

    elseif kind_sym == :ThreeCompIVBolus
        CL = get(param_map, :CL, 1.0)
        V1 = get(param_map, :V1, 1.0)
        Q2 = get(param_map, :Q2, 1.0)
        V2 = get(param_map, :V2, 1.0)
        Q3 = get(param_map, :Q3, 1.0)
        V3 = get(param_map, :V3, 1.0)
        params = ThreeCompIVBolusParams(CL, V1, Q2, V2, Q3, V3)
        return ModelSpec(ThreeCompIVBolus(), name, params, doses)

    elseif kind_sym == :MichaelisMentenElimination
        Vmax = get(param_map, :VM, get(param_map, :Vmax, 1.0))
        Km = get(param_map, :KM, get(param_map, :Km, 1.0))
        V = get(param_map, :V, 1.0)
        params = MichaelisMentenEliminationParams(Vmax, Km, V)
        return ModelSpec(MichaelisMentenElimination(), name, params, doses)
    end

    push!(warnings, "Model kind $(kind_sym) not yet supported for conversion")
    return nothing
end

"""
Convert OMEGA blocks to IIV specification.
"""
function _convert_omega_to_iiv(
    omegas::Vector{OMEGABlock},
    param_symbols::Vector{Symbol},
    warnings::Vector{String}
)::Union{Nothing,IIVSpec}
    if isempty(omegas)
        return nothing
    end

    # Collect all omega values
    all_omegas = Float64[]
    for block in omegas
        if block.structure == :diagonal
            append!(all_omegas, block.values)
        else
            # For block structure, extract diagonals (simplified)
            append!(all_omegas, block.values)
        end
    end

    if isempty(all_omegas)
        return nothing
    end

    # Map omegas to parameters
    omegas_dict = Dict{Symbol,Float64}()
    for (i, omega) in enumerate(all_omegas)
        if i <= length(param_symbols)
            # Convert variance to SD (omega in NONMEM is variance)
            omegas_dict[param_symbols[i]] = sqrt(omega)
        end
    end

    if isempty(omegas_dict)
        return nothing
    end

    return IIVSpec(LogNormalIIV(), omegas_dict, UInt64(12345), 1)
end

"""
Convert SIGMA to residual error specification.
"""
function _convert_sigma_to_error(
    sigmas::Vector{SIGMABlock},
    warnings::Vector{String}
)::Union{Nothing,ResidualErrorSpec}
    if isempty(sigmas)
        return nothing
    end

    # Get first sigma value
    sigma_values = sigmas[1].values
    if isempty(sigma_values)
        return nothing
    end

    if length(sigma_values) == 1
        # Single sigma - assume proportional error
        # NONMEM stores variance, convert to SD
        sigma = sqrt(sigma_values[1])
        return ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(sigma),
            :conc,
            UInt64(12345)
        )
    elseif length(sigma_values) >= 2
        # Two sigmas - assume combined error (additive + proportional)
        sigma_add = sqrt(sigma_values[1])
        sigma_prop = sqrt(sigma_values[2])
        return ResidualErrorSpec(
            CombinedError(),
            CombinedErrorParams(sigma_add, sigma_prop),
            :conc,
            UInt64(12345)
        )
    end

    return nothing
end
