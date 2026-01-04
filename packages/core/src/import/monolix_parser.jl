# Monolix Project File Parser
# Parses .mlxtran files into MonolixProject structures

export parse_monolix_project, read_monolix_project

"""
Read and parse a Monolix project file from disk.

Arguments:
- filepath: Path to .mlxtran file

Returns:
- MonolixProject structure
"""
function read_monolix_project(filepath::AbstractString)::MonolixProject
    text = read(filepath, String)
    return parse_monolix_project(text)
end

"""
Parse Monolix project file text into a structured representation.

Arguments:
- text: Raw .mlxtran file text

Returns:
- MonolixProject structure
"""
function parse_monolix_project(text::String)::MonolixProject
    # Extract blocks
    blocks = _extract_monolix_blocks(text)

    # Parse each section
    description = _parse_monolix_description(get(blocks, "DESCRIPTION", ""))
    data = _parse_monolix_data(get(blocks, "DATAFILE", ""), get(blocks, "DATARELATION", ""))
    model = _parse_monolix_model(get(blocks, "STRUCTURAL_MODEL", ""))
    parameters = _parse_monolix_parameters(get(blocks, "PARAMETER", ""), get(blocks, "POPULATION", ""))
    observations = _parse_monolix_observations(get(blocks, "OBSERVATION_MODEL", ""))
    estimation = _parse_monolix_estimation(get(blocks, "FIT", ""))

    return MonolixProject(
        description,
        data,
        model,
        parameters,
        observations,
        estimation,
        text
    )
end

"""
Extract blocks from Monolix project file.
Returns Dict mapping block name to block content.
"""
function _extract_monolix_blocks(text::String)::Dict{String,String}
    blocks = Dict{String,String}()
    current_block = nothing
    current_content = String[]

    for line in split(text, '\n')
        stripped = strip(line)

        # Check for block start
        m = match(r"^\[([A-Z_]+)\]", stripped)
        if m !== nothing
            # Save previous block
            if current_block !== nothing
                blocks[current_block] = join(current_content, "\n")
            end
            current_block = m.captures[1]
            current_content = String[]
        elseif current_block !== nothing
            push!(current_content, stripped)
        end
    end

    # Save last block
    if current_block !== nothing
        blocks[current_block] = join(current_content, "\n")
    end

    return blocks
end

"""
Parse [DESCRIPTION] block.
"""
function _parse_monolix_description(text::String)::String
    return strip(text)
end

"""
Parse [DATAFILE] and [DATARELATION] blocks.
"""
function _parse_monolix_data(datafile_text::String, relation_text::String)::Union{Nothing,MonolixDataset}
    if isempty(datafile_text)
        return nothing
    end

    filename = ""
    header_types = Dict{String,String}()

    for line in split(datafile_text, '\n')
        line = strip(line)
        if isempty(line)
            continue
        end

        # Parse file = 'path'
        m = match(r"file\s*=\s*'([^']+)'", line)
        if m !== nothing
            filename = String(m.captures[1])
            continue
        end

        # Parse header = {...}
        m = match(r"header\s*=\s*\{([^}]+)\}", line)
        if m !== nothing
            # Parse column types
            for part in split(m.captures[1], ',')
                part = strip(part)
                if contains(part, '=')
                    name, type = split(part, '=', limit=2)
                    header_types[String(strip(name))] = String(strip(type))
                end
            end
        end
    end

    # Determine standard columns from header_types
    id_col = "ID"
    time_col = "TIME"
    obs_col = "DV"
    dose_col = "AMT"
    rate_col = "RATE"

    for (name, type) in header_types
        type_upper = uppercase(type)
        if type_upper == "ID"
            id_col = name
        elseif type_upper == "TIME"
            time_col = name
        elseif type_upper in ["OBSERVATION", "DV"]
            obs_col = name
        elseif type_upper in ["AMOUNT", "AMT", "DOSE"]
            dose_col = name
        elseif type_upper == "RATE"
            rate_col = name
        end
    end

    return MonolixDataset(
        filename;
        header_types=header_types,
        id_column=id_col,
        time_column=time_col,
        observation_column=obs_col,
        dose_column=dose_col,
        rate_column=rate_col
    )
end

"""
Parse [STRUCTURAL_MODEL] block.
"""
function _parse_monolix_model(text::String)::Union{Nothing,MonolixStructuralModel}
    if isempty(text)
        return nothing
    end

    model_type = nothing
    admin_type = "iv"
    n_compartments = 1
    elimination = "linear"
    absorption = "bolus"
    has_lag = false
    has_bio = false

    for line in split(text, '\n')
        line = strip(line)
        if isempty(line)
            continue
        end

        # Parse lib = pklib:model_name
        m = match(r"lib\s*=\s*(\w+):(\S+)", line)
        if m !== nothing
            model_type = MonolixModelType(String(m.captures[1]), String(m.captures[2]))

            # Infer properties from model name
            model_name = lowercase(String(m.captures[2]))
            if occursin("oral", model_name)
                admin_type = "oral"
                absorption = "firstOrder"
            end
            if occursin("1cpt", model_name)
                n_compartments = 1
            elseif occursin("2cpt", model_name)
                n_compartments = 2
            elseif occursin("3cpt", model_name)
                n_compartments = 3
            end
            if occursin("tlag", model_name)
                has_lag = true
            end
            if occursin("tbio", model_name) || occursin("bio", model_name)
                has_bio = true
            end
            if occursin("mm", model_name) || occursin("michaelis", model_name)
                elimination = "mm"
            end
            continue
        end

        # Parse administration = ...
        m = match(r"administration\s*=\s*(\w+)", line)
        if m !== nothing
            admin_type = lowercase(String(m.captures[1]))
            continue
        end

        # Parse absorption = ...
        m = match(r"absorption\s*=\s*(\w+)", line)
        if m !== nothing
            absorption = lowercase(String(m.captures[1]))
            continue
        end

        # Parse elimination = ...
        m = match(r"elimination\s*=\s*(\w+)", line)
        if m !== nothing
            elimination = lowercase(String(m.captures[1]))
            continue
        end
    end

    return MonolixStructuralModel(
        model_type,
        admin_type,
        n_compartments,
        elimination,
        absorption,
        has_lag,
        has_bio
    )
end

"""
Parse [PARAMETER] and [POPULATION] blocks.
"""
function _parse_monolix_parameters(param_text::String, pop_text::String)::Vector{MonolixParameter}
    parameters = MonolixParameter[]

    # Parse population parameters first for IIV info
    omega_values = Dict{String,Float64}()
    distributions = Dict{String,String}()

    for line in split(pop_text, '\n')
        line = strip(line)
        if isempty(line)
            continue
        end

        # Parse parameter definitions with omega
        # e.g., ka = {value=1.0, method=MLE}
        m = match(r"(\w+)\s*=\s*\{([^}]+)\}", line)
        if m !== nothing
            name = String(m.captures[1])
            props = String(m.captures[2])

            # Extract omega_xxx
            omega_match = match(r"omega\s*=\s*([0-9.eE+-]+)", props)
            if omega_match !== nothing
                omega_values[name] = parse(Float64, omega_match.captures[1])
            end

            # Extract distribution
            dist_match = match(r"distribution\s*=\s*(\w+)", props)
            if dist_match !== nothing
                distributions[name] = String(dist_match.captures[1])
            end
        end
    end

    # Parse individual parameters
    for line in split(param_text, '\n')
        line = strip(line)
        if isempty(line)
            continue
        end

        # Parse parameter = value or parameter = {props}
        m = match(r"(\w+)\s*=\s*([0-9.eE+-]+)", line)
        if m !== nothing
            name = String(m.captures[1])
            value = parse(Float64, m.captures[2])
            omega = get(omega_values, name, 0.0)
            dist = get(distributions, name, "logNormal")
            has_iiv = omega > 0.0

            push!(parameters, MonolixParameter(
                name, value;
                distribution=dist,
                omega=omega,
                has_iiv=has_iiv
            ))
            continue
        end

        # Parse parameter = {value=..., method=...}
        m = match(r"(\w+)\s*=\s*\{([^}]+)\}", line)
        if m !== nothing
            name = String(m.captures[1])
            props = String(m.captures[2])

            # Extract value
            val_match = match(r"value\s*=\s*([0-9.eE+-]+)", props)
            value = val_match !== nothing ? parse(Float64, val_match.captures[1]) : 1.0

            # Check if fixed
            fixed = occursin("method=FIXED", uppercase(props))

            omega = get(omega_values, name, 0.0)
            dist = get(distributions, name, "logNormal")
            has_iiv = omega > 0.0

            push!(parameters, MonolixParameter(
                name, value;
                fixed=fixed,
                distribution=dist,
                omega=omega,
                has_iiv=has_iiv
            ))
        end
    end

    return parameters
end

"""
Parse [OBSERVATION_MODEL] block.
"""
function _parse_monolix_observations(text::String)::Vector{MonolixObservation}
    observations = MonolixObservation[]

    if isempty(text)
        return observations
    end

    for line in split(text, '\n')
        line = strip(line)
        if isempty(line)
            continue
        end

        # Parse observation definition
        # e.g., y1 = {type=continuous, prediction=Cc, error=combined1}
        m = match(r"(\w+)\s*=\s*\{([^}]+)\}", line)
        if m !== nothing
            name = String(m.captures[1])
            props = String(m.captures[2])

            obs_type = "continuous"
            error_model = "combined"
            error_params = Float64[]

            # Extract type
            type_match = match(r"type\s*=\s*(\w+)", props)
            if type_match !== nothing
                obs_type = lowercase(String(type_match.captures[1]))
            end

            # Extract error model
            error_match = match(r"error\s*=\s*(\w+)", props)
            if error_match !== nothing
                error_str = lowercase(String(error_match.captures[1]))
                if occursin("constant", error_str) || occursin("additive", error_str)
                    error_model = "additive"
                elseif occursin("proportional", error_str)
                    error_model = "proportional"
                elseif occursin("combined", error_str)
                    error_model = "combined"
                elseif occursin("exponential", error_str)
                    error_model = "exponential"
                end
            end

            push!(observations, MonolixObservation(
                name;
                type=obs_type,
                error_model=error_model,
                error_params=error_params
            ))
        end
    end

    return observations
end

"""
Parse [FIT] or estimation method blocks.
"""
function _parse_monolix_estimation(text::String)::String
    if isempty(text)
        return "SAEM"  # Default Monolix estimation method
    end

    text_upper = uppercase(text)

    if occursin("SAEM", text_upper)
        return "SAEM"
    elseif occursin("FOCE", text_upper)
        return "FOCE"
    elseif occursin("LAPLACE", text_upper)
        return "LAPLACIAN"
    end

    return "SAEM"
end
