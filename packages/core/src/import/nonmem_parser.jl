# NONMEM Control File Parser
# Parses .ctl files into NONMEMControlFile structures

export parse_nonmem_control, read_nonmem_control

"""
Read and parse a NONMEM control file from disk.

Arguments:
- filepath: Path to .ctl file

Returns:
- NONMEMControlFile structure
"""
function read_nonmem_control(filepath::AbstractString)::NONMEMControlFile
    text = read(filepath, String)
    return parse_nonmem_control(text)
end

"""
Parse NONMEM control file text into a structured representation.

Arguments:
- text: Raw control file text

Returns:
- NONMEMControlFile structure
"""
function parse_nonmem_control(text::String)::NONMEMControlFile
    # Remove comments (lines starting with ; or text after ;)
    lines = String[]
    for line in split(text, '\n')
        # Remove inline comments
        comment_idx = findfirst(';', line)
        if comment_idx !== nothing
            line = line[1:comment_idx-1]
        end
        push!(lines, rstrip(line))
    end

    # Find all $SECTION blocks
    sections = _extract_sections(lines)

    # Parse each section
    problem = _parse_problem(get(sections, "PROBLEM", String[]))
    data = _parse_data(get(sections, "DATA", String[]))
    input = _parse_input(get(sections, "INPUT", String[]))
    subroutines = _parse_subroutines(get(sections, "SUBROUTINES", String[]))
    thetas = _parse_theta(get(sections, "THETA", String[]))
    omegas = _parse_omega(get(sections, "OMEGA", String[]))
    sigmas = _parse_sigma(get(sections, "SIGMA", String[]))
    pk_code = get(sections, "PK", String[])
    error_code = get(sections, "ERROR", String[])
    estimation = _parse_estimation(get(sections, "ESTIMATION", String[]))
    tables = _parse_tables(get(sections, "TABLE", String[]))

    return NONMEMControlFile(
        problem,
        data,
        input,
        subroutines,
        thetas,
        omegas,
        sigmas,
        pk_code,
        error_code,
        estimation,
        tables,
        text
    )
end

"""
Extract sections from control file lines.
Returns Dict mapping section name to lines within that section.
"""
function _extract_sections(lines::Vector{String})::Dict{String,Vector{String}}
    sections = Dict{String,Vector{String}}()
    current_section = nothing
    current_lines = String[]

    for line in lines
        stripped = strip(line)
        if isempty(stripped)
            continue
        end

        # Check if this is a section header
        if startswith(stripped, '$')
            # Save previous section
            if current_section !== nothing
                sections[current_section] = current_lines
            end

            # Parse new section name
            parts = split(stripped)
            section_name = uppercase(replace(parts[1], '$' => ""))

            # Handle $SUB as alias for $SUBROUTINES
            if section_name == "SUB"
                section_name = "SUBROUTINES"
            end

            current_section = section_name
            # Include rest of line if present
            if length(parts) > 1
                current_lines = [join(parts[2:end], " ")]
            else
                current_lines = String[]
            end
        elseif current_section !== nothing
            push!(current_lines, stripped)
        end
    end

    # Save last section
    if current_section !== nothing
        sections[current_section] = current_lines
    end

    return sections
end

"""
Parse \$PROBLEM section.
"""
function _parse_problem(lines::Vector{String})::String
    return join(lines, " ")
end

"""
Parse \$DATA section.
"""
function _parse_data(lines::Vector{String})::Union{Nothing,DataSpec}
    if isempty(lines)
        return nothing
    end

    text = join(lines, " ")
    parts = split(text)

    if isempty(parts)
        return nothing
    end

    filename = String(parts[1])
    ignore = String[]
    accept = String[]

    for i in 2:length(parts)
        part = uppercase(String(parts[i]))
        if startswith(part, "IGNORE=")
            push!(ignore, replace(part, "IGNORE=" => ""))
        elseif startswith(part, "ACCEPT=")
            push!(accept, replace(part, "ACCEPT=" => ""))
        end
    end

    return DataSpec(filename, ignore, accept)
end

"""
Parse \$INPUT section.
"""
function _parse_input(lines::Vector{String})::Vector{InputColumn}
    if isempty(lines)
        return InputColumn[]
    end

    text = join(lines, " ")
    # Handle continuation and multiple spaces
    text = replace(text, r"\s+" => " ")

    columns = InputColumn[]
    for part in split(text)
        part = strip(part)
        if isempty(part)
            continue
        end

        # Check for DROP
        if uppercase(part) == "DROP" || uppercase(part) == "SKIP"
            continue
        end

        # Check for NAME=DROP or NAME=ALIAS format
        if contains(part, "=")
            name, value = split(part, "=", limit=2)
            name = String(strip(name))
            value = String(strip(uppercase(value)))

            if value == "DROP" || value == "SKIP"
                push!(columns, InputColumn(name, true))
            else
                push!(columns, InputColumn(name, false, value))
            end
        else
            push!(columns, InputColumn(String(part), false))
        end
    end

    return columns
end

"""
Parse \$SUBROUTINES section.
"""
function _parse_subroutines(lines::Vector{String})::Union{Nothing,SubroutineSpec}
    if isempty(lines)
        return nothing
    end

    text = uppercase(join(lines, " "))
    text = replace(text, r"\s+" => " ")

    advan = 0
    trans = 1
    other = String[]

    for part in split(text)
        part = strip(part)
        if startswith(part, "ADVAN")
            advan_str = replace(part, "ADVAN" => "")
            advan = parse(Int, advan_str)
        elseif startswith(part, "TRANS")
            trans_str = replace(part, "TRANS" => "")
            trans = parse(Int, trans_str)
        elseif !isempty(part)
            push!(other, part)
        end
    end

    if advan == 0
        return nothing
    end

    return SubroutineSpec(advan, trans, other)
end

"""
Parse \$THETA section.
"""
function _parse_theta(lines::Vector{String})::Vector{THETASpec}
    if isempty(lines)
        return THETASpec[]
    end

    text = join(lines, " ")
    thetas = THETASpec[]

    # Handle parenthesized specifications
    # Formats: (LOWER, INIT, UPPER), (INIT), INIT, (INIT) FIX
    i = 1
    while i <= length(text)
        # Skip whitespace
        while i <= length(text) && isspace(text[i])
            i += 1
        end
        if i > length(text)
            break
        end

        if text[i] == '('
            # Find matching closing paren
            j = findnext(')', text, i)
            if j === nothing
                break
            end

            inner = text[i+1:j-1]
            parts = [strip(p) for p in split(inner, ',')]

            # Check for FIX after closing paren
            rest = j < length(text) ? text[j+1:end] : ""
            fixed = occursin(r"^\s*FIX"i, rest)

            if length(parts) == 1
                init = parse(Float64, parts[1])
                push!(thetas, THETASpec(-Inf, init, Inf, fixed))
            elseif length(parts) == 2
                # Format: (lower, init) - no upper bound
                lower = parts[1] == "" || lowercase(parts[1]) == "-inf" ? -Inf : parse(Float64, parts[1])
                init = parse(Float64, parts[2])
                push!(thetas, THETASpec(lower, init, Inf, fixed))
            elseif length(parts) == 3
                lower = parts[1] == "" || lowercase(parts[1]) == "-inf" ? -Inf : parse(Float64, parts[1])
                init = parse(Float64, parts[2])
                upper = parts[3] == "" || lowercase(parts[3]) == "inf" ? Inf : parse(Float64, parts[3])
                push!(thetas, THETASpec(lower, init, upper, fixed))
            end

            # Move past FIX if present
            if fixed
                i = j + 1
                while i <= length(text) && (isspace(text[i]) || uppercase(text[i]) in ['F','I','X'])
                    i += 1
                end
            else
                i = j + 1
            end
        else
            # Simple number format
            j = i
            while j <= length(text) && (isdigit(text[j]) || text[j] in ['.', '-', '+', 'E', 'e'])
                j += 1
            end
            if j > i
                num_str = text[i:j-1]
                try
                    init = parse(Float64, num_str)
                    push!(thetas, THETASpec(init))
                catch
                end
            end
            i = j
        end
    end

    return thetas
end

"""
Parse \$OMEGA section.
"""
function _parse_omega(lines::Vector{String})::Vector{OMEGABlock}
    if isempty(lines)
        return OMEGABlock[]
    end

    text = join(lines, " ")
    text = replace(text, r"\s+" => " ")

    blocks = OMEGABlock[]
    structure = :diagonal
    fixed = false

    # Check for BLOCK or DIAGONAL
    if occursin(r"BLOCK"i, text)
        structure = :block
        # Extract block size if specified
        m = match(r"BLOCK\s*\((\d+)\)"i, text)
        if m !== nothing
            # block_size = parse(Int, m.captures[1])
        end
        text = replace(text, r"BLOCK\s*\(\d+\)"i => "")
    end

    if occursin(r"DIAGONAL"i, text)
        structure = :diagonal
        text = replace(text, r"DIAGONAL\s*\(\d+\)"i => "")
        text = replace(text, r"DIAGONAL"i => "")
    end

    if occursin(r"\bFIX\b"i, text)
        fixed = true
        text = replace(text, r"\bFIX\b"i => "")
    end

    # Extract numeric values
    values = Float64[]
    for m in eachmatch(r"[\d.eE+-]+", text)
        try
            push!(values, parse(Float64, m.match))
        catch
        end
    end

    if !isempty(values)
        push!(blocks, OMEGABlock(values, structure, fixed))
    end

    return blocks
end

"""
Parse \$SIGMA section.
"""
function _parse_sigma(lines::Vector{String})::Vector{SIGMABlock}
    if isempty(lines)
        return SIGMABlock[]
    end

    text = join(lines, " ")
    text = replace(text, r"\s+" => " ")

    blocks = SIGMABlock[]
    structure = :diagonal
    fixed = false

    if occursin(r"BLOCK"i, text)
        structure = :block
        text = replace(text, r"BLOCK\s*\(\d+\)"i => "")
    end

    if occursin(r"\bFIX\b"i, text)
        fixed = true
        text = replace(text, r"\bFIX\b"i => "")
    end

    # Extract numeric values
    values = Float64[]
    for m in eachmatch(r"[\d.eE+-]+", text)
        try
            push!(values, parse(Float64, m.match))
        catch
        end
    end

    if !isempty(values)
        push!(blocks, SIGMABlock(values, structure, fixed))
    end

    return blocks
end

"""
Parse \$ESTIMATION section.
"""
function _parse_estimation(lines::Vector{String})::Dict{String,Any}
    if isempty(lines)
        return Dict{String,Any}()
    end

    text = uppercase(join(lines, " "))
    result = Dict{String,Any}()

    # Common estimation options
    if occursin("METHOD=0", text) || occursin("METHOD=ZERO", text)
        result["method"] = "FO"
    elseif occursin("METHOD=1", text) || occursin("METHOD=COND", text)
        result["method"] = "FOCE"
    elseif occursin("LAPLACIAN", text) || occursin("LAPLACE", text)
        result["method"] = "LAPLACIAN"
    elseif occursin("SAEM", text)
        result["method"] = "SAEM"
    end

    if occursin("INTERACTION", text) || occursin("INTER", text)
        result["interaction"] = true
    end

    # Parse MAXEVAL
    m = match(r"MAXEVAL\s*=\s*(\d+)", text)
    if m !== nothing
        result["maxeval"] = parse(Int, m.captures[1])
    end

    return result
end

"""
Parse \$TABLE section.
"""
function _parse_tables(lines::Vector{String})::Vector{Dict{String,Any}}
    if isempty(lines)
        return Dict{String,Any}[]
    end

    # For now, just store raw table specifications
    tables = Dict{String,Any}[]
    push!(tables, Dict{String,Any}("raw" => join(lines, " ")))

    return tables
end
