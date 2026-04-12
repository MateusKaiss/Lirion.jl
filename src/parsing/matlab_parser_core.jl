# ==============================================================================
# Lirion.jl — Matlab/Matpower minimal parser (adapted from InfrastructureModels)
#
# LICENSE NOTE:
# Portions of this parser are adapted from InfrastructureModels.jl (LANL/Triad).
# Keep a copy of the original license text in:
#   licenses/InfrastructureModels_LICENSE.txt
# and clearly mark modifications (this file is modified/adapted).
# ==============================================================================

# ------------------------------------------------------------------------------
# Result type
# ------------------------------------------------------------------------------

struct MatlabParseResult
    data::Dict{String, Any}
    function_name::Union{Nothing, String}
    column_names::Dict{String, Any}
end

# ------------------------------------------------------------------------------
# Base conversions
# ------------------------------------------------------------------------------

function check_type(typ, value)
    if isa(value, typ)
        return value
    elseif isa(value, String) || isa(value, SubString)
        try
            return parse(typ, String(value))
        catch e
            @error "Parsing error: string \"$(value)\" cannot be parsed to $(typ)"
            rethrow(e)
        end
    else
        try
            return typ(value)
        catch e
            @error "Parsing error: value $(value) of type $(typeof(value)) cannot be converted to $(typ)"
            rethrow(e)
        end
    end
end

# ------------------------------------------------------------------------------
# Assignment parsing + scalar typing
# ------------------------------------------------------------------------------

# Extract `(name, value)` from a MATLAB assignment like `foo.bar = 123;`.
# Uses substring views + findfirst instead of `split`, avoiding per-line Vector
# allocations.
function _extract_matlab_assignment(s::AbstractString)
    # Trim at first ';' (if present).
    semi = findfirst(==(';'), s)
    statement = semi === nothing ? s : SubString(s, firstindex(s), prevind(s, semi))

    # Split at first '='.
    eq = findfirst(==('='), statement)
    @assert eq !== nothing "expected '=' in assignment: $(s)"

    name = strip(SubString(statement, firstindex(statement), prevind(statement, eq)))
    value = strip(SubString(statement, nextind(statement, eq), lastindex(statement)))
    return name, value
end

@inline function _is_quoted_literal(s::AbstractString)
    isempty(s) && return false
    c1 = s[begin]
    return (c1 == '\'' || c1 == '"') && lastindex(s) > firstindex(s) && s[end] == c1
end

@inline function _unquote_matlab(s::AbstractString)
    c = s[begin]
    inner = s[nextind(s, firstindex(s)):prevind(s, lastindex(s))]
    if c == '\''
        return replace(String(inner), "''" => "'")
    else
        return replace(String(inner), "\"\"" => "\"")
    end
end

function _type_value(value_string::AbstractString)
    s = strip(value_string)

    if _is_quoted_literal(s)
        return _unquote_matlab(s)
    end

    if s == "Inf" || s == "+Inf"
        return Inf
    elseif s == "-Inf"
        return -Inf
    elseif s == "NaN"
        return NaN
    end

    if occursin(".", s) || occursin("e", lowercase(s))
        return check_type(Float64, s)
    else
        return check_type(Int, s)
    end
end

# ------------------------------------------------------------------------------
# Line normalization helpers
# ------------------------------------------------------------------------------

function _strip_comment_matlab(s::AbstractString)
    in_quote = false
    quote_char = '\0'
    i = firstindex(s)

    while i <= lastindex(s)
        c = s[i]

        if in_quote
            if c == quote_char
                j = nextind(s, i)
                if j <= lastindex(s) && s[j] == quote_char
                    i = nextind(s, j)
                    continue
                end
                in_quote = false
            end
        else
            if c == '\'' || c == '"'
                in_quote = true
                quote_char = c
            elseif c == '%'
                return s[firstindex(s):prevind(s, i)]
            end
        end

        i = nextind(s, i)
    end

    return s
end

function _clean_line(raw::AbstractString)
    line = strip(replace(raw, "\r" => ""))
    isempty(line) && return nothing
    startswith(line, "%") && return nothing

    line = strip(_strip_comment_matlab(line))
    isempty(line) && return nothing

    return line
end

function _add_line_delimiter(mp_line::AbstractString, start_char, end_char)
    line = mp_line
    strip(line) == string(start_char) && return line

    if !occursin(";", line) && !occursin(string(end_char), line)
        line = "$(line);"
    end

    if occursin(string(end_char), line)
        prefix = strip(split(line, end_char)[1])
        if !isempty(prefix) && !occursin(";", prefix)
            line = replace(line, end_char => ";$(end_char)")
        end
    end

    return line
end

# ------------------------------------------------------------------------------
# Tokenization helpers (matrix rows)
# ------------------------------------------------------------------------------

function _split_line(mp_line::AbstractString)
    tokens = String[]
    curr = IOBuffer()
    in_quote = false
    quote_char = '\0'
    i = firstindex(mp_line)

    @inline isquote(c) = (c == '\'' || c == '"')

    function flush!()
        s = strip(String(take!(curr)))
        !isempty(s) && push!(tokens, s)
        return nothing
    end

    while i <= lastindex(mp_line)
        c = mp_line[i]

        if in_quote
            if c == quote_char
                j = nextind(mp_line, i)
                if j <= lastindex(mp_line) && mp_line[j] == quote_char
                    write(curr, c)
                    write(curr, mp_line[j])
                    i = nextind(mp_line, j)
                    continue
                else
                    write(curr, c)
                    in_quote = false
                    i = nextind(mp_line, i)
                    continue
                end
            else
                write(curr, c)
                i = nextind(mp_line, i)
                continue
            end
        else
            if isquote(c)
                in_quote = true
                quote_char = c
                write(curr, c)
                i = nextind(mp_line, i)
                continue
            end

            if isspace(c)
                flush!()
                while i <= lastindex(mp_line) && isspace(mp_line[i])
                    i = nextind(mp_line, i)
                end
                continue
            else
                write(curr, c)
                i = nextind(mp_line, i)
                continue
            end
        end
    end

    flush!()
    return tokens
end

# ------------------------------------------------------------------------------
# Matrix/cell block parsing engine
# ------------------------------------------------------------------------------

function _parse_matlab_data(lines, index, start_char, end_char)
    last_index = length(lines)
    line_count = 0
    columns = -1

    @assert occursin("=", lines[index+line_count])
    first_line = _clean_line(lines[index+line_count])
    @assert first_line !== nothing

    matrix_assignment = strip(first_line)

    @assert occursin(".", matrix_assignment)
    matrix_assignment_parts = split(matrix_assignment, '=')
    matrix_name = strip(matrix_assignment_parts[1])

    matrix_assignment_rhs = ""
    if length(matrix_assignment_parts) > 1
        matrix_assignment_rhs = strip(matrix_assignment_parts[2])
    end

    line_count += 1
    matrix_body_lines = String[matrix_assignment_rhs]
    found_close = occursin(string(end_char), matrix_assignment_rhs)

    while index + line_count <= last_index && !found_close
        line = _clean_line(lines[index+line_count])
        if line === nothing
            line_count += 1
            continue
        end

        occursin(string(end_char), line) && (found_close = true)
        push!(matrix_body_lines, line)
        line_count += 1
    end

    matrix_body_lines =
        [_add_line_delimiter(line, start_char, end_char) for line in matrix_body_lines]

    matrix_body = join(matrix_body_lines, ' ')
    matrix_body = strip(matrix_body)
    matrix_body = strip(strip(matrix_body), start_char)
    matrix_body = replace(matrix_body, "$(end_char);" => "")
    matrix_body = strip(matrix_body)

    matrix_body_rows = split(matrix_body, ';')
    while !isempty(matrix_body_rows) && isempty(strip(matrix_body_rows[end]))
        pop!(matrix_body_rows)
    end

    if !isempty(matrix_body_rows)
        lastrow = matrix_body_rows[end]
        if occursin(string(end_char), lastrow)
            matrix_body_rows[end] = strip(replace(lastrow, string(end_char) => ""))
            isempty(strip(matrix_body_rows[end])) && pop!(matrix_body_rows)
        end
    end

    matrix = Vector{Vector{String}}()
    for row in matrix_body_rows
        row_s = strip(row)
        isempty(row_s) && continue
        row_items = _split_line(row_s)
        push!(matrix, row_items)
        if columns < 0
            columns = length(row_items)
        elseif columns != length(row_items)
            @error "Matrix parsing error: inconsistent number of items in each row\n$(row_s)"
        end
    end

    if isempty(matrix)
        return Dict(
            "name" => matrix_name,
            "data" => Vector{Vector{String}}(),
            "line_count" => line_count,
        )
    end

    matrix_dict = Dict("name" => matrix_name, "data" => matrix, "line_count" => line_count)

    if index > 1 && occursin("%column_names%", lines[index-1])
        col_line = replace(lines[index-1], "%column_names%" => "")
        col_line = strip(col_line)
        colnames = split(col_line)

        if !isempty(matrix) && length(matrix[1]) != length(colnames)
            @error "Column name parsing error: data cols $(length(matrix[1])) vs column names $(length(colnames))"
        end

        matrix_dict["column_names"] = colnames
    end

    return matrix_dict
end

_parse_matlab_cells(lines, index) = _parse_matlab_data(lines, index, '{', '}')
_parse_matlab_matrix(lines, index) = _parse_matlab_data(lines, index, '[', ']')

# ------------------------------------------------------------------------------
# Dispatch helpers for main loop
# ------------------------------------------------------------------------------

@inline _is_function_decl(line::AbstractString) = startswith(line, "function")
@inline _is_assignment(line::AbstractString) = occursin('=', line)

@inline function _rhs_starts_matrix(line::AbstractString)
    p = findfirst(==('='), line)
    p === nothing && return false
    q = findnext(==('['), line, nextind(line, p))
    return q !== nothing
end

@inline function _rhs_starts_cells(line::AbstractString)
    p = findfirst(==('='), line)
    p === nothing && return false
    q = findnext(==('{'), line, nextind(line, p))
    return q !== nothing
end

@inline function _is_unexpected_struct_assignment(
    line::AbstractString,
    struct_name::Union{Nothing, String},
)
    struct_name === nothing && return false
    occursin('.', line) || return false
    return !startswith(line, string(struct_name, '.'))
end

# ------------------------------------------------------------------------------
# Handlers
# ------------------------------------------------------------------------------

function _handle_function_line!(
    line::AbstractString,
    struct_name_ref::Base.RefValue{Union{Nothing, String}},
    function_name_ref::Base.RefValue{Union{Nothing, String}},
)
    lhs, rhs = _extract_matlab_assignment(line)
    struct_name_ref[] = strip(replace(lhs, "function" => ""))
    function_name_ref[] = strip(rhs)
    return nothing
end

function _handle_scalar_assignment!(line::AbstractString, data::Dict{String, Any})
    name, value = _extract_matlab_assignment(line)
    data[name] = _type_value(value)
    return nothing
end

function _handle_matrix_assignment!(
    lines,
    i::Int,
    data::Dict{String, Any},
    column_names::Dict{String, Any},
)
    md = _parse_matlab_matrix(lines, i)
    data[md["name"]] = md["data"]
    if haskey(md, "column_names")
        column_names[md["name"]] = md["column_names"]
    end
    return md["line_count"]
end

function _handle_cells_assignment!(
    lines,
    i::Int,
    data::Dict{String, Any},
    column_names::Dict{String, Any},
)
    cd = _parse_matlab_cells(lines, i)
    data[cd["name"]] = cd["data"]
    if haskey(cd, "column_names")
        column_names[cd["name"]] = cd["column_names"]
    end
    return cd["line_count"]
end

function parse_matlab_file(path::String; kwargs...)
    open(path, "r") do io
        return parse_matlab_string(read(io, String); kwargs...)
    end
end

# ============================================================
# Main entrypoint
# ============================================================

function parse_matlab_string(data_string::AbstractString; meta::Bool = false)
    lines = split(data_string, '\n')

    data = Dict{String, Any}()
    column_names = Dict{String, Any}()

    struct_name = Ref{Union{Nothing, String}}(nothing)
    function_name = Ref{Union{Nothing, String}}(nothing)

    i = 1
    last = length(lines)

    while i <= last
        line = _clean_line(lines[i])
        if line === nothing
            i += 1
            continue
        end

        if _is_function_decl(line)
            _handle_function_line!(line, struct_name, function_name)
            i += 1
            continue
        end

        if !_is_assignment(line)
            @debug "Skipping line in matlab parser: $line"
            i += 1
            continue
        end

        if _is_unexpected_struct_assignment(line, struct_name[])
            @debug "Expected assignments to $(struct_name[]), got: $line"
        end

        if _rhs_starts_matrix(line)
            consumed = _handle_matrix_assignment!(lines, i, data, column_names)
            i += consumed
            continue
        end

        if _rhs_starts_cells(line)
            consumed = _handle_cells_assignment!(lines, i, data, column_names)
            i += consumed
            continue
        end

        _handle_scalar_assignment!(line, data)
        i += 1
    end

    return meta ? MatlabParseResult(data, function_name[], column_names) : data
end
