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
# Public API
# ------------------------------------------------------------------------------

function parse_matlab_string(data_string::String; extended::Bool=false)
    data_lines = split(data_string, '\n')

    matlab_dict = Dict{String,Any}()
    struct_name = nothing
    function_name = nothing
    column_names = Dict{String,Any}()

    last_index = length(data_lines)
    index = 1
    while index <= last_index
        line = strip(replace(data_lines[index], "\r" => ""))
        isempty(line) && (index += 1; continue)
        startswith(line, "%") && (index += 1; continue)

        # strip inline comment
        line_nocom = strip(split(line, "%")[1])
        isempty(line_nocom) && (index += 1; continue)

        if occursin("function", line_nocom)
            func, value = _extract_matlab_assignment(line_nocom)
            struct_name = strip(replace(func, "function" => ""))
            function_name = value
        elseif occursin("=", line_nocom)
            if struct_name !== nothing && !occursin("$(struct_name).", line_nocom)
                @warn "Assignments expected to be made to \"$(struct_name)\" but got: $(line_nocom)"
            end

            if occursin("[", line_nocom)
                matrix_dict = _parse_matlab_matrix(data_lines, index)
                matlab_dict[matrix_dict["name"]] = matrix_dict["data"]
                if haskey(matrix_dict, "column_names")
                    column_names[matrix_dict["name"]] = matrix_dict["column_names"]
                end
                index = index + matrix_dict["line_count"] - 1

            elseif occursin("{", line_nocom)
                cell_dict = _parse_matlab_cells(data_lines, index)
                matlab_dict[cell_dict["name"]] = cell_dict["data"]
                if haskey(cell_dict, "column_names")
                    column_names[cell_dict["name"]] = cell_dict["column_names"]
                end
                index = index + cell_dict["line_count"] - 1

            else
                name, value = _extract_matlab_assignment(line_nocom)
                matlab_dict[name] = _type_value(value)
            end
        else
            @warn "Matlab parser skipping line: $(line_nocom)"
        end

        index += 1
    end

    return extended ? (matlab_dict, function_name, column_names) : matlab_dict
end

function parse_matlab_file(path::String; kwargs...)
    open(path, "r") do io
        return parse_matlab_string(read(io, String); kwargs...)
    end
end

# ------------------------------------------------------------------------------
# Internals (minimal subset)
# ------------------------------------------------------------------------------

function _extract_matlab_assignment(string::AbstractString)
    statement = split(string, ';')[1]
    statement_parts = split(statement, '=')
    @assert length(statement_parts) == 2
    name = strip(statement_parts[1])
    value = strip(statement_parts[2])
    return name, value
end

function _type_value(value_string::AbstractString)
    value_string = strip(value_string)

    if occursin("'", value_string) # string
        return strip(value_string, '\'')
    end

    # handle Inf/-Inf/NaN (Matlab tokens)
    if value_string == "Inf"
        return Inf
    elseif value_string == "-Inf"
        return -Inf
    elseif value_string == "NaN"
        return NaN
    end

    # float or int
    if occursin(".", value_string) || occursin("e", lowercase(value_string))
        return check_type(Float64, value_string)
    else
        return check_type(Int, value_string)
    end
end

function _type_array(string_array::Vector{T}) where {T<:AbstractString}
    ss = [strip(s) for s in string_array]

    if any(occursin("'", s) for s in ss)
        return [strip(s, '\'') for s in ss]
    elseif any(occursin(".", s) || occursin("e", lowercase(s)) || occursin("Inf", s) || occursin("NaN", s) for s in ss)
        return [check_type(Float64, s) for s in ss]
    else
        return [check_type(Int, s) for s in ss]
    end
end

_parse_matlab_cells(lines, index)  = _parse_matlab_data(lines, index, '{', '}')
_parse_matlab_matrix(lines, index) = _parse_matlab_data(lines, index, '[', ']')

function _parse_matlab_data(lines, index, start_char, end_char)
    last_index = length(lines)
    line_count = 0
    columns = -1

    @assert occursin("=", lines[index + line_count])
    matrix_assignment = split(lines[index + line_count], '%')[1] |> strip

    @assert occursin(".", matrix_assignment) # expects mpc.bus, mpc.gen, etc.
    matrix_assignment_parts = split(matrix_assignment, '=')
    matrix_name = strip(matrix_assignment_parts[1])

    matrix_assignment_rhs = ""
    if length(matrix_assignment_parts) > 1
        matrix_assignment_rhs = strip(matrix_assignment_parts[2])
    end

    line_count += 1
    matrix_body_lines = String[matrix_assignment_rhs]
    found_close = occursin(string(end_char), matrix_assignment_rhs)

    while index + line_count < last_index && !found_close
        line = strip(replace(lines[index + line_count], "\r" => ""))
        if isempty(line) || startswith(line, "%")
            line_count += 1
            continue
        end

        line = strip(split(line, "%")[1])
        occursin(string(end_char), line) && (found_close = true)
        push!(matrix_body_lines, line)
        line_count += 1
    end

    matrix_body_lines = [_add_line_delimiter(line, start_char, end_char) for line in matrix_body_lines]
    matrix_body = join(matrix_body_lines, ' ')
    matrix_body = strip(replace(strip(strip(matrix_body), start_char), "$(end_char);" => ""))
    matrix_body_rows = split(matrix_body, ';')
    matrix_body_rows = matrix_body_rows[1:(length(matrix_body_rows)-1)]

    matrix = Vector{Vector{Any}}()
    for row in matrix_body_rows
        row_items = split_line(strip(row))
        push!(matrix, row_items)
        if columns < 0
            columns = length(row_items)
        elseif columns != length(row_items)
            @error "Matrix parsing error: inconsistent number of items in each row\n$(row)"
        end
    end

    rows = length(matrix)
    typed_columns = [_type_array([string(matrix[r][c]) for r in 1:rows]) for c in 1:columns]
    for r in 1:rows
        matrix[r] = Any[typed_columns[c][r] for c in 1:columns]
    end

    matrix_dict = Dict("name" => matrix_name, "data" => matrix, "line_count" => line_count)

    # optional column names line: "%column_names% ..."
    if index > 1 && occursin("%column_names%", lines[index-1])
        column_names_string = replace(lines[index-1], "%column_names%" => "")
        colnames = split(column_names_string)
        if length(matrix[1]) != length(colnames)
            @error "Column name parsing error: data cols $(length(matrix[1])) vs column names $(length(colnames))"
        end
        matrix_dict["column_names"] = colnames
    end

    return matrix_dict
end

function split_line(mp_line::AbstractString)
    tokens = Any[]
    curr_token = ""
    is_curr_token_quote = false
    isquote(c) = (c == '\'' || c == '"')

    function _push_curr_token()
        curr_token_local = strip(curr_token)
        if !isempty(curr_token_local)
            push!(tokens, curr_token_local)
        end
        curr_token = ""
        curr_pos += 1
        is_curr_token_quote = false
        return nothing
    end

    function _push_curr_char()
        curr_token *= mp_line[curr_pos]
        curr_pos = nextind(mp_line, curr_pos)
        return nothing
    end

    curr_pos = 1
    while curr_pos <= lastindex(mp_line)
        if is_curr_token_quote
            if mp_line[curr_pos] == curr_token[1]
                if mp_line[prevind(mp_line, curr_pos)] == '\\'
                    _push_curr_char()
                elseif curr_pos < lastindex(mp_line) && mp_line[nextind(mp_line, curr_pos)] == curr_token[1]
                    curr_token *= mp_line[curr_pos]
                    curr_pos = nextind(mp_line, nextind(mp_line, curr_pos))
                else
                    _push_curr_token()
                end
            else
                _push_curr_char()
            end
        else
            if isspace(mp_line[curr_pos]) && curr_pos < lastindex(mp_line) && !isspace(mp_line[nextind(mp_line, curr_pos)])
                _push_curr_token()
            elseif isquote(mp_line[curr_pos])
                curr_token = strip(curr_token * mp_line[curr_pos])
                is_curr_token_quote = true
                curr_pos = nextind(mp_line, curr_pos)
            else
                _push_curr_char()
            end
        end
    end
    # push last
    curr_token_local = strip(curr_token)
    !isempty(curr_token_local) && push!(tokens, curr_token_local)
    return tokens
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