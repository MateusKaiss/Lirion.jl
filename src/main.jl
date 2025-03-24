function main(matpower_database_path::String, model_type::Type{<:OPFModel})
    matpower_data = open(matpower_database_path) do io
        _parse_matpower_string(read(io, String))
        return nothing
    end

    return results
end
