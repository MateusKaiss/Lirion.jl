function main(matpower_database_path::String, model_type::Type{<:OPFModel})
    matpower_data = open(matpower_database_path) do io
        return _parse_matpower_string(read(io, String))
    end

    results = build_and_optimize_model(matpower_data, model_type)

    return results
end
