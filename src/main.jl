function main(database_path::String, model_type::Type{<:OPFModel})
    data = nothing

    if endswith(lowercase(database_path), ".m")
        data = open(database_path) do io
            return _parse_matpower_string(read(io, String))
        end
    elseif isdir(database_path)
        data = parse_case_from_csv(database_path)
    else
        error(
            "Unrecognized input format. Please provide a .m file or a directory with .csv files.",
        )
    end

    result_model = build_and_optimize_model(data, model_type)
    result_dict = build_result_dict(result_model, data)

    return result_dict
end
