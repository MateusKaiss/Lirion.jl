function main(
    database_path::String,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
)
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

    # Pré-processamento
    add_case_index_maps!(data)

    # Escolha da função de otimização
    result_model = if model_type == SOCP
        try_with_fallback_scale(data, model_type; optimizer = optimizer)
    else
        build_and_optimize_model(data, model_type; optimizer = optimizer)
    end

    result_dict = build_result_dict(result_model, data)

    return result_dict, result_model, data
end