function main(
    database_path::String,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
)
    # Etapa 1: Carregamento dos dados
    data = if endswith(lowercase(database_path), ".m")
        open(database_path) do io
            return parse_matpower_string_silent(read(io, String))
        end
    elseif isdir(database_path)
        parse_case_from_csv_silent(database_path)
    else
        error("Unrecognized input format. Provide a .m file or a directory with .csv files.")
    end

    # Etapa 2: Pré-processamento
    add_case_index_maps!(data)

    # Etapa 3: Construção e otimização
    result_model =
        model_type == SOCP ?
        try_with_fallback_scale(data, model_type; optimizer = optimizer) :
        build_and_optimize_model(data, model_type; optimizer = optimizer)

    result_dict = build_result_dict(result_model, data)

    return result_dict, result_model, data
end
