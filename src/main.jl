function main(ARGS)
    args = Args(ARGS)

    @info "Caso: $(args.database_path)"
    @info "Modelo: $(args.model_type)"
    @info "Otimizador: $(args.optimizer)"

    data = if endswith(lowercase(args.database_path), ".m")
        open(args.database_path) do io
            return parse_matpower_string_silent(read(io, String))
        end
    elseif isdir(args.database_path)
        parse_case_from_csv_silent(args.database_path)
    else
        error("Unrecognized input format. Provide a .m file or a directory with .csv files.")
    end

    add_case_index_maps!(data)

    result_model = args.model_type == SOCP ?
        try_with_fallback_scale(data, args.model_type; optimizer = args.optimizer) :
        build_and_optimize_model(data, args.model_type; optimizer = args.optimizer)

    result_dict = build_result_dict(result_model, data)

    return result_dict, result_model, data
end

function julia_main()::Cint
    main(ARGS)
    return 0
end