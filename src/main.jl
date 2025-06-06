function main(ARGS)
    args = Args(ARGS)

    @info "Case: $(args.database_path)"
    @info "Model: $(args.model_type)"
    @info "Optimizer: $(args.optimizer)"

    data = load_input_data(args.database_path)
    preprocess_case_data!(data)

    result_model = build_model(data, args.model_type; optimizer = args.optimizer)

    result_dict = build_result_dict(result_model, data)

    return result_dict, result_model, data
end

function load_input_data(path::String)
    if endswith(lowercase(path), ".m")
        return open(path) do io
            parse_matpower_string_silent(read(io, String))
        end
    elseif isdir(path)
        return parse_case_from_csv_silent(path)
    else
        error("Unrecognized input format. Use .m or folder with .csv files.")
    end
end


function julia_main()::Cint
    main(ARGS)
    return 0
end