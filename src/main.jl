function main(ARGS)
    args = Args(ARGS)

    @info "Case: $(args.database_path)"
    @info "Model: $(args.model_type)"
    @info "Optimizer: $(args.optimizer)"

    # --- Load and preprocess input data ---
    data = load_input_data(args.database_path)
    preprocess_case_data!(data, OPFModel)

    if args.model_type isa DistributedOPFModel
        runner = get_method_runner(args.model_type)
        result = runner(args.model_type, data)
        return result
    end

    # --- Build and solve the optimization model ---
    result_model = build_model(data, args.model_type; optimizer = args.optimizer)

    # --- Post-process and organize results ---
    result_dict = build_result_dict(result_model, data)

    return result_dict, result_model, data
end

function julia_main()::Cint
    main(ARGS)
    return 0
end
