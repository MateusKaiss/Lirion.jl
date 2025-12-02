function run_lirion(inputs::LirionInputs)
    @info "Case: $(inputs.database_path)"
    @info "Model: $(inputs.model_type)"
    @info "Optimizer (raw): $(inputs.optimizer)"

    if is_distributed_model(inputs.model_type)
        result = run_distributed(inputs)
        return result
    else
        result_model = solve_model(inputs)
        sys  = inputs.system
        data = sys.raw

        result_dict = build_result_dict(result_model, data)

        if inputs.write_outputs
            save_result_to_csv(result_dict, data, inputs.database_path)
        end

        return result_dict, result_model, sys
    end
end

function main(ARGS)
    args   = Args(ARGS)
    inputs = build_inputs(args)
    return run_lirion(inputs)
end

function julia_main()::Cint
    main(ARGS)
    return 0
end