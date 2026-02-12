function run_lirion(inputs::LirionInputs)::LirionRunResult
    @info "Case: $(inputs.database_path)"
    @info "Model: $(inputs.model_type)"
    @info "Optimizer (raw): $(inputs.optimizer)"

    if is_distributed_model(inputs.model_type)
        return run_distributed(inputs)
    else
        result_model = solve_model(inputs)

        sys = inputs.system
        data = sys.raw

        result_dict = build_result_dict(result_model, data)

        if inputs.write_outputs
            save_result_to_csv(result_dict, data, inputs.database_path)
        end

        return LirionRunResult(;
            result = result_dict,
            model = result_model,
            system = sys,
            diagnostics = NoDiagnostics(),
        )
    end
end

function solve(
    path::AbstractString;
    model,
    algorithm = centralized(),
    optimizer = nothing,
    outputs::Bool = false,
    plots::Bool = false,
)::LirionRunResult

    path_str = abspath(String(path))
    if !(isfile(path_str) || isdir(path_str))
        error("Invalid path: $path_str is not a file or directory.")
    end

    inputs = _build_inputs(
        path_str;
        model = model,
        algorithm = algorithm,
        optimizer = optimizer,
        outputs = outputs,
    )

    out = run_lirion(inputs)

    _maybe_build_plots!(plots, path_str, inputs.model_type, out)

    return out
end

# centralized 

# out = Lirion.solve(
#     file;
#     model = Lirion.ac(),              # ou dc(), socp()
#     algorithm = Lirion.centralized(), # opcional (padrão)
#     optimizer = Ipopt.Optimizer,
#     outputs = true,
#     plots = false,
# )
# obj = out.result["objective"]
