# ============================================================
# Run + Solve
# ============================================================

# -----------------------------
# Helpers
# -----------------------------

@inline function _assert_valid_path(path::AbstractString)::String
    p = abspath(String(path))
    (isfile(p) || isdir(p)) ||
        throw(ArgumentError("Invalid path: $p is not a file or directory."))
    return p
end

function run_lirion(inputs::LirionInputs)::LirionRunResult
    @info "Case: $(inputs.database_path)"
    @info "Model: $(inputs.model_type)"
    @info "Optimizer (case_raw): $(inputs.optimizer)"

    if is_distributed_model(inputs.model_type)
        return run_distributed(inputs)
    else
        result_model = solve_model(inputs)

        sys = inputs.system
        @assert sys !== nothing

        data = sys.case_raw
        result_dict = build_result_dict(result_model, sys)

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
    model::Type{<:OPFModel},
    algorithm::AbstractAlgorithm = centralized(),
    optimizer = nothing,
    outputs::Bool = false,
    plots::Bool = false,
)::LirionRunResult
    path_str = _assert_valid_path(path)

    inputs = _build_inputs(
        path_str,
        algorithm;
        model = model,
        optimizer = optimizer,
        outputs = outputs,
    )

    out = run_lirion(inputs)
    _maybe_build_plots!(plots, path_str, inputs.model_type, out)
    return out
end
