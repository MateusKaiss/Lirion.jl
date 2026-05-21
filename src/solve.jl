@inline function _assert_valid_path(path::AbstractString)::String
    p = abspath(String(path))
    (isfile(p) || isdir(p)) ||
        throw(ArgumentError("Invalid path: $p is not a file or directory."))
    return p
end

function run_lirion(inputs::CentralizedInputs)::LirionResult
    @info "Case: $(inputs.database_path)"
    @info "Model: $(inputs.model_type)"
    @info "Optimizer (raw): $(inputs.optimizer)"

    t0 = time_ns()
    result_model = solve_model(inputs)
    sys = inputs.system
    return build_lirion_result_centralized(
        result_model, sys, inputs.model_type, t0,
    )
end

function run_lirion(inputs::DistributedInputs)::LirionResult
    @info "Case: $(inputs.database_path)"
    @info "Model: $(inputs.model_type)"
    @info "Optimizer (raw): $(inputs.optimizer)"

    t0 = time_ns()
    partitions, ws, iter_done = solve_distributed(inputs)
    return build_lirion_result_distributed(
        partitions, ws, inputs.model_type, iter_done, t0, inputs.system,
    )
end

function solve(
    path::AbstractString;
    model::Type{<:OPFModel},
    algorithm::AbstractAlgorithm = centralized(),
    optimizer = nothing,
    outputs::Bool = false,
    plots::Bool = false,
)::LirionResult
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

function solve(
    case::CaseRaw;
    model::Type{<:OPFModel},
    algorithm::AbstractAlgorithm = centralized(),
    optimizer = nothing,
    outputs::Bool = false,
)::LirionResult
    inputs = _build_inputs(
        case,
        algorithm;
        model = model,
        optimizer = optimizer,
        outputs = outputs,
    )
    return run_lirion(inputs)
end
