# ============================================================
# Input validation
# ============================================================

@inline function _assert_valid_path(path::AbstractString)::String
    p = abspath(String(path))
    (isfile(p) || isdir(p)) ||
        throw(ArgumentError("Invalid path: $p is not a file or directory."))
    return p
end

# ============================================================
# Internal runners (centralized / distributed dispatch)
# ============================================================

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

# ============================================================
# Public API — solve
# ============================================================

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
    plots::Bool = false,
)::LirionResult
    inputs = _build_inputs(
        case,
        algorithm;
        model = model,
        optimizer = optimizer,
        outputs = outputs,
    )
    out = run_lirion(inputs)
    # No source path here: plots are written relative to the case name,
    # i.e. into ./plots_outputs in the current working directory.
    _maybe_build_plots!(plots, inputs.database_path, inputs.model_type, out)
    return out
end
