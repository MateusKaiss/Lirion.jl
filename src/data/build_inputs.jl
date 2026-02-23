# ==============================================================================
# Inputs
# ==============================================================================

struct LirionInputs
    database_path::String
    model_type::Union{Type{<:OPFModel}, DistributedOPFModel}
    optimizer::Any
    write_outputs::Bool
    system::Union{Nothing, SystemData}
    admm_system::Union{Nothing, DistributedSystem}
end

# ============================================================
# IO helpers
# ============================================================

@inline detect_input(path::AbstractString)::InputFormat =
    endswith(lowercase(path), ".m") ? MatpowerM() :
    isdir(path) ? CsvFolder() :
    throw(
        ArgumentError(
            "Unrecognized input format for path=$(abspath(String(path))). Use .m or a folder with .csv files.",
        ),
    )

@inline _load(path::AbstractString, matpower::MatpowerM)::CaseRaw =
    open(String(path)) do io
        return _parse_case(read(io, String), matpower)
    end

@inline _load(path::AbstractString, csv::CsvFolder)::CaseRaw =
    _parse_case(String(path), csv)

function load_input_data(path::AbstractString)::CaseRaw
    p = abspath(String(path))
    return _load(p, detect_input(p))
end

# ============================================================
# Inputs builder
# ============================================================

function _build_inputs(
    path::String,
    algorithm::CentralizedAlgorithm;
    model::Type{<:OPFModel},
    optimizer::Any,
    outputs::Bool,
)::LirionInputs
    raw = load_input_data(path)
    preprocess_case_data!(raw, OPFModel)

    full_sys = build_system_data(raw)
    return LirionInputs(path, model, optimizer, outputs, full_sys, nothing)
end

function _build_inputs(
    path::String,
    algorithm::ADMMAlgorithmSpec;
    model::Type{<:OPFModel},
    optimizer::Any,
    outputs::Bool,
)::LirionInputs
    raw = load_input_data(path)
    preprocess_case_data!(raw, OPFModel)

    cfg = DistributedOPFModel(
        model, algorithm.coupling, algorithm.sectioning, algorithm.max_iter;
        algorithm.kwargs...,
    )

    admm_sys = build_distributed_system(raw, algorithm.sectioning)

    return LirionInputs(path, cfg, optimizer, outputs, nothing, admm_sys)
end
