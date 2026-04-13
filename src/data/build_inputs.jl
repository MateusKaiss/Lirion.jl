# ==============================================================================
# Inputs
# ==============================================================================

abstract type AbstractLirionInputs end

struct CentralizedInputs <: AbstractLirionInputs
    database_path::String
    model_type::Type{<:OPFModel}
    optimizer::Any
    write_outputs::Bool
    system::SystemData
end

struct DistributedInputs <: AbstractLirionInputs
    database_path::String
    model_type::AbstractDistributedConfig
    optimizer::Any
    write_outputs::Bool
    system::DistributedSystem
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
    params = abspath(String(path))
    return _load(params, detect_input(params))
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
)::CentralizedInputs
    raw_case = load_input_data(path)
    preprocess_case_data!(raw_case, OPFModel)

    full_sys = build_system_data(raw_case)

    return CentralizedInputs(
        path,
        model,
        optimizer,
        outputs,
        full_sys,
    )
end

function _build_inputs(
    path::String,
    algorithm::ADMMAlgorithm;
    model::Type{<:OPFModel},
    optimizer::Any,
    outputs::Bool,
)::DistributedInputs
    raw_case = load_input_data(path)
    preprocess_case_data!(raw_case, OPFModel)

    cfg = DistributedOPFModel(
        model, algorithm.coupling, algorithm.sectioning, algorithm.max_iter;
        algorithm.kwargs...,
    )

    dist_sys = build_distributed_system(raw_case, algorithm.sectioning)

    return DistributedInputs(
        path,
        cfg,
        optimizer,
        outputs,
        dist_sys,
    )
end

function _build_inputs(
    path::String,
    algorithm::ALADINAlgorithm;
    model::Type{<:OPFModel},
    optimizer::Any,
    outputs::Bool,
)::DistributedInputs
    raw_case = load_input_data(path)
    preprocess_case_data!(raw_case, OPFModel)

    cfg = ALADINConfig(
        model, algorithm.coupling, algorithm.sectioning, algorithm.max_iter;
        algorithm.kwargs...,
    )

    dist_sys = build_distributed_system(raw_case, algorithm.sectioning)

    return DistributedInputs(
        path,
        cfg,
        optimizer,
        outputs,
        dist_sys,
    )
end
