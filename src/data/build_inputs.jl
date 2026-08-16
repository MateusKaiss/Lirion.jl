# ============================================================
# Inputs
# ============================================================

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

function load_inputs(path::AbstractString)::CaseRaw
    params = abspath(String(path))
    return _load(params, detect_input(params))
end

# ============================================================
# Inputs builder
# ============================================================

# ------------------------------------------------------------
# Path entry (shared)
# ------------------------------------------------------------

function _build_inputs(
    path::String,
    algorithm::AbstractAlgorithm;
    model::Type{<:OPFModel},
    optimizer::Any,
    outputs::Bool,
)::AbstractLirionInputs
    return _build_inputs(
        load_inputs(path), algorithm;
        model, optimizer, outputs, label = path,
    )
end

# ------------------------------------------------------------
# Centralized
# ------------------------------------------------------------

function _build_inputs(
    raw_case::CaseRaw,
    ::CentralizedAlgorithm;
    model::Type{<:OPFModel},
    optimizer::Any,
    outputs::Bool,
    label::String = raw_case.name,
)::CentralizedInputs
    preprocess_case_data!(raw_case, model)
    return CentralizedInputs(label, model, optimizer, outputs, build_system_data(raw_case))
end

# ------------------------------------------------------------
# Distributed
# ------------------------------------------------------------

_make_config(algorithm::ADMMAlgorithm, model::Type{<:OPFModel}) =
    ADMMConfig(
        model,
        algorithm.coupling,
        algorithm.sectioning,
        algorithm.max_iter;
        algorithm.kwargs...,
    )
_make_config(algorithm::ALADINAlgorithm, model::Type{<:OPFModel}) =
    ALADINConfig(
        model,
        algorithm.coupling,
        algorithm.sectioning,
        algorithm.max_iter;
        algorithm.kwargs...,
    )

function _build_inputs(
    raw_case::CaseRaw,
    algorithm::AbstractDistributedAlgorithm;
    model::Type{<:OPFModel},
    optimizer::Any,
    outputs::Bool,
    label::String = raw_case.name,
)::DistributedInputs
    preprocess_case_data!(raw_case, model)
    cfg = _make_config(algorithm, model)
    dist_sys = build_distributed_system(raw_case, algorithm.sectioning)
    return DistributedInputs(label, cfg, optimizer, outputs, dist_sys)
end
