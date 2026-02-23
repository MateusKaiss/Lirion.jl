# ==============================================================================
# Distributed system
# ==============================================================================

struct DistributedSystem
    zones::Vector{Int}
    parts::Vector{SystemData}

    zone_branch_idx_in_full::Vector{Vector{Int}}
    branch_idx_full_to_zone::Vector{Dict{Int, Int}}
    shared_buses_zone::Vector{Vector{Int}}

    shared_buses::Vector{Int}
    interzonal_branches::Vector{Int}
end

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
        model,
        algorithm.coupling,
        algorithm.sectioning,
        algorithm.max_iter;
        algorithm.kwargs...,
    )

    parts_dict = prepare_distributed_partitions(raw, algorithm.sectioning)
    interzonal_branches = _find_interzonal_branches(raw, algorithm.sectioning)

    zones = sort!(collect(keys(parts_dict)))
    n_zones = length(zones)

    parts = Vector{SystemData}(undef, n_zones)
    zone_branch_idx_in_full = Vector{Vector{Int}}(undef, n_zones)
    branch_idx_full_to_zone = Vector{Dict{Int, Int}}(undef, n_zones)
    shared_buses_zone = Vector{Vector{Int}}(undef, n_zones)

    @inbounds for i in 1:n_zones
        z = zones[i]
        zone_raw = parts_dict[z]::CaseRaw

        parts[i] = build_system_data(zone_raw)

        zone_branch_idx_in_full[i] = zone_raw.dist.zone_branch_idx_in_full
        branch_idx_full_to_zone[i] = zone_raw.dist.branch_idx_full_to_zone
        shared_buses_zone[i] = zone_raw.dist.zone_shared_bus_ids
    end

    shared_buses = raw.dist.full_shared_bus_ids

    admm_sys = DistributedSystem(
        zones,
        parts,
        zone_branch_idx_in_full,
        branch_idx_full_to_zone,
        shared_buses_zone,
        shared_buses,
        interzonal_branches,
    )

    return LirionInputs(path, cfg, optimizer, outputs, nothing, admm_sys)
end
