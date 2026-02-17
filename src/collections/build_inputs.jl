# ============================================================
# Core structs
# ============================================================

struct ADMMSystem
    zones::Vector{Int}
    zone_to_idx::Dict{Int, Int}
    parts::Vector{SystemData}
    shared_buses::Vector{Int}
    interzonal_branches::Vector{Int}
end

@inline get_zone(sys::ADMMSystem, zone::Int) = @inbounds sys.parts[sys.zone_to_idx[zone]]

struct LirionInputs
    database_path::String
    model_type::Union{Type{<:OPFModel}, DistributedOPFModel}
    optimizer::Any
    write_outputs::Bool
    system::Union{Nothing, SystemData}
    admm_system::Union{Nothing, ADMMSystem}
end

# ============================================================
# IO helpers
# ============================================================

function parse_case_from_csv_silent(dir::String)
    with_logger(NullLogger()) do
        return parse_case_from_csv(dir)
    end
end

function parse_matpower_string_silent(data_string::String)
    with_logger(NullLogger()) do
        return _parse_matpower_string(data_string)
    end
end

function load_input_data(path::String)
    result = if endswith(lowercase(path), ".m")
        open(path) do io
            return parse_matpower_string_silent(read(io, String))
        end
    elseif isdir(path)
        parse_case_from_csv_silent(path)
    else
        error("Unrecognized input format. Use .m or folder with .csv files.")
    end
    global_logger(ConsoleLogger(stderr, Logging.Info))
    return result
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
    data = load_input_data(path)
    preprocess_case_data!(data, OPFModel)

    full_sys = build_system_data(data)
    return LirionInputs(path, model, optimizer, outputs, full_sys, nothing)
end

function _build_inputs(
    path::String,
    algorithm::ADMMAlgorithmSpec;
    model::Type{<:OPFModel},
    optimizer::Any,
    outputs::Bool,
)::LirionInputs
    data = load_input_data(path)
    preprocess_case_data!(data, OPFModel)

    cfg = DistributedOPFModel(
        model,
        algorithm.coupling,
        algorithm.sectioning,
        algorithm.max_iter;
        algorithm.kwargs...,
    )

    parts_dict, shared_buses = prepare_admm_partitions(data)
    interzonal_branches = _find_interzonal_branches(data)

    zones = sort!(collect(keys(parts_dict)))
    zone_to_idx = Dict{Int, Int}(z => i for (i, z) in pairs(zones))

    parts = Vector{SystemData}(undef, length(zones))
    @inbounds for (i, z) in pairs(zones)
        zone_data = parts_dict[z]::Dict{String, Any}

        admm_meta = ADMMZoneMeta(
            zone_data["branch_local_indices"]::Vector{Int},
            zone_data["branch_global_to_local"]::Dict{Int, Int},
            zone_data["shared_buses_zone"]::Vector{Int},
        )

        parts[i] = build_system_data(zone_data; admm_meta = admm_meta)
    end

    admm_sys = ADMMSystem(zones, zone_to_idx, parts, shared_buses, interzonal_branches)
    return LirionInputs(path, cfg, optimizer, outputs, nothing, admm_sys)
end
