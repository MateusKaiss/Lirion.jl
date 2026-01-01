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

@inline zoneindex(sys::ADMMSystem, zone::Int) = sys.zone_to_idx[zone]
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

function build_inputs(args::Args)::LirionInputs
    data = load_input_data(args.database_path)
    preprocess_case_data!(data, OPFModel)

    if is_centralized_model(args.model_type)
        full_sys = build_system_data(data)
        return LirionInputs(
            args.database_path,
            args.model_type,
            args.optimizer,
            args.write_outputs,
            full_sys,
            nothing,
        )
    end

    parts_dict, shared_buses = prepare_admm_partitions(data)
    interzonal_branches = _find_interzonal_branches(data)

    zones = sort!(collect(keys(parts_dict)))
    zone_to_idx = Dict{Int, Int}(z => i for (i, z) in pairs(zones))
    parts = Vector{SystemData}(undef, length(zones))
    for (i, z) in pairs(zones)
        parts[i] = build_system_data(parts_dict[z])
    end

    admm_sys = ADMMSystem(zones, zone_to_idx, parts, shared_buses, interzonal_branches)

    return LirionInputs(
        args.database_path,
        args.model_type,
        args.optimizer,
        args.write_outputs,
        nothing,
        admm_sys,
    )
end
