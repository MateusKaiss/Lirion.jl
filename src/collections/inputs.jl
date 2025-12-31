# ============================================================
# Core structs
# ============================================================

struct ADMMSystem
    parts::Dict{Int, SystemData}
    shared_buses::Vector{Int}
    cut_edges::Vector{Int}
end

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
    cut_edges = _find_cut_edges(data)

    parts_system = Dict{Int, SystemData}()
    for (zone, subdict) in parts_dict
        parts_system[zone] = build_system_data(subdict)
    end

    admm_sys = ADMMSystem(parts_system, shared_buses, cut_edges)

    return LirionInputs(
        args.database_path,
        args.model_type,
        args.optimizer,
        args.write_outputs,
        nothing,
        admm_sys,
    )
end
