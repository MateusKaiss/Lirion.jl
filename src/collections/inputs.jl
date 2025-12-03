struct ADMMSystem
    full::SystemData     
    parts::Dict{Int, SystemData}  
    shared_buses::Vector{Int}         
    cut_edges::Vector{Int} 
end

struct LirionInputs
    database_path::String
    model_type::Union{Type{<:OPFModel}, DistributedOPFModel}
    optimizer::Any
    write_outputs::Bool
    system::SystemData
    admm_system::Union{Nothing, ADMMSystem}
end

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

function _find_cut_edges(full::Dict{String, Any})
    zones    = full["bus"]["zone"]
    bus_ids  = full["bus"]["bus_i"]
    bus2zone = Dict(bus_ids[i] => zones[i] for i in eachindex(bus_ids))

    f_bus = full["branch"]["f_bus"]
    t_bus = full["branch"]["t_bus"]
    br_status =
        haskey(full["branch"], "br_status") ? full["branch"]["br_status"] :
        ones(Int, length(f_bus))

    cut = Int[]
    for i in eachindex(f_bus)
        if br_status[i] == 1 && bus2zone[f_bus[i]] != bus2zone[t_bus[i]]
            push!(cut, i)
        end
    end
    return cut
end

function prepare_admm_partitions(data::Dict{String, Any})
    bus      = data["bus"]
    bus_ids  = bus["bus_i"]
    zones    = bus["zone"]
    unique_zones = unique(zones)

    partitions = Dict{Int, Dict{String, Any}}()

    has_branch = haskey(data, "branch")
    if has_branch
        f_bus_full = data["branch"]["f_bus"]
        t_bus_full = data["branch"]["t_bus"]
        br_status  =
            haskey(data["branch"], "br_status") ?
            data["branch"]["br_status"] :
            ones(Int, length(f_bus_full))
        @assert length(br_status) == length(f_bus_full)
        active_mask_full = [br_status[i] == 1 for i in eachindex(f_bus_full)]
    end

    for zone in unique_zones
        zone_mask    = findall(==(zone), zones)
        zone_bus_ids = @views bus_ids[zone_mask]
        zone_bus_set = Set(zone_bus_ids)

        extra = Set{Int}()
        if has_branch
            @inbounds for i in eachindex(f_bus_full)
                f_in = in(f_bus_full[i], zone_bus_set)
                t_in = in(t_bus_full[i], zone_bus_set)
                if f_in ⊻ t_in
                    push!(extra, f_in ? t_bus_full[i] : f_bus_full[i])
                end
            end
        end

        all_bus_set = union(zone_bus_set, extra)
        bus_mask    = findall(b -> in(b, all_bus_set), bus_ids)

        sub = Dict{String, Any}()
        sub["baseMVA"]        = data["baseMVA"]
        sub["T"]              = get(data, "T", 1)
        sub["source_type"]    = get(data, "source_type", "matpower")
        sub["source_version"] = get(data, "source_version", "")
        sub["name"]           = "$(get(data,"name","case"))_zone_$zone"

        sub["bus"] = Dict(
            k => v[bus_mask] for (k, v) in data["bus"] if length(v) == length(bus_ids)
        )

        zone_bus_set = Set(zone_bus_ids)
        is_virtual = [!in(b, zone_bus_set) for b in sub["bus"]["bus_i"]]
        sub["bus"]["is_virtual"] = is_virtual
        @inbounds for (i, b) in enumerate(sub["bus"]["bus_i"])
            if is_virtual[i]
                sub["bus"]["pd"][i]       = 0.0
                sub["bus"]["qd"][i]       = 0.0
                sub["bus"]["bus_type"][i] = 1
            end
        end

        for key in ("gen", "prosumer", "consumer", "storage")
            if haskey(data, key)
                tbl = data[key]
                bus_key =
                    haskey(tbl, "bus")         ? "bus"         :
                    haskey(tbl, "gen_bus")     ? "gen_bus"     :
                    haskey(tbl, "storage_bus") ? "storage_bus" : nothing
                if bus_key !== nothing
                    keep = findall(
                        i -> in(tbl[bus_key][i], zone_bus_set),
                        eachindex(tbl[bus_key]),
                    )
                    sub[key] =
                        !isempty(keep) ?
                        Dict(
                            k => v[keep] for
                            (k, v) in tbl if length(v) == length(tbl[bus_key])
                        ) :
                        Dict(k => Any[] for k in keys(tbl))
                else
                    sub[key] = Dict(k => Any[] for k in keys(tbl))
                end
            end
        end

        if has_branch
            keep_all = findall(
                i ->
                    (active_mask_full[i]) &&
                    (in(f_bus_full[i], zone_bus_set) || in(t_bus_full[i], zone_bus_set)),
                eachindex(f_bus_full),
            )

            sub["branch"] = Dict(
                k => v[keep_all] for
                (k, v) in data["branch"] if length(v) == length(f_bus_full)
            )

            sub["branch_full"]          = data["branch"]
            sub["branch_local_indices"] = keep_all

            edge_to_local = Dict{Int, Int}()
            @inbounds for (loc, glob) in enumerate(keep_all)
                edge_to_local[glob] = loc
            end
            sub["edge_to_local"] = edge_to_local
        end

        for key in ("bus_timeseries", "prosumer_timeseries", "consumer_timeseries")
            if haskey(data, key)
                ts   = data[key]
                keep = findall(i -> in(ts["bus"][i], zone_bus_set), eachindex(ts["bus"]))
                sub[key] =
                    Dict(k => v[keep] for (k, v) in ts if length(v) == length(ts["bus"]))
            end
        end

        if haskey(data, "gen")
            if haskey(data, "gencost")
                gtbl = data["gen"]
                bkey = haskey(gtbl, "bus") ? "bus" : "gen_bus"
                keep = findall(i -> in(gtbl[bkey][i], zone_bus_set), eachindex(gtbl[bkey]))
                sub["gencost"] = Dict(
                    k => v[keep] for
                    (k, v) in data["gencost"] if length(v) == length(gtbl[bkey])
                )
            end
            if haskey(data, "gencost_timeseries")
                gtbl = data["gen"]
                bkey = haskey(gtbl, "bus") ? "bus" : "gen_bus"
                keep = findall(i -> in(gtbl[bkey][i], zone_bus_set), eachindex(gtbl[bkey]))
                sub["gencost_timeseries"] = Dict(
                    k => v[keep] for
                    (k, v) in data["gencost_timeseries"] if length(v) == length(gtbl[bkey])
                )
            end
        else
            sub["gen"]                = Dict{String, Vector{Any}}()
            sub["gencost"]            = Dict{String, Vector{Any}}()
            sub["gencost_timeseries"] = Dict{String, Vector{Any}}()
        end

        add_case_index_maps!(sub)

        sub["bus_set"] = Set(sub["bus"]["bus_i"])
        partitions[zone] = sub
    end

    bus_zone_appearance = Dict{Int, Set{Int}}()
    for (zone_id, sub) in partitions, b in sub["bus"]["bus_i"]
        push!(get!(bus_zone_appearance, b, Set{Int}()), zone_id)
    end
    shared_buses = [b for (b, zs) in bus_zone_appearance if length(zs) > 1]
    shared_set   = Set(shared_buses)

    for (_, sub) in partitions
        local_b = sub["bus"]["bus_i"]
        sub["shared_buses_local"] = [b for b in local_b if in(b, shared_set)]
    end

    return partitions, shared_buses
end

function build_inputs(args::Args)::LirionInputs
    # inputs
    data = load_input_data(args.database_path)
    preprocess_case_data!(data, OPFModel)

    # structs
    full_sys = build_system_data(data)

    if is_centralized_model(args.model_type)
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

    admm_sys = ADMMSystem(
        full_sys,
        parts_system,
        shared_buses,
        cut_edges,
    )

    return LirionInputs(
        args.database_path,
        args.model_type,
        args.optimizer,
        args.write_outputs,
        full_sys,
        admm_sys,
    )
end