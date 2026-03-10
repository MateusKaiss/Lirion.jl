function map_items_by_bus(bus_vec::AbstractVector{<:Integer})
    mapping = Dict{Int, Vector{Int}}()
    @inbounds for (i, b) in enumerate(bus_vec)
        push!(get!(mapping, Int(b), Int[]), i)
    end
    return mapping
end

function add_case_index_maps!(raw_case::CaseRaw)
    m = raw_case.maps

    empty!(m.bus_idx_map)
    if raw_case.bus !== nothing
        bus_ids = raw_case.bus.bus_i
        @inbounds for (i, b) in enumerate(bus_ids)
            m.bus_idx_map[Int(b)] = i
        end
    end

    empty!(m.bus_to_gens)
    empty!(m.gen_idx_map)
    if raw_case.gen !== nothing
        gen_bus = raw_case.gen.gen_bus
        resize!(m.gen_idx_map, length(gen_bus))
        @inbounds for (i, b) in enumerate(gen_bus)
            m.gen_idx_map[i] = Int(b)
        end
        m.bus_to_gens = map_items_by_bus(gen_bus)
    end

    empty!(m.bus_to_branches_from)
    empty!(m.bus_to_branches_to)
    empty!(m.bus_to_neighbors)
    if raw_case.branch !== nothing
        f_bus = raw_case.branch.f_bus
        t_bus = raw_case.branch.t_bus

        m.bus_to_branches_from = map_items_by_bus(f_bus)
        m.bus_to_branches_to = map_items_by_bus(t_bus)

        bus_to_neighbors = Dict{Int, Vector{Int}}()
        @inbounds for (f, t) in zip(f_bus, t_bus)
            push!(get!(bus_to_neighbors, Int(f), Int[]), Int(t))
            push!(get!(bus_to_neighbors, Int(t), Int[]), Int(f))
        end
        m.bus_to_neighbors = bus_to_neighbors
    end

    empty!(m.bus_to_prosumers)
    if raw_case.prosumer !== nothing
        m.bus_to_prosumers = map_items_by_bus(raw_case.prosumer.bus)
    end

    empty!(m.bus_to_batteries)
    if raw_case.storage !== nothing
        m.bus_to_batteries = map_items_by_bus(raw_case.storage.storage_bus)
    end

    return raw_case
end

function convert_angles_to_rad!(raw_case::CaseRaw)
    deg2rad = π / 180

    if raw_case.bus !== nothing
        va = raw_case.bus.va
        @inbounds for i in eachindex(va)
            va[i] *= deg2rad
        end
    end

    if raw_case.branch !== nothing
        br = raw_case.branch
        vshift = br.shift
        v_min = br.angmin
        v_max = br.angmax

        @inbounds for i in eachindex(vshift)
            vshift[i] *= deg2rad
        end
        @inbounds for i in eachindex(v_min)
            v_min[i] *= deg2rad
        end
        @inbounds for i in eachindex(v_max)
            v_max[i] *= deg2rad
        end
    end

    return raw_case
end

function convert_to_pu!(raw_case::CaseRaw)
    Sbase = raw_case.baseMVA

    if raw_case.bus !== nothing
        b = raw_case.bus
        b.pd ./= Sbase
        b.qd ./= Sbase
        b.gs ./= Sbase
        b.bs ./= Sbase
    end

    if raw_case.gen !== nothing
        g = raw_case.gen
        g.pg ./= Sbase
        g.qg ./= Sbase
        g.pmax ./= Sbase
        g.pmin ./= Sbase
        g.qmax ./= Sbase
        g.qmin ./= Sbase
    end

    if raw_case.branch !== nothing
        br = raw_case.branch
        br.rate_a ./= Sbase
        br.rate_b ./= Sbase
        br.rate_c ./= Sbase
    end

    if raw_case.storage !== nothing
        s = raw_case.storage
        s.ps ./= Sbase
        s.qs ./= Sbase
        s.energy ./= Sbase
        s.energy_rating ./= Sbase
        s.charge_rating ./= Sbase
        s.discharge_rating ./= Sbase
        s.thermal_rating ./= Sbase
        s.qmin ./= Sbase
        s.qmax ./= Sbase
        s.q_loss ./= Sbase
        s.p_loss ./= Sbase
    end

    return raw_case
end

function convert_timeseries_to_pu!(raw_case::CaseRaw)
    Sbase = raw_case.baseMVA

    if raw_case.bus_timeseries !== nothing
        ts = raw_case.bus_timeseries
        ts.pd ./= Sbase
        ts.qd ./= Sbase
        ts.gs ./= Sbase
        ts.bs ./= Sbase
    end

    if raw_case.prosumer_timeseries !== nothing
        ts = raw_case.prosumer_timeseries
        ts.p_load ./= Sbase
        ts.p_pv ./= Sbase
        ts.p_curtail ./= Sbase
    end

    if raw_case.consumer_timeseries !== nothing
        ts = raw_case.consumer_timeseries
        ts.demand ./= Sbase
        ts.q_demand ./= Sbase
    end

    return raw_case
end

# ---------------------------------------------------------------------
# NOTE / TODO (Timeseries review needed)
#
# This whole timeseries preprocessing pipeline still needs a full review.
# ---------------------------------------------------------------------

function preprocess_case_data!(
    raw_case::CaseRaw,
    ::Type{T} = OPFModel,
) where {T <: OPFModel}
    if is_centralized_model(T)
        all_time_steps = get_model_times(raw_case)
        time_horizon = isempty(all_time_steps) ? 1 : maximum(all_time_steps)
        raw_case.T = time_horizon
    end

    time_horizon = raw_case.T

    check_time_alignment(raw_case)

    ensure_bus_timeseries!(raw_case, time_horizon)
    ensure_gencost_timeseries!(raw_case, time_horizon)
    ensure_prosumer_timeseries!(raw_case, time_horizon)
    ensure_consumer_timeseries!(raw_case, time_horizon)

    extend_bus_timeseries!(raw_case, time_horizon)
    extend_gencost_timeseries!(raw_case, time_horizon)
    extend_prosumer_timeseries!(raw_case, time_horizon)
    extend_consumer_timeseries!(raw_case, time_horizon)

    sort_timeseries_by_time_and_id!(raw_case)

    convert_angles_to_rad!(raw_case)
    add_case_index_maps!(raw_case)

    convert_to_pu!(raw_case)
    convert_timeseries_to_pu!(raw_case)

    return nothing
end
