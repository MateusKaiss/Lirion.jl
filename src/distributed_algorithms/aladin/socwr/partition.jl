function _build_shared_bus_pairs_socwr(
    partitions::Partitions{PartitionSOCWR},
    zones::Vector{AladinZoneDataSOCWR},
)::Vector{SharedBusPairSOCWR}
    acc = Dict{Int, Vector{Tuple{Int, Int}}}()  # cidx -> [(zpos, v2_col)]
    bus_of = Dict{Int, Int}()

    for (zpos, part) in enumerate(partitions.parts)
        zdata = zones[zpos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            cidx = part.shared_bus_consensus_idx[i]
            cidx == 0 && continue
            col = zdata.shared_v2_cols[i]
            v = get!(acc, cidx, Tuple{Int, Int}[])
            push!(v, (zpos, col))
            bus_of[cidx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPairSOCWR[]
    for cidx in sort!(collect(keys(acc)))
        entries = acc[cidx]
        length(entries) == 2 || continue
        zi, ci = entries[1]
        zj, cj = entries[2]
        push!(pairs, SharedBusPairSOCWR(cidx, bus_of[cidx], zi, ci, zj, cj))
    end
    return pairs
end

function _build_shared_edge_pairs_socwr(
    partitions::Partitions{PartitionSOCWR},
    interzonal::Vector{Int},
)::Vector{SharedEdgePairSOCWR}
    pairs = SharedEdgePairSOCWR[]
    for (cidx, e) in enumerate(interzonal)
        entries = Tuple{Int, Int}[]
        for (zpos, part) in enumerate(partitions.parts)
            part.branch_idx_full_to_zone[e] == 0 && continue
            local_k = findfirst(==(e), part.zone_edges)
            local_k === nothing && continue
            push!(entries, (zpos, local_k))
        end
        length(entries) == 2 || continue
        zi, li = entries[1]
        zj, lj = entries[2]
        push!(pairs, SharedEdgePairSOCWR(cidx, e, zi, li, zj, lj))
    end
    return pairs
end

function _build_zone_data_vector_socwr(
    partitions::Partitions{PartitionSOCWR};
    include_edges::Bool,
)::Vector{AladinZoneDataSOCWR}
    parts = partitions.parts
    zones = Vector{AladinZoneDataSOCWR}(undef, length(parts))
    @inbounds for i in eachindex(parts)
        zones[i] = build_aladin_zone_data_socwr(parts[i]; include_edges = include_edges)
    end
    return zones
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{State, SOCWR},
)
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    build_local_model!(partitions, sys, State, optimizer_factory, inter)

    zones = _build_zone_data_vector_socwr(partitions; include_edges = true)
    bus_pairs = _build_shared_bus_pairs_socwr(partitions, zones)
    edge_pairs = _build_shared_edge_pairs_socwr(partitions, inter)

    n_bus = length(bus_pairs)
    n_edge = length(edge_pairs)

    return SOCWRStateWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        edge_pairs = edge_pairs,
        vars = (
            v2 = LambdaState(n_bus),
            r = LambdaState(n_edge),
            t = LambdaState(n_edge),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
    )
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Power, SOCWR},
)
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    build_local_model!(partitions, sys, Power, optimizer_factory, inter)

    zones = _build_zone_data_vector_socwr(partitions; include_edges = true)
    edge_pairs = _build_shared_edge_pairs_socwr(partitions, inter)
    n_edge = length(edge_pairs)

    return SOCWRPowerWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        edge_pairs = edge_pairs,
        vars = (
            p_flow = LambdaState(n_edge),
            q_flow = LambdaState(n_edge),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
    )
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Hybrid, SOCWR},
)
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    build_local_model!(partitions, sys, Hybrid, optimizer_factory, inter)

    zones = _build_zone_data_vector_socwr(partitions; include_edges = true)
    bus_pairs = _build_shared_bus_pairs_socwr(partitions, zones)
    edge_pairs = _build_shared_edge_pairs_socwr(partitions, inter)

    n_bus = length(bus_pairs)
    n_edge = length(edge_pairs)

    return SOCWRHybridWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        edge_pairs = edge_pairs,
        vars = (
            v2 = LambdaState(n_bus),
            p_flow = LambdaState(n_edge),
            q_flow = LambdaState(n_edge),
            r = LambdaState(n_edge),
            t = LambdaState(n_edge),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
    )
end
