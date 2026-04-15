function _build_shared_bus_pairs(
    partitions::Partitions{PartitionDC},
    zones::Vector{AladinZoneData},
)::Vector{SharedBusPair}
    acc = Dict{Int, Vector{Tuple{Int, Int}}}()  # consensus_idx -> [(zone_pos, col)]
    bus_of = Dict{Int, Int}()                   # consensus_idx -> bus_id

    for (zpos, part) in enumerate(partitions.parts)
        zdata = zones[zpos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            cidx = part.shared_bus_consensus_idx[i]
            cidx == 0 && continue
            col = zdata.shared_theta_cols[i]
            v = get!(acc, cidx, Tuple{Int, Int}[])
            push!(v, (zpos, col))
            bus_of[cidx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPair[]
    for cidx in sort!(collect(keys(acc)))
        entries = acc[cidx]
        length(entries) == 2 || continue
        zi, ci = entries[1]
        zj, cj = entries[2]
        push!(pairs, SharedBusPair(cidx, bus_of[cidx], zi, ci, zj, cj))
    end
    return pairs
end

function _build_shared_edge_pairs(
    partitions::Partitions{PartitionDC},
    interzonal::Vector{Int},
)::Vector{SharedEdgePair}
    pairs = SharedEdgePair[]

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
        push!(pairs, SharedEdgePair(cidx, e, zi, li, zj, lj))
    end
    return pairs
end

function _build_zone_data_vector(
    partitions::Partitions{PartitionDC};
    include_edges::Bool,
)::Vector{AladinZoneData}
    parts = partitions.parts
    zones = Vector{AladinZoneData}(undef, length(parts))
    @inbounds for i in eachindex(parts)
        zones[i] = build_aladin_zone_data(parts[i]; include_edges = include_edges)
    end
    return zones
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{State, DC},
)
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones

    build_local_model!(partitions, sys, State, optimizer_factory)

    zones = _build_zone_data_vector(partitions; include_edges = false)
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    n_pairs = length(bus_pairs)

    return DCStateWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        zones = zones,
        bus_pairs = bus_pairs,
        vars = (theta = LambdaState(n_pairs),),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_theta_other = cfg.sigma_theta_other,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
    )
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Power, DC},
)
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    build_local_model!(partitions, sys, Power, optimizer_factory, inter)

    zones = _build_zone_data_vector(partitions; include_edges = true)
    edge_pairs = _build_shared_edge_pairs(partitions, inter)
    n_pairs = length(edge_pairs)

    return DCPowerWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        edge_pairs = edge_pairs,
        vars = (p_flow = LambdaState(n_pairs),),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_theta_other = cfg.sigma_theta_other,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
    )
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Hybrid, DC},
)
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    build_local_model!(partitions, sys, Hybrid, optimizer_factory, inter)

    zones = _build_zone_data_vector(partitions; include_edges = true)
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    edge_pairs = _build_shared_edge_pairs(partitions, inter)
    n_bus = length(bus_pairs)
    n_edge = length(edge_pairs)

    return DCHybridWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        edge_pairs = edge_pairs,
        vars = (
            theta = LambdaState(n_bus),
            p_flow = LambdaState(n_edge),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho_bus = Ref(cfg.rho),
        rho_edge = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_theta_other = cfg.sigma_theta_other,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
    )
end
