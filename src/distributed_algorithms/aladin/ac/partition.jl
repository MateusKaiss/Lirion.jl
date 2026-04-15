function _build_shared_bus_pairs_ac(
    partitions::Partitions{PartitionAC},
    zones::Vector{AladinZoneDataAC},
)::Vector{SharedBusPairAC}
    acc = Dict{Int, Vector{Tuple{Int, Int, Int}}}()  # cidx → [(zpos, θ_col, vm_col)]
    bus_of = Dict{Int, Int}()

    for (zpos, part) in enumerate(partitions.parts)
        zdata = zones[zpos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            cidx = part.shared_bus_consensus_idx[i]
            cidx == 0 && continue
            tcol = zdata.shared_theta_cols[i]
            vcol = zdata.shared_vm_cols[i]
            v = get!(acc, cidx, Tuple{Int, Int, Int}[])
            push!(v, (zpos, tcol, vcol))
            bus_of[cidx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPairAC[]
    for cidx in sort!(collect(keys(acc)))
        entries = acc[cidx]
        length(entries) == 2 || continue
        zi, ti, vi = entries[1]
        zj, tj, vj = entries[2]
        push!(pairs, SharedBusPairAC(cidx, bus_of[cidx], zi, ti, vi, zj, tj, vj))
    end
    return pairs
end

function _build_shared_edge_pairs_ac(
    partitions::Partitions{PartitionAC},
    interzonal::Vector{Int},
)::Vector{SharedEdgePairAC}
    pairs = SharedEdgePairAC[]
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
        push!(pairs, SharedEdgePairAC(cidx, e, zi, li, zj, lj))
    end
    return pairs
end

function _build_zone_data_vector_ac(
    partitions::Partitions{PartitionAC};
    include_edges::Bool,
)::Vector{AladinZoneDataAC}
    parts = partitions.parts
    zones = Vector{AladinZoneDataAC}(undef, length(parts))
    @inbounds for i in eachindex(parts)
        zones[i] = build_aladin_zone_data_ac(parts[i]; include_edges = include_edges)
    end
    return zones
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{State, AC},
)
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    build_local_model!(partitions, sys, State, optimizer_factory)

    zones = _build_zone_data_vector_ac(partitions; include_edges = false)
    bus_pairs = _build_shared_bus_pairs_ac(partitions, zones)
    n_bus = length(bus_pairs)

    return ACStateWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        vars = (
            theta = LambdaState(n_bus),
            v_mag = LambdaState(n_bus),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
    )
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Power, AC},
)
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    build_local_model!(partitions, sys, Power, optimizer_factory, inter)

    zones = _build_zone_data_vector_ac(partitions; include_edges = true)
    edge_pairs = _build_shared_edge_pairs_ac(partitions, inter)
    n_edge = length(edge_pairs)

    return ACPowerWorkspaceALADIN(;
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
    )
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Hybrid, AC},
)
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    build_local_model!(partitions, sys, Hybrid, optimizer_factory, inter)

    zones = _build_zone_data_vector_ac(partitions; include_edges = true)
    bus_pairs = _build_shared_bus_pairs_ac(partitions, zones)
    edge_pairs = _build_shared_edge_pairs_ac(partitions, inter)
    n_bus = length(bus_pairs)
    n_edge = length(edge_pairs)

    return ACHybridWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        edge_pairs = edge_pairs,
        vars = (
            theta = LambdaState(n_bus),
            v_mag = LambdaState(n_bus),
            p_flow = LambdaState(n_edge),
            q_flow = LambdaState(n_edge),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
    )
end
