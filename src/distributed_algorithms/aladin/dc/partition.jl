# ============================================================
# ALADIN-DC — Per-partition setup
# ============================================================

# ============================================================
# Per-partition initialisers
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function init_state_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    cfg::ALADINConfig,
)
    build_local_model!(part, sys_zone, State, optimizer_factory)
    zone = build_aladin_zone_data(part; include_edges = false)
    _build_quad_static_primal_dc!(
        zone, cfg.rho;
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
    )
    return zone
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function init_power_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    cfg::ALADINConfig,
    interzonal_branches::Vector{Int},
)
    build_local_model!(part, sys_zone, Power, optimizer_factory, interzonal_branches)
    zone = build_aladin_zone_data(part; include_edges = true)
    if cfg.proximal_mode === :pflow
        _build_quad_static_pflow_dc!(
            zone, cfg.rho;
            sigma_pflow = cfg.sigma_pflow,
            sigma_pg = cfg.sigma_pg,
            sigma_other = cfg.sigma_other,
        )
    else
        _build_quad_static_primal_dc!(
            zone, cfg.rho;
            sigma_theta_shared = cfg.sigma_theta_shared,
            sigma_pg = cfg.sigma_pg,
            sigma_other = cfg.sigma_other,
        )
    end
    return zone
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    cfg::ALADINConfig,
    ::TieLine,
    interzonal_branches::Vector{Int},
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    zone = build_aladin_zone_data(part; include_edges = true)
    _build_quad_static_primal_dc!(
        zone, cfg.rho;
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
    )
    return zone
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    cfg::ALADINConfig,
    ::DummyBus,
    interzonal_branches::Vector{Int},
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    zone = build_aladin_zone_data(part; include_edges = true)
    _build_quad_static_primal_dc!(
        zone, cfg.rho;
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
    )
    return zone
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    cfg::ALADINConfig,
    ::SplitBus,
    interzonal_branches::Vector{Int},
    sys_bus_pairs::Vector{NTuple{3, Int}},
    p_inj_store::Dict{Int, JuMP.VariableRef},
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    _attach_p_inj_variables!(part, sys_bus_pairs, p_inj_store)
    zone = build_aladin_zone_data(part; include_edges = true)
    _build_quad_static_primal_dc!(
        zone, cfg.rho;
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
    )
    return zone
end

# ============================================================
# Batch entry points
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function init_state_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    cfg::ALADINConfig,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    zones = Vector{AladinZoneData}(undef, n)
    @inbounds for i in 1:n
        zones[i] = init_state_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory, cfg,
        )
    end
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    return zones, bus_pairs
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function init_power_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    cfg::ALADINConfig,
)
    inter = sys.interzonal_branches
    sys_branch_pairs = sys.branch_pairs

    n = length(partitions.parts)
    @assert length(sys.parts) == n
    zones = Vector{AladinZoneData}(undef, n)
    @inbounds for i in 1:n
        zones[i] = init_power_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory, cfg, inter,
        )
    end
    edge_pairs = _build_shared_edge_pairs(partitions, inter)
    dummy_pairs = _build_dummy_branch_pairs(partitions, zones, sys_branch_pairs)
    return zones, edge_pairs, dummy_pairs
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    cfg::ALADINConfig,
    sectioning::TieLine,
)
    inter = sys.interzonal_branches

    n = length(partitions.parts)
    @assert length(sys.parts) == n
    zones = Vector{AladinZoneData}(undef, n)
    @inbounds for i in 1:n
        zones[i] = init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory, cfg,
            sectioning, inter,
        )
    end
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    edge_pairs = _build_shared_edge_pairs(partitions, inter)
    return zones, bus_pairs, edge_pairs
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    cfg::ALADINConfig,
    sectioning::DummyBus,
)
    inter = sys.interzonal_branches
    sys_branch_pairs = sys.branch_pairs

    n = length(partitions.parts)
    @assert length(sys.parts) == n
    zones = Vector{AladinZoneData}(undef, n)
    @inbounds for i in 1:n
        zones[i] = init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory, cfg,
            sectioning, inter,
        )
    end
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    edge_pairs = _build_shared_edge_pairs(partitions, inter)
    dummy_pairs = _build_dummy_branch_pairs(partitions, zones, sys_branch_pairs)
    return zones, bus_pairs, edge_pairs, dummy_pairs
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    cfg::ALADINConfig,
    sectioning::SplitBus,
)
    inter = sys.interzonal_branches
    sys_bus_pairs = sys.bus_pairs

    n_zones = length(partitions.parts)
    @assert length(sys.parts) == n_zones
    p_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    zones = Vector{AladinZoneData}(undef, n_zones)
    @inbounds for i in 1:n_zones
        zones[i] = init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory, cfg,
            sectioning, inter, sys_bus_pairs, p_inj_vars[i],
        )
    end
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    edge_pairs = _build_shared_edge_pairs(partitions, inter)
    split_pairs = _build_split_bus_pairs(partitions, zones, sys_bus_pairs, p_inj_vars)
    return zones, bus_pairs, edge_pairs, split_pairs, p_inj_vars
end
