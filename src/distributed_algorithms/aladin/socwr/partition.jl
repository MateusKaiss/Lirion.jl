# ============================================================
# ALADIN-SOCWR — Per-partition setup
# ============================================================

# ============================================================
# Per-partition initialisers
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function init_state_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    cfg::ALADINConfig,
    interzonal_branches::Vector{Int},
)
    build_local_model!(part, sys_zone, State, optimizer_factory, interzonal_branches)
    return build_aladin_zone_data(part; include_edges = true)
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function init_power_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    cfg::ALADINConfig,
    interzonal_branches::Vector{Int},
)
    build_local_model!(part, sys_zone, Power, optimizer_factory, interzonal_branches)
    return build_aladin_zone_data(part; include_edges = true)
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    cfg::ALADINConfig,
    ::TieLine,
    interzonal_branches::Vector{Int},
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    return build_aladin_zone_data(part; include_edges = true)
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    cfg::ALADINConfig,
    ::DummyBus,
    interzonal_branches::Vector{Int},
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    return build_aladin_zone_data(part; include_edges = true)
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    cfg::ALADINConfig,
    ::SplitBus,
    interzonal_branches::Vector{Int},
    sys_bus_pairs::Vector{NTuple{3, Int}},
    p_inj_store::Dict{Int, JuMP.VariableRef},
    q_inj_store::Dict{Int, JuMP.VariableRef},
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    _attach_pq_inj_variables_socwr!(part, sys_bus_pairs, p_inj_store, q_inj_store)
    return build_aladin_zone_data(part; include_edges = true)
end

# ============================================================
# Batch entry points
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function init_state_partition!(
    partitions::Partitions{PartitionSOCWR},
    sys::DistributedSystem,
    optimizer_factory,
    cfg::ALADINConfig,
)
    inter = sys.interzonal_branches

    n = length(partitions.parts)
    @assert length(sys.parts) == n
    zones = Vector{AladinZoneDataSOCWR}(undef, n)
    @inbounds for i in 1:n
        zones[i] = init_state_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory, cfg, inter,
        )
    end
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    interzonal_line_pairs = _build_interzonal_line_pairs(partitions, inter)
    return zones, bus_pairs, interzonal_line_pairs
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function init_power_partition!(
    partitions::Partitions{PartitionSOCWR},
    sys::DistributedSystem,
    optimizer_factory,
    cfg::ALADINConfig,
)
    inter = sys.interzonal_branches
    sys_branch_pairs = sys.branch_pairs

    n = length(partitions.parts)
    @assert length(sys.parts) == n
    zones = Vector{AladinZoneDataSOCWR}(undef, n)
    @inbounds for i in 1:n
        zones[i] = init_power_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory, cfg, inter,
        )
    end
    interzonal_line_pairs = _build_interzonal_line_pairs(partitions, inter)
    dummy_pairs = _build_dummy_line_pairs(partitions, zones, sys_branch_pairs)
    return zones, interzonal_line_pairs, dummy_pairs
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionSOCWR},
    sys::DistributedSystem,
    optimizer_factory,
    cfg::ALADINConfig,
    sectioning::TieLine,
)
    inter = sys.interzonal_branches

    n = length(partitions.parts)
    @assert length(sys.parts) == n
    zones = Vector{AladinZoneDataSOCWR}(undef, n)
    @inbounds for i in 1:n
        zones[i] = init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory, cfg,
            sectioning, inter,
        )
    end
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    interzonal_line_pairs = _build_interzonal_line_pairs(partitions, inter)
    return zones, bus_pairs, interzonal_line_pairs
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionSOCWR},
    sys::DistributedSystem,
    optimizer_factory,
    cfg::ALADINConfig,
    sectioning::DummyBus,
)
    inter = sys.interzonal_branches
    sys_branch_pairs = sys.branch_pairs

    n = length(partitions.parts)
    @assert length(sys.parts) == n
    zones = Vector{AladinZoneDataSOCWR}(undef, n)
    @inbounds for i in 1:n
        zones[i] = init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory, cfg,
            sectioning, inter,
        )
    end
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    interzonal_line_pairs = _build_interzonal_line_pairs(partitions, inter)
    dummy_pairs = _build_dummy_line_pairs(partitions, zones, sys_branch_pairs)
    return zones, bus_pairs, interzonal_line_pairs, dummy_pairs
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionSOCWR},
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
    q_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    zones = Vector{AladinZoneDataSOCWR}(undef, n_zones)
    @inbounds for i in 1:n_zones
        zones[i] = init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory, cfg,
            sectioning, inter, sys_bus_pairs, p_inj_vars[i], q_inj_vars[i],
        )
    end
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    interzonal_line_pairs = _build_interzonal_line_pairs(partitions, inter)
    split_pairs = _build_split_bus_pairs(zones, sys.bus_groups, p_inj_vars, q_inj_vars)
    return zones, bus_pairs, interzonal_line_pairs, split_pairs, p_inj_vars, q_inj_vars
end
