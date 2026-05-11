# ============================================================
# Per-partition initializers
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function init_state_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
)
    build_local_model!(part, sys_zone, State, optimizer_factory)
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function init_power_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    interzonal_branches,
)
    build_local_model!(part, sys_zone, Power, optimizer_factory, interzonal_branches)
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    ::TieLine,
    interzonal_branches,
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    ::DummyBus,
    interzonal_branches,
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    ::SplitBus,
    interzonal_branches,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_store::Dict{Int, JuMP.VariableRef},
    q_inj_store::Dict{Int, JuMP.VariableRef},
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    _attach_pq_inj_variables_ac!(
        part, sys_zone, bus_pairs, p_inj_store, q_inj_store,
    )
    return nothing
end

# ============================================================
# Batch entry points
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function init_state_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_state_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function init_power_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_power_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            interzonal_branches,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    sectioning::TieLine,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            sectioning, interzonal_branches,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    sectioning::DummyBus,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            sectioning, interzonal_branches,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    sectioning::SplitBus,
    interzonal_branches,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @assert length(p_inj_vars) == n
    @assert length(q_inj_vars) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            sectioning, interzonal_branches,
            bus_pairs, p_inj_vars[i], q_inj_vars[i],
        )
    end
    return nothing
end
