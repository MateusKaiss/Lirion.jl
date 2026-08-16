# ============================================================
# Helpers
# ============================================================

@inline function _zone_cut_edges(part, interzonal_branches::AbstractVector{Int})
    m = part.branch_idx_full_to_zone
    return [e for e in interzonal_branches if m[e] != 0]
end

@inline function _setup_partition_model!(
    part::AbstractPartition,
    sys_zone::SystemData,
    ::Type{M},
    optimizer_factory,
) where {M <: OPFModel}
    prepare_hotpaths!(part)
    model = build_model(sys_zone, M; optimizer = optimizer_factory)
    part.local_model = model
    part.obj_base = model[:generation_cost]
    return model
end

# ============================================================
# DC
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function build_local_model!(
    part::PartitionDC,
    sys_zone::SystemData,
    ::Type{State},
    optimizer_factory,
)
    model = _setup_partition_model!(part, sys_zone, DC, optimizer_factory)
    part.theta_var = model[:theta]
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_local_model!(
    part::PartitionDC,
    sys_zone::SystemData,
    ::Type{Power},
    optimizer_factory,
    interzonal_branches::Vector{Int},
)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)
    model = _setup_partition_model!(part, sys_zone, DC, optimizer_factory)
    part.theta_var = model[:theta]
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    return nothing
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function build_local_model!(
    part::PartitionDC,
    sys_zone::SystemData,
    ::Type{Hybrid},
    optimizer_factory,
    interzonal_branches::Vector{Int},
)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)
    model = _setup_partition_model!(part, sys_zone, DC, optimizer_factory)
    part.theta_var = model[:theta]
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    return nothing
end

# ============================================================
# AC
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function build_local_model!(
    part::PartitionAC,
    sys_zone::SystemData,
    ::Type{State},
    optimizer_factory,
)
    model = _setup_partition_model!(part, sys_zone, AC, optimizer_factory)
    part.theta_var = model[:theta]
    part.v_mag_var = model[:v_magnitude]
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    for I in keys(part.v_mag_var)
        set_start_value(part.v_mag_var[I], 1.0)
    end
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_local_model!(
    part::PartitionAC,
    sys_zone::SystemData,
    ::Type{Power},
    optimizer_factory,
    interzonal_branches,
)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)
    model = _setup_partition_model!(part, sys_zone, AC, optimizer_factory)
    part.theta_var = model[:theta]
    part.v_mag_var = model[:v_magnitude]
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    for I in keys(part.v_mag_var)
        set_start_value(part.v_mag_var[I], 1.0)
    end
    return nothing
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function build_local_model!(
    part::PartitionAC,
    sys_zone::SystemData,
    ::Type{Hybrid},
    optimizer_factory,
    interzonal_branches,
)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)
    model = _setup_partition_model!(part, sys_zone, AC, optimizer_factory)
    part.theta_var = model[:theta]
    part.v_mag_var = model[:v_magnitude]
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    for I in keys(part.v_mag_var)
        set_start_value(part.v_mag_var[I], 1.0)
    end
    return nothing
end

# ============================================================
# SOCWR
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function build_local_model!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    ::Type{State},
    optimizer_factory,
    interzonal_branches,
)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)
    model = _setup_partition_model!(part, sys_zone, SOCWR, optimizer_factory)
    part.v2_var = model[:v_magnitude_sq]
    for I in keys(part.v2_var)
        set_start_value(part.v2_var[I], 1.0)
    end
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_local_model!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    ::Type{Power},
    optimizer_factory,
    interzonal_branches,
)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)
    _setup_partition_model!(part, sys_zone, SOCWR, optimizer_factory)
    return nothing
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function build_local_model!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    ::Type{Hybrid},
    optimizer_factory,
    interzonal_branches,
)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)
    model = _setup_partition_model!(part, sys_zone, SOCWR, optimizer_factory)
    part.v2_var = model[:v_magnitude_sq]
    for I in keys(part.v2_var)
        set_start_value(part.v2_var[I], 1.0)
    end
    return nothing
end

# ============================================================
# Batch entry point
# ============================================================

function build_local_model!(
    partitions::Partitions{P},
    sys::DistributedSystem,
    ::Type{C},
    optimizer_factory,
    args...,
) where {P <: AbstractPartition, C <: ADMMVariant}
    n = length(partitions.parts)
    @assert length(sys.parts) == n "DistributedSystem.parts must be aligned with partitions.parts"
    @inbounds for i in 1:n
        build_local_model!(partitions.parts[i], sys.parts[i], C, optimizer_factory, args...)
    end
    return nothing
end
