# ============================================================
# DC PARTITION INITIALIZATION
# ============================================================

function init_state_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
)
    model = _setup_partition_model!(part, sys_zone, DC, optimizer_factory)

    part.theta_var = model[:theta]
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end

    build_quad_bus_expr!(part, rho_theta)
    @objective(model, Min, part.obj_base + part.quad_bus_expr)

    return nothing
end

function init_power_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)
    model = _setup_partition_model!(part, sys_zone, DC, optimizer_factory)

    part.theta_var = model[:theta]
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end

    build_quad_branch_expr!(part, rho, part.zone_edges)
    @objective(model, Min, part.obj_base + part.quad_branch_expr)

    return nothing
end

function init_hybrid_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)
    model = _setup_partition_model!(part, sys_zone, DC, optimizer_factory)

    part.theta_var = model[:theta]
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end

    build_quad_bus_expr!(part, rho_theta)
    build_quad_branch_expr!(part, rho, part.zone_edges)

    @objective(model, Min, part.obj_base + part.quad_bus_expr + part.quad_branch_expr)

    return nothing
end

# ============================================================
# BATCH
# ============================================================

function init_state_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n "DistributedSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n
        init_state_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_theta,
        )
    end
    return nothing
end

function init_power_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n "DistributedSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n
        init_power_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho,
            interzonal_branches,
        )
    end
    return nothing
end

function init_hybrid_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n "DistributedSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_theta,
            rho,
            interzonal_branches,
        )
    end
    return nothing
end
