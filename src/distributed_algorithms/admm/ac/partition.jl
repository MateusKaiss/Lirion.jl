# ============================================================
# AC PARTITION INITIALIZATION — ADMM
# ------------------------------------------------------------
# See admm/dc/partition.jl header.
# ============================================================

function init_state_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    build_local_model!(part, sys_zone, State, optimizer_factory)
    build_quad_bus_expr!(part, rho_theta, rho_v_mag)
    @objective(part.model, Min, part.obj_base + part.quad_bus_expr)
    return nothing
end

function init_power_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
)
    build_local_model!(part, sys_zone, Power, optimizer_factory, interzonal_branches)
    build_quad_branch_expr!(part, rho_p, rho_q, part.zone_edges)
    @objective(part.model, Min, part.obj_base + part.quad_branch_expr)
    return nothing
end

function init_hybrid_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    build_quad_bus_expr!(part, rho_theta, rho_v_mag)
    build_quad_branch_expr!(part, rho_p, rho_q, part.zone_edges)
    @objective(
        part.model,
        Min,
        part.obj_base + part.quad_bus_expr + part.quad_branch_expr,
    )
    return nothing
end

# ============================================================
# BATCH
# ============================================================

function init_state_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_state_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_theta, rho_v_mag,
        )
    end
    return nothing
end

function init_power_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_power_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_p, rho_q, interzonal_branches,
        )
    end
    return nothing
end

function init_hybrid_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_theta, rho_v_mag, rho_p, rho_q, interzonal_branches,
        )
    end
    return nothing
end
