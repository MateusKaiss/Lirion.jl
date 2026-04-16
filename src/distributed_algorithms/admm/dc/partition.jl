# ============================================================
# DC PARTITION INITIALIZATION — ADMM
# ------------------------------------------------------------

function init_state_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
)
    build_local_model!(part, sys_zone, State, optimizer_factory)
    build_quad_bus_expr!(part, rho_theta)
    part.obj_static = part.obj_base + part.quad_bus_expr
    @objective(part.model, Min, part.obj_static)
    return nothing
end

function init_power_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho::Float64,
    interzonal_branches::Vector{Int},
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    build_local_model!(part, sys_zone, Power, optimizer_factory, interzonal_branches)
    build_quad_branch_expr!(part, rho, part.zone_edges, branch_pairs)
    part.obj_static = part.obj_base + part.quad_branch_expr
    @objective(part.model, Min, part.obj_static)
    return nothing
end

function init_hybrid_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    interzonal_branches::Vector{Int},
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    build_quad_bus_expr!(part, rho_theta)
    build_quad_branch_expr!(part, rho, part.zone_edges, branch_pairs)
    part.obj_static = part.obj_base + part.quad_bus_expr + part.quad_branch_expr
    @objective(part.model, Min, part.obj_static)
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
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_state_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory, rho_theta,
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
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_power_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho, interzonal_branches, branch_pairs,
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
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_theta, rho, interzonal_branches, branch_pairs,
        )
    end
    return nothing
end
