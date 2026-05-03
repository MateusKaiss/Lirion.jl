# ============================================================
# Per-partition initializers
# ============================================================
#
# `_attach_p_inj_variables!` lives in
# `distributed_algorithms/shared/sectioning_helpers.jl` (shared with
# ALADIN — pure JuMP plumbing, algorithm-agnostic).

# ------------------------------------------------------------
# State
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

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    ::TieLine,
    interzonal_branches::Vector{Int},
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    build_quad_bus_expr!(part, rho_theta)
    build_quad_branch_expr!(part, rho, part.zone_edges, NTuple{4, Int}[])
    part.obj_static =
        part.obj_base + part.quad_bus_expr + part.quad_branch_expr
    @objective(part.model, Min, part.obj_static)
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    ::DummyBus,
    interzonal_branches::Vector{Int},
    branch_pairs::Vector{NTuple{4, Int}},
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    build_quad_bus_expr!(part, rho_theta)
    build_quad_branch_expr!(part, rho, part.zone_edges, branch_pairs)
    part.obj_static =
        part.obj_base + part.quad_bus_expr + part.quad_branch_expr
    @objective(part.model, Min, part.obj_static)
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    ::SplitBus,
    interzonal_branches::Vector{Int},
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_store::Dict{Int, JuMP.VariableRef},
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    _attach_p_inj_variables!(part, bus_pairs, p_inj_store)
    build_quad_bus_expr!(part, rho_theta)
    build_quad_p_inj_expr!(part, rho, p_inj_store)
    part.obj_static =
        part.obj_base + part.quad_bus_expr + part.quad_p_inj_expr
    @objective(part.model, Min, part.obj_static)
    return nothing
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

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    sectioning::TieLine,
    interzonal_branches::Vector{Int},
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_theta, rho, sectioning, interzonal_branches,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    sectioning::DummyBus,
    interzonal_branches::Vector{Int},
    branch_pairs::Vector{NTuple{4, Int}},
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_theta, rho, sectioning, interzonal_branches, branch_pairs,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    sectioning::SplitBus,
    interzonal_branches::Vector{Int},
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @assert length(p_inj_vars) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_theta, rho, sectioning, interzonal_branches,
            bus_pairs, p_inj_vars[i],
        )
    end
    return nothing
end
