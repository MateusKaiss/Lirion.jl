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
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    build_local_model!(part, sys_zone, State, optimizer_factory)
    build_quad_bus_expr!(part, rho_theta, rho_v_mag)
    @objective(part.model, Min, part.obj_base + part.quad_bus_expr)
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function init_power_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    build_local_model!(part, sys_zone, Power, optimizer_factory, interzonal_branches)
    build_quad_branch_expr!(part, rho_p, rho_q, part.zone_edges, branch_pairs)
    @objective(part.model, Min, part.obj_base + part.quad_branch_expr)
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
    rho_p::Float64,
    rho_q::Float64,
    ::TieLine,
    interzonal_branches,
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    build_quad_bus_expr!(part, rho_theta, rho_v_mag)
    build_quad_branch_expr!(part, rho_p, rho_q, part.zone_edges, NTuple{4, Int}[])
    @objective(
        part.model, Min,
        part.obj_base + part.quad_bus_expr + part.quad_branch_expr,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
    rho_p::Float64,
    rho_q::Float64,
    ::DummyBus,
    interzonal_branches,
    branch_pairs::Vector{NTuple{4, Int}},
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    build_quad_bus_expr!(part, rho_theta, rho_v_mag)
    build_quad_branch_expr!(part, rho_p, rho_q, part.zone_edges, branch_pairs)
    @objective(
        part.model, Min,
        part.obj_base + part.quad_bus_expr + part.quad_branch_expr,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
    rho_p::Float64,
    rho_q::Float64,
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
    build_quad_bus_expr!(part, rho_theta, rho_v_mag)
    build_quad_p_inj_expr_ac!(part, rho_p, p_inj_store)
    build_quad_q_inj_expr_ac!(part, rho_q, q_inj_store)
    @objective(
        part.model, Min,
        part.obj_base + part.quad_bus_expr +
        part.quad_p_inj_expr + part.quad_q_inj_expr,
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

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function init_power_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_power_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_p, rho_q, interzonal_branches, branch_pairs,
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
    rho_theta::Float64,
    rho_v_mag::Float64,
    rho_p::Float64,
    rho_q::Float64,
    sectioning::TieLine,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_theta, rho_v_mag, rho_p, rho_q, sectioning, interzonal_branches,
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
    rho_theta::Float64,
    rho_v_mag::Float64,
    rho_p::Float64,
    rho_q::Float64,
    sectioning::DummyBus,
    interzonal_branches,
    branch_pairs::Vector{NTuple{4, Int}},
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_theta, rho_v_mag, rho_p, rho_q, sectioning,
            interzonal_branches, branch_pairs,
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
    rho_theta::Float64,
    rho_v_mag::Float64,
    rho_p::Float64,
    rho_q::Float64,
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
            rho_theta, rho_v_mag, rho_p, rho_q, sectioning,
            interzonal_branches, bus_pairs, p_inj_vars[i], q_inj_vars[i],
        )
    end
    return nothing
end
