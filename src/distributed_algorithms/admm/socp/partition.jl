# ============================================================
# SOCP PARTITION INITIALIZATION
# ============================================================

function init_state_partition!(
    part::PartitionSOCP,
    sys_zone::SystemData,
    optimizer_factory,
    rho_v2::Float64,
    rho_r::Float64,
    rho_t::Float64,
    interzonal_branches,
)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)
    model = _setup_partition_model!(part, sys_zone, SOCP, optimizer_factory)

    part.v2_var = model[:v_magnitude_sq]
    for I in keys(part.v2_var)
        set_start_value(part.v2_var[I], 1.0)
    end

    build_quad_bus_expr!(part, rho_v2)
    build_quad_branch_rt_expr!(part, rho_r, rho_t, part.zone_edges)

    part.obj_static = part.obj_base + part.quad_bus_expr_socp + part.quad_rt_expr_socp
    @objective(model, Min, part.obj_static)

    return nothing
end

function init_power_partition!(
    part::PartitionSOCP,
    sys_zone::SystemData,
    optimizer_factory,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)
    model = _setup_partition_model!(part, sys_zone, SOCP, optimizer_factory)

    if is_gurobi(model)
        soft_anchor_branch_expr!(
            part;
            kappa_v2 = 1e-6,
            kappa_pq = 1e-8,
            use_bus = true,
            use_pq = true,
        )
    else
        part.soft_anchor_expr = zero(JuMP.QuadExpr)
    end

    build_quad_branch_pq_expr!(part, rho_p, rho_q, part.zone_edges)

    part.obj_static = part.obj_base + part.quad_pq_expr_socp + part.soft_anchor_expr
    @objective(model, Min, part.obj_static)

    return nothing
end

function init_hybrid_partition!(
    part::PartitionSOCP,
    sys_zone::SystemData,
    optimizer_factory,
    rho_v2::Float64,
    rho_p::Float64,
    rho_q::Float64,
    rho_r::Float64,
    rho_t::Float64,
    interzonal_branches,
)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)
    model = _setup_partition_model!(part, sys_zone, SOCP, optimizer_factory)

    part.v2_var = model[:v_magnitude_sq]
    for I in keys(part.v2_var)
        set_start_value(part.v2_var[I], 1.0)
    end

    build_quad_bus_expr!(part, rho_v2)
    build_quad_branch_pq_expr!(part, rho_p, rho_q, part.zone_edges)
    build_quad_branch_rt_expr!(part, rho_r, rho_t, part.zone_edges)

    part.obj_static =
        part.obj_base + part.quad_bus_expr_socp + part.quad_pq_expr_socp +
        part.quad_rt_expr_socp

    @objective(model, Min, part.obj_static)

    return nothing
end

# ============================================================
# BATCH
# ============================================================

function init_state_partition!(
    partitions::Partitions{PartitionSOCP},
    sys::DistributedSystem,
    optimizer_factory,
    rho_v2::Float64,
    rho_r::Float64,
    rho_t::Float64,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n "DistributedSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n
        init_state_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_v2,
            rho_r,
            rho_t,
            interzonal_branches,
        )
    end
    return nothing
end

function init_power_partition!(
    partitions::Partitions{PartitionSOCP},
    sys::DistributedSystem,
    optimizer_factory,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n "DistributedSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n
        init_power_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_p,
            rho_q,
            interzonal_branches,
        )
    end
    return nothing
end

function init_hybrid_partition!(
    partitions::Partitions{PartitionSOCP},
    sys::DistributedSystem,
    optimizer_factory,
    rho_v2::Float64,
    rho_p::Float64,
    rho_q::Float64,
    rho_r::Float64,
    rho_t::Float64,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n "DistributedSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_v2,
            rho_p,
            rho_q,
            rho_r,
            rho_t,
            interzonal_branches,
        )
    end
    return nothing
end

# ============================================================
# SOCP EXTRAS — soft anchoring (Gurobi stabilization)
# ============================================================

function soft_anchor_branch_expr!(
    part::PartitionSOCP;
    kappa_v2::Float64 = 1e-6,
    kappa_pq::Float64 = 1e-8,
    use_bus::Bool = true,
    use_pq::Bool = true,
)
    m = part.model::JuMP.Model
    quad = zero(JuMP.QuadExpr)

    if use_bus && haskey(m, :v_magnitude_sq)
        v2 = m[:v_magnitude_sq]
        part.v2_var = v2

        for I in keys(v2)
            x = v2[I]
            x_ref = (JuMP.start_value(x) === nothing) ? 1.0 : JuMP.start_value(x)
            quad += 0.5 * kappa_v2 * (x - x_ref)^2
        end
    end

    if use_pq && !isempty(part.zone_edges)
        @inbounds for e in part.zone_edges
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue

            p, q = _pq_flow_from_exprs(m, l)
            quad += 0.5 * kappa_pq * (p^2 + q^2)
        end
    end

    part.soft_anchor_expr = quad
    return quad
end
