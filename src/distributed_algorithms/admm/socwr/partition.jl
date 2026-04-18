# ============================================================
# SOCWR PARTITION INITIALIZATION — ADMM

function init_state_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    rho_v2::Float64,
    rho_r::Float64,
    rho_t::Float64,
    interzonal_branches,
)
    build_local_model!(part, sys_zone, State, optimizer_factory, interzonal_branches)
    build_quad_bus_expr!(part, rho_v2)
    build_quad_branch_rt_expr!(part, rho_r, rho_t, part.zone_edges)
    part.obj_static = part.obj_base + part.quad_bus_expr + part.quad_rt_expr
    @objective(part.model, Min, part.obj_static)
    return nothing
end

function init_power_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    build_local_model!(part, sys_zone, Power, optimizer_factory, interzonal_branches)
    if is_gurobi(part.model)
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
    build_quad_branch_pq_expr!(part, rho_p, rho_q, part.zone_edges, branch_pairs)
    part.obj_static = part.obj_base + part.quad_pq_expr + part.soft_anchor_expr
    @objective(part.model, Min, part.obj_static)
    return nothing
end

function init_hybrid_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    rho_v2::Float64,
    rho_p::Float64,
    rho_q::Float64,
    rho_r::Float64,
    rho_t::Float64,
    interzonal_branches,
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
    bus_pairs::Vector{NTuple{3, Int}} = NTuple{3, Int}[],
    p_inj_store::Union{Nothing, Dict{Int, JuMP.VariableRef}} = nothing,
    q_inj_store::Union{Nothing, Dict{Int, JuMP.VariableRef}} = nothing,
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)

    # SplitBus (Method 3): introduce free p_inj[m] and q_inj[m] at every split
    # bus in this zone, folded into the active / reactive balance. SOCWR uses
    # `v_magnitude_sq` (a single variable, not v²), so both balances are AFFINE
    # in the new injection — `set_normalized_coefficient` is sufficient.
    if !isempty(bus_pairs) && p_inj_store !== nothing && q_inj_store !== nothing
        _attach_pq_inj_variables_socwr!(
            part, bus_pairs, p_inj_store, q_inj_store,
        )
    end

    build_quad_bus_expr!(part, rho_v2)
    build_quad_branch_pq_expr!(part, rho_p, rho_q, part.zone_edges, branch_pairs)
    build_quad_branch_rt_expr!(part, rho_r, rho_t, part.zone_edges)
    build_quad_p_inj_expr_socwr!(part, rho_p, p_inj_store)
    build_quad_q_inj_expr_socwr!(part, rho_q, q_inj_store)

    part.obj_static =
        part.obj_base + part.quad_bus_expr + part.quad_pq_expr +
        part.quad_rt_expr + part.quad_p_inj_expr + part.quad_q_inj_expr
    @objective(part.model, Min, part.obj_static)
    return nothing
end

# ------------------------------------------------------------
# Create p_inj / q_inj variables and fold them into the AFFINE
# SOCWR active/reactive balance constraints at split buses.
# ------------------------------------------------------------
function _attach_pq_inj_variables_socwr!(
    part::PartitionSOCWR,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_store::Dict{Int, JuMP.VariableRef},
    q_inj_store::Dict{Int, JuMP.VariableRef},
)
    model = part.model::JuMP.Model
    pos = part.position

    zone_split_bus_ids = Int[]
    for (fbus, za, zb) in bus_pairs
        if za == pos || zb == pos
            push!(zone_split_bus_ids, fbus)
        end
    end
    isempty(zone_split_bus_ids) && return nothing

    bal_p = model[:active_power_balance]
    bal_q = model[:reactive_power_balance]

    for b in zone_split_bus_ids
        vp = @variable(model, base_name = "p_inj_$(b)")
        JuMP.set_start_value(vp, 0.0)
        p_inj_store[b] = vp
        JuMP.set_normalized_coefficient(bal_p[b, 1], vp, 1.0)

        vq = @variable(model, base_name = "q_inj_$(b)")
        JuMP.set_start_value(vq, 0.0)
        q_inj_store[b] = vq
        JuMP.set_normalized_coefficient(bal_q[b, 1], vq, 1.0)
    end
    return nothing
end

function build_quad_p_inj_expr_socwr!(
    part::PartitionSOCWR,
    rho::Float64,
    p_inj_store::Union{Nothing, Dict{Int, JuMP.VariableRef}},
)
    quad = zero(JuMP.QuadExpr)
    if p_inj_store !== nothing
        for (_, v) in p_inj_store
            quad += 0.5 * rho * v^2
        end
    end
    part.quad_p_inj_expr = quad
    return quad
end

function build_quad_q_inj_expr_socwr!(
    part::PartitionSOCWR,
    rho::Float64,
    q_inj_store::Union{Nothing, Dict{Int, JuMP.VariableRef}},
)
    quad = zero(JuMP.QuadExpr)
    if q_inj_store !== nothing
        for (_, v) in q_inj_store
            quad += 0.5 * rho * v^2
        end
    end
    part.quad_q_inj_expr = quad
    return quad
end

# ============================================================
# BATCH
# ============================================================

function init_state_partition!(
    partitions::Partitions{PartitionSOCWR},
    sys::DistributedSystem,
    optimizer_factory,
    rho_v2::Float64,
    rho_r::Float64,
    rho_t::Float64,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_state_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_v2, rho_r, rho_t, interzonal_branches,
        )
    end
    return nothing
end

function init_power_partition!(
    partitions::Partitions{PartitionSOCWR},
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

function init_hybrid_partition!(
    partitions::Partitions{PartitionSOCWR},
    sys::DistributedSystem,
    optimizer_factory,
    rho_v2::Float64,
    rho_p::Float64,
    rho_q::Float64,
    rho_r::Float64,
    rho_t::Float64,
    interzonal_branches,
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
    bus_pairs::Vector{NTuple{3, Int}} = NTuple{3, Int}[],
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Dict{Int, JuMP.VariableRef}[],
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Dict{Int, JuMP.VariableRef}[],
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    have_p_inj = length(p_inj_vars) == n
    have_q_inj = length(q_inj_vars) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_v2, rho_p, rho_q, rho_r, rho_t, interzonal_branches,
            branch_pairs,
            bus_pairs,
            have_p_inj ? p_inj_vars[i] : nothing,
            have_q_inj ? q_inj_vars[i] : nothing,
        )
    end
    return nothing
end

# ============================================================
# SOCWR EXTRAS — soft anchoring (Gurobi stabilization)
# ============================================================

function soft_anchor_branch_expr!(
    part::PartitionSOCWR;
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
