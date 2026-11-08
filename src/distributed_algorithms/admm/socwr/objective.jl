# ============================================================
# Quadratic terms
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function build_quad_bus_expr!(part::PartitionSOCWR, rho_v2::Float64)
    v2 = part.v2_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    quad = zero(JuMP.QuadExpr)

    @inbounds for b in part.zone_shared_bus_ids
        quad += 0.5 * rho_v2 * v2[b, 1]^2
    end

    part.quad_bus_expr = quad
    return quad
end

function build_quad_branch_rt_expr!(
    part::PartitionSOCWR,
    rho_r::Float64,
    rho_t::Float64,
    zone_edges::Vector{Int},
)
    m = part.local_model::JuMP.Model
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue
        r, t = _rt_from_exprs(m, l)
        quad += 0.5 * rho_r * r^2 + 0.5 * rho_t * t^2
    end

    part.quad_rt_expr = quad
    return quad
end

# ------------------------------------------------------------
# Power — SOCWRNL lifted flow penalty (opt-in)
# ------------------------------------------------------------

const _SOCWRNL_LIFTED_FLOW_ENV = "LIRION_SOCWRNL_LIFTED_FLOW"
const _SOCWRNL_LIFT_KEY = :socwrnl_lifted_pq

@inline function _socwrnl_lifted_flow_env_enabled()::Bool
    v = lowercase(strip(get(ENV, _SOCWRNL_LIFTED_FLOW_ENV, "")))
    return v == "1" || v == "true" || v == "yes" || v == "on"
end

# Lifted carrier store: (zone branch l, side) -> (vp, vq), kept in model.ext.
const _LiftStore = Dict{Tuple{Int, Symbol}, Tuple{JuMP.VariableRef, JuMP.VariableRef}}

# Serial (quadratic builder only); creates the store on first use.
@inline function _lift_store!(part::PartitionSOCWR)::Union{Nothing, _LiftStore}
    m = part.local_model::JuMP.Model
    s = get(m.ext, _SOCWRNL_LIFT_KEY, nothing)
    s === nothing || return s::_LiftStore
    (part.zone_model === SOCWRNL && _socwrnl_lifted_flow_env_enabled()) || return nothing
    s = _LiftStore()
    m.ext[_SOCWRNL_LIFT_KEY] = s
    return s
end

# Read-only lookup.
@inline function _lift_store(part::PartitionSOCWR)::Union{Nothing, _LiftStore}
    s = get((part.local_model::JuMP.Model).ext, _SOCWRNL_LIFT_KEY, nothing)
    return s === nothing ? nothing : s::_LiftStore
end

@inline _pq_flow_side_exprs(m::JuMP.Model, l::Int, side::Symbol) =
    side === :to ? _pq_flow_to_exprs(m, l) : _pq_flow_from_exprs(m, l)

# Penalty carriers for the QUADRATIC builder: may create the lifted pair.
function _pq_penalty_carriers!(part::PartitionSOCWR, l::Int, side::Symbol)
    m = part.local_model::JuMP.Model
    store = _lift_store!(part)
    store === nothing && return _pq_flow_side_exprs(m, l, side)
    key = (l, side)
    vv = get(store, key, nothing)
    vv === nothing || return vv
    p, q = _pq_flow_side_exprs(m, l, side)
    vp = JuMP.@variable(m, base_name = "lift_p[$(l),$(side)]")
    vq = JuMP.@variable(m, base_name = "lift_q[$(l),$(side)]")
    # Linear equalities tying the carriers to the Ohm's-law flow of the zone state.
    JuMP.@constraint(m, vp == p)
    JuMP.@constraint(m, vq == q)
    # explicit start (Ipopt starts from VariablePrimalStart)
    JuMP.set_start_value(vp, 0.0)
    JuMP.set_start_value(vq, 0.0)
    store[key] = (vp, vq)
    return (vp, vq)
end

# Penalty carriers for the LINEAR builders: read-only, fail-closed.
function _pq_penalty_carriers(part::PartitionSOCWR, l::Int, side::Symbol)
    m = part.local_model::JuMP.Model
    store = _lift_store(part)
    store === nothing && return _pq_flow_side_exprs(m, l, side)
    vv = get(store, (l, side), nothing)
    vv === nothing && error(
        "SOCWRNL lifted flow: no lifted carrier for zone $(part.id) branch l=$(l) " *
        "side=$(side) — the quadratic builder must run before the linear term.",
    )
    return vv
end

# True when the lift is active on at least one zone of `parts`.
function _socwrnl_lifted_flow_active(parts::Partitions{PartitionSOCWR})::Bool
    for part in parts.parts
        _lift_store(part) === nothing || return true
    end
    return false
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_quad_branch_pq_expr!(
    part::PartitionSOCWR,
    rho_p::Float64,
    rho_q::Float64,
    zone_edges::Vector{Int},
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    quad = zero(JuMP.QuadExpr)

    @inbounds for e in zone_edges
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue
        p, q = _pq_penalty_carriers!(part, l, :from)
        quad += 0.5 * rho_p * p^2 + 0.5 * rho_q * q^2
    end
    pos = part.position
    @inbounds for (ga, za, gb, zb) in branch_pairs
        if za == pos
            l = part.branch_idx_full_to_zone[ga]
            if l != 0
                # half-A: dummy bus is its TO end -> TO-side flow (see _pq_flow_to_exprs).
                p, q = _pq_penalty_carriers!(part, l, :to)
                quad += 0.5 * rho_p * p^2 + 0.5 * rho_q * q^2
            end
        end
        if zb == pos
            l = part.branch_idx_full_to_zone[gb]
            if l != 0
                # half-B: dummy bus is its FROM end -> from-side is correct.
                p, q = _pq_penalty_carriers!(part, l, :from)
                quad += 0.5 * rho_p * p^2 + 0.5 * rho_q * q^2
            end
        end
    end

    part.quad_pq_expr = quad
    return quad
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

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
# Linear terms
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function build_linear_bus_term(
    part::PartitionSOCWR,
    v2_consensus::Vector{Float64},
    v2_duals::AbstractVector{Float64},
    rho_v2::Float64,
)::JuMP.AffExpr
    v2 = part.v2_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)

    @inbounds for i in eachindex(part.zone_shared_bus_ids)
        b = part.zone_shared_bus_ids[i]
        j = part.shared_bus_consensus_idx[i]
        lin += (v2_duals[j] - rho_v2 * v2_consensus[j]) * v2[b, 1]
    end

    return lin
end

function build_linear_branch_rt_term(
    part::PartitionSOCWR,
    zone_edges::Vector{Int},
    r_consensus::Vector{Float64},
    t_consensus::Vector{Float64},
    duals_r::AbstractVector{Float64},
    duals_t::AbstractVector{Float64},
    interzonal::Vector{Int},
    rho_r::Float64,
    rho_t::Float64,
)::JuMP.AffExpr
    m = part.local_model::JuMP.Model
    lin = zero(JuMP.AffExpr)

    @inbounds for (k, e) in enumerate(interzonal)
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue
        r, t = _rt_from_exprs(m, l)

        r_cons = r_consensus[k]
        t_cons = t_consensus[k]
        r_dual = duals_r[k]
        t_dual = duals_t[k]

        lin += (r_dual - rho_r * r_cons) * r + (t_dual - rho_t * t_cons) * t
    end

    return lin
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_linear_branch_pq_term(
    part::PartitionSOCWR,
    p_consensus::Vector{Float64},
    q_consensus::Vector{Float64},
    duals_p::AbstractVector{Float64},
    duals_q::AbstractVector{Float64},
    interzonal::Vector{Int},
    rho_p::Float64,
    rho_q::Float64,
)::JuMP.AffExpr
    lin = zero(JuMP.AffExpr)

    @inbounds for (k, e) in enumerate(interzonal)
        l = part.branch_idx_full_to_zone[e]
        l == 0 && continue
        p, q = _pq_penalty_carriers(part, l, :from)

        p_cons = p_consensus[k]
        q_cons = q_consensus[k]
        p_dual = duals_p[k]
        q_dual = duals_q[k]

        lin += (p_dual - rho_p * p_cons) * p + (q_dual - rho_q * q_cons) * q
    end

    return lin
end

function build_linear_dummy_pq_term(
    part::PartitionSOCWR,
    p_consensus::Vector{Float64},
    q_consensus::Vector{Float64},
    duals_p::AbstractVector{Float64},
    duals_q::AbstractVector{Float64},
    branch_pairs::Vector{NTuple{4, Int}},
    rho_p::Float64,
    rho_q::Float64,
)::JuMP.AffExpr
    lin = zero(JuMP.AffExpr)
    isempty(branch_pairs) && return lin

    pos = part.position
    @inbounds for (k, (ga, za, gb, zb)) in enumerate(branch_pairs)
        p_cons = p_consensus[k]
        q_cons = q_consensus[k]
        p_dual = duals_p[k]
        q_dual = duals_q[k]
        if za == pos
            l = part.branch_idx_full_to_zone[ga]
            if l != 0
                # half-A: dummy bus is its TO end -> TO-side flow (see _pq_flow_to_exprs).
                p, q = _pq_penalty_carriers(part, l, :to)
                lin += (p_dual - rho_p * p_cons) * p +
                       (q_dual - rho_q * q_cons) * q
            end
        end
        if zb == pos
            l = part.branch_idx_full_to_zone[gb]
            if l != 0
                # half-B: dummy bus is its FROM end -> from-side is correct.
                p, q = _pq_penalty_carriers(part, l, :from)
                lin += (p_dual - rho_p * p_cons) * p +
                       (q_dual - rho_q * q_cons) * q
            end
        end
    end

    return lin
end

# ============================================================
# Update objective
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function update_objective!(part::PartitionSOCWR, ws::SOCWRStateWorkspace)
    p = part.position

    lin_bus = build_linear_bus_term(
        part,
        ws.vars.v2.consensus,
        view(ws.vars.v2.duals, :, p),
        ws.vars.v2.rho[],
    )

    lin_rt = build_linear_branch_rt_term(
        part,
        part.zone_edges,
        ws.vars.r.consensus,
        ws.vars.t.consensus,
        view(ws.vars.r.duals, :, p),
        view(ws.vars.t.duals, :, p),
        ws.interzonal,
        ws.vars.r.rho[],
        ws.vars.t.rho[],
    )

    @objective(part.local_model, Min, part.obj_static + lin_bus + lin_rt)
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function update_objective!(part::PartitionSOCWR, ws::SOCWRPowerWorkspace)
    p = part.position

    lin_pq = build_linear_branch_pq_term(
        part,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        view(ws.vars.p_flow.duals, :, p),
        view(ws.vars.q_flow.duals, :, p),
        ws.interzonal,
        ws.vars.p_flow.rho[],
        ws.vars.q_flow.rho[],
    )
    lin_dummy = build_linear_dummy_pq_term(
        part,
        ws.vars.p_dummy.consensus,
        ws.vars.q_dummy.consensus,
        view(ws.vars.p_dummy.duals, :, p),
        view(ws.vars.q_dummy.duals, :, p),
        ws.branch_pairs,
        ws.vars.p_dummy.rho[],
        ws.vars.q_dummy.rho[],
    )

    @objective(part.local_model, Min, part.obj_static + lin_pq + lin_dummy)
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function update_objective!(
    part::PartitionSOCWR,
    ws::SOCWRHybridWorkspace{<:Union{TieLine, DummyBus}},
)
    p = part.position

    lin_bus = build_linear_bus_term(
        part,
        ws.vars.v2.consensus,
        view(ws.vars.v2.duals, :, p),
        ws.vars.v2.rho[],
    )

    lin_pq = build_linear_branch_pq_term(
        part,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        view(ws.vars.p_flow.duals, :, p),
        view(ws.vars.q_flow.duals, :, p),
        ws.interzonal,
        ws.vars.p_flow.rho[],
        ws.vars.q_flow.rho[],
    )
    lin_dummy = build_linear_dummy_pq_term(
        part,
        ws.vars.p_dummy.consensus,
        ws.vars.q_dummy.consensus,
        view(ws.vars.p_dummy.duals, :, p),
        view(ws.vars.q_dummy.duals, :, p),
        ws.branch_pairs,
        ws.vars.p_dummy.rho[],
        ws.vars.q_dummy.rho[],
    )

    lin_rt = build_linear_branch_rt_term(
        part,
        part.zone_edges,
        ws.vars.r.consensus,
        ws.vars.t.consensus,
        view(ws.vars.r.duals, :, p),
        view(ws.vars.t.duals, :, p),
        ws.interzonal,
        ws.vars.r.rho[],
        ws.vars.t.rho[],
    )

    @objective(
        part.local_model,
        Min,
        part.obj_static + lin_bus + lin_pq + lin_dummy + lin_rt,
    )
    return nothing
end

function update_objective!(
    part::PartitionSOCWR,
    ws::SOCWRHybridWorkspace{SplitBus},
)
    p = part.position

    lin_bus = build_linear_bus_term(
        part,
        ws.vars.v2.consensus,
        view(ws.vars.v2.duals, :, p),
        ws.vars.v2.rho[],
    )

    lin_p_inj = build_linear_inj_term_groups(
        p,
        ws.inj_groups,
        ws.inj_copy_zone,
        ws.p_inj_vars,
        ws.vars.p_inj.consensus,
        view(ws.vars.p_inj.duals, :, p),
        ws.vars.p_inj.rho[],
    )
    lin_q_inj = build_linear_inj_term_groups(
        p,
        ws.inj_groups,
        ws.inj_copy_zone,
        ws.q_inj_vars,
        ws.vars.q_inj.consensus,
        view(ws.vars.q_inj.duals, :, p),
        ws.vars.q_inj.rho[],
    )

    @objective(
        part.local_model, Min,
        part.obj_static + lin_bus + lin_p_inj + lin_q_inj,
    )
    return nothing
end
