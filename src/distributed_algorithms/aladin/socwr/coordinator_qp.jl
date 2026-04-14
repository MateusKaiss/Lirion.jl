# ============================================================
# ALADIN-SOCWR — Phase 2 Coordinator QP (canonical HFD 2016)
# ------------------------------------------------------------
# ============================================================

@inline function _add_zone_cost!(
    qp_obj::JuMP.QuadExpr,
    zone::AladinZoneDataSOCWR,
    dy::Vector{JuMP.VariableRef},
)
    @inbounds for k in 1:zone.n_vars
        gk = zone.g[k]
        if gk != 0.0
            JuMP.add_to_expression!(qp_obj, gk, dy[k])
        end
        hk = zone.hess_diag[k]
        if hk != 0.0
            JuMP.add_to_expression!(qp_obj, 0.5 * hk, dy[k], dy[k])
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Consensus row helpers
#   (y_i + Δy_i) - (y_j + Δy_j) = s
# ------------------------------------------------------------

@inline function _bus_v2_expr_socwr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneDataSOCWR},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SharedBusPairSOCWR,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    yi = zi.y[pair.v2_col_i]
    yj = zj.y[pair.v2_col_j]
    JuMP.add_to_expression!(expr, yi - yj)
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_i][pair.v2_col_i])
    JuMP.add_to_expression!(expr, -1.0, dy_all[pair.zone_j][pair.v2_col_j])
    return expr
end

@inline function _edge_r_expr_socwr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneDataSOCWR},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SharedEdgePairSOCWR,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    ri_col = zi.edge_r_cols[pair.edge_local_i]
    rj_col = zj.edge_r_cols[pair.edge_local_j]
    yi = zi.y[ri_col]
    yj = zj.y[rj_col]
    JuMP.add_to_expression!(expr, yi - yj)
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_i][ri_col])
    JuMP.add_to_expression!(expr, -1.0, dy_all[pair.zone_j][rj_col])
    return expr
end

@inline function _edge_t_expr_socwr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneDataSOCWR},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SharedEdgePairSOCWR,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    ti_col = zi.edge_t_cols[pair.edge_local_i]
    tj_col = zj.edge_t_cols[pair.edge_local_j]
    yi = zi.y[ti_col]
    yj = zj.y[tj_col]
    JuMP.add_to_expression!(expr, yi - yj)
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_i][ti_col])
    JuMP.add_to_expression!(expr, -1.0, dy_all[pair.zone_j][tj_col])
    return expr
end

# ------------------------------------------------------------
# Flow value / affine linearisation helpers (p and q)
# ------------------------------------------------------------
#
# SOCWR flow = Σ coef_k · x_k + const. At the current y this is a
# scalar; adding `coef_k · Δx_k` gives the flow at y + Δy, which the
# coordinator QP uses as the linearised consensus row. Because the
# flow is already affine, the linearisation is exact.

@inline function _socwr_flow_value_at_y(
    zone::AladinZoneDataSOCWR,
    le::Int,
    which::Symbol,
)::Float64
    cols = zone.edge_flow_cols[le]
    coefs = which === :p ? zone.edge_pflow_coefs[le] : zone.edge_qflow_coefs[le]
    const_val = which === :p ? zone.edge_pflow_const[le] : zone.edge_qflow_const[le]
    acc = const_val
    y = zone.y
    @inbounds for i in eachindex(cols)
        acc += coefs[i] * y[cols[i]]
    end
    return acc
end

@inline function _add_socwr_flow_grad!(
    expr::JuMP.AffExpr,
    dy::Vector{JuMP.VariableRef},
    zone::AladinZoneDataSOCWR,
    le::Int,
    sign::Float64,
    which::Symbol,
)
    cols = zone.edge_flow_cols[le]
    coefs = which === :p ? zone.edge_pflow_coefs[le] : zone.edge_qflow_coefs[le]
    @inbounds for i in eachindex(cols)
        JuMP.add_to_expression!(expr, sign * coefs[i], dy[cols[i]])
    end
    return nothing
end

@inline function _pq_coupling_expr_socwr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneDataSOCWR},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SharedEdgePairSOCWR,
    which::Symbol,   # :p or :q
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]

    fi = _socwr_flow_value_at_y(zi, pair.edge_local_i, which)
    fj = _socwr_flow_value_at_y(zj, pair.edge_local_j, which)
    JuMP.add_to_expression!(expr, fi - fj)

    _add_socwr_flow_grad!(expr, dy_all[pair.zone_i], zi, pair.edge_local_i, +1.0, which)
    _add_socwr_flow_grad!(expr, dy_all[pair.zone_j], zj, pair.edge_local_j, -1.0, which)
    return expr
end

# ------------------------------------------------------------
# QP skeleton: Δy variables, zone costs, active-set preservation
# ------------------------------------------------------------

const _SOCWRCoordWorkspace = Union{
    SOCWRStateWorkspaceALADIN,
    SOCWRPowerWorkspaceALADIN,
    SOCWRHybridWorkspaceALADIN,
}

function _build_socwr_qp_skeleton(ws::_SOCWRCoordWorkspace)
    model = JuMP.Model(ws.qp_optimizer)
    n_zones = length(ws.zones)

    dy_all = Vector{Vector{JuMP.VariableRef}}(undef, n_zones)
    @inbounds for i in 1:n_zones
        dy_all[i] =
            JuMP.@variable(model, [1:ws.zones[i].n_vars], base_name = "dy$(i)")
    end

    qp_obj = zero(JuMP.QuadExpr)
    @inbounds for i in 1:n_zones
        _add_zone_cost!(qp_obj, ws.zones[i], dy_all[i])
    end

    @inbounds for i in 1:n_zones
        _add_active_constraints!(model, ws.zones[i], dy_all[i])
        _add_active_bound_constraints!(model, ws.zones[i], dy_all[i])
    end

    return model, dy_all, qp_obj
end

function _extract_delta_socwr!(ws::_SOCWRCoordWorkspace, dy_all)
    @inbounds for i in eachindex(ws.zones)
        delta = ws.zones[i].delta
        v = dy_all[i]
        for k in 1:ws.zones[i].n_vars
            delta[k] = JuMP.value(v[k])
        end
    end
    return nothing
end

# ------------------------------------------------------------
# STATE — v² consensus + R/T consensus
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::SOCWRStateWorkspaceALADIN)
    model, dy_all, qp_obj = _build_socwr_qp_skeleton(ws)

    n_bus = length(ws.bus_pairs)
    n_edge = length(ws.edge_pairs)

    s_v2 = JuMP.@variable(model, [1:n_bus], base_name = "sv2")
    s_r = JuMP.@variable(model, [1:n_edge], base_name = "sr")
    s_t = JuMP.@variable(model, [1:n_edge], base_name = "st")

    mu = ws.mu[]
    @inbounds for k in 1:n_bus
        JuMP.add_to_expression!(qp_obj, ws.lambda_v2[k], s_v2[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_v2[k], s_v2[k])
    end
    @inbounds for k in 1:n_edge
        JuMP.add_to_expression!(qp_obj, ws.lambda_r[k], s_r[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_r[k], s_r[k])
        JuMP.add_to_expression!(qp_obj, ws.lambda_t[k], s_t[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_t[k], s_t[k])
    end
    JuMP.@objective(model, Min, qp_obj)

    cons_v2 = Vector{JuMP.ConstraintRef}(undef, n_bus)
    cons_r = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_t = Vector{JuMP.ConstraintRef}(undef, n_edge)

    @inbounds for k in 1:n_bus
        ev = zero(JuMP.AffExpr)
        _bus_v2_expr_socwr!(ev, ws.zones, dy_all, ws.bus_pairs[k])
        cons_v2[k] = JuMP.@constraint(model, ev == s_v2[k])
    end
    @inbounds for k in 1:n_edge
        er = zero(JuMP.AffExpr)
        _edge_r_expr_socwr!(er, ws.zones, dy_all, ws.edge_pairs[k])
        cons_r[k] = JuMP.@constraint(model, er == s_r[k])

        et = zero(JuMP.AffExpr)
        _edge_t_expr_socwr!(et, ws.zones, dy_all, ws.edge_pairs[k])
        cons_t[k] = JuMP.@constraint(model, et == s_t[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-state/SOCWR coordinator QP")

    _extract_delta_socwr!(ws, dy_all)

    nu_v2 = Vector{Float64}(undef, n_bus)
    nu_r = Vector{Float64}(undef, n_edge)
    nu_t = Vector{Float64}(undef, n_edge)
    @inbounds for k in 1:n_bus
        nu_v2[k] = -JuMP.dual(cons_v2[k])
    end
    @inbounds for k in 1:n_edge
        nu_r[k] = -JuMP.dual(cons_r[k])
        nu_t[k] = -JuMP.dual(cons_t[k])
    end
    return (nu_v2 = nu_v2, nu_r = nu_r, nu_t = nu_t)
end

# ------------------------------------------------------------
# POWER — p_flow / q_flow consensus
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::SOCWRPowerWorkspaceALADIN)
    model, dy_all, qp_obj = _build_socwr_qp_skeleton(ws)

    n_edge = length(ws.edge_pairs)

    s_p = JuMP.@variable(model, [1:n_edge], base_name = "sp")
    s_q = JuMP.@variable(model, [1:n_edge], base_name = "sq")

    mu = ws.mu[]
    @inbounds for k in 1:n_edge
        JuMP.add_to_expression!(qp_obj, ws.lambda_p[k], s_p[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_p[k], s_p[k])
        JuMP.add_to_expression!(qp_obj, ws.lambda_q[k], s_q[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_q[k], s_q[k])
    end
    JuMP.@objective(model, Min, qp_obj)

    cons_p = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_q = Vector{JuMP.ConstraintRef}(undef, n_edge)
    @inbounds for k in 1:n_edge
        ep = zero(JuMP.AffExpr)
        _pq_coupling_expr_socwr!(ep, ws.zones, dy_all, ws.edge_pairs[k], :p)
        cons_p[k] = JuMP.@constraint(model, ep == s_p[k])

        eq = zero(JuMP.AffExpr)
        _pq_coupling_expr_socwr!(eq, ws.zones, dy_all, ws.edge_pairs[k], :q)
        cons_q[k] = JuMP.@constraint(model, eq == s_q[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-power/SOCWR coordinator QP")

    _extract_delta_socwr!(ws, dy_all)

    nu_p = Vector{Float64}(undef, n_edge)
    nu_q = Vector{Float64}(undef, n_edge)
    @inbounds for k in 1:n_edge
        nu_p[k] = -JuMP.dual(cons_p[k])
        nu_q[k] = -JuMP.dual(cons_q[k])
    end
    return (nu_p = nu_p, nu_q = nu_q)
end

# ------------------------------------------------------------
# HYBRID — v²/R/T + p_flow/q_flow
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::SOCWRHybridWorkspaceALADIN)
    model, dy_all, qp_obj = _build_socwr_qp_skeleton(ws)

    n_bus = length(ws.bus_pairs)
    n_edge = length(ws.edge_pairs)

    s_v2 = JuMP.@variable(model, [1:n_bus], base_name = "sv2")
    s_r = JuMP.@variable(model, [1:n_edge], base_name = "sr")
    s_t = JuMP.@variable(model, [1:n_edge], base_name = "st")
    s_p = JuMP.@variable(model, [1:n_edge], base_name = "sp")
    s_q = JuMP.@variable(model, [1:n_edge], base_name = "sq")

    mu = ws.mu[]
    @inbounds for k in 1:n_bus
        JuMP.add_to_expression!(qp_obj, ws.lambda_v2[k], s_v2[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_v2[k], s_v2[k])
    end
    @inbounds for k in 1:n_edge
        JuMP.add_to_expression!(qp_obj, ws.lambda_r[k], s_r[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_r[k], s_r[k])
        JuMP.add_to_expression!(qp_obj, ws.lambda_t[k], s_t[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_t[k], s_t[k])
        JuMP.add_to_expression!(qp_obj, ws.lambda_p[k], s_p[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_p[k], s_p[k])
        JuMP.add_to_expression!(qp_obj, ws.lambda_q[k], s_q[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_q[k], s_q[k])
    end
    JuMP.@objective(model, Min, qp_obj)

    cons_v2 = Vector{JuMP.ConstraintRef}(undef, n_bus)
    cons_r = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_t = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_p = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_q = Vector{JuMP.ConstraintRef}(undef, n_edge)

    @inbounds for k in 1:n_bus
        ev = zero(JuMP.AffExpr)
        _bus_v2_expr_socwr!(ev, ws.zones, dy_all, ws.bus_pairs[k])
        cons_v2[k] = JuMP.@constraint(model, ev == s_v2[k])
    end
    @inbounds for k in 1:n_edge
        er = zero(JuMP.AffExpr)
        _edge_r_expr_socwr!(er, ws.zones, dy_all, ws.edge_pairs[k])
        cons_r[k] = JuMP.@constraint(model, er == s_r[k])

        et = zero(JuMP.AffExpr)
        _edge_t_expr_socwr!(et, ws.zones, dy_all, ws.edge_pairs[k])
        cons_t[k] = JuMP.@constraint(model, et == s_t[k])

        ep = zero(JuMP.AffExpr)
        _pq_coupling_expr_socwr!(ep, ws.zones, dy_all, ws.edge_pairs[k], :p)
        cons_p[k] = JuMP.@constraint(model, ep == s_p[k])

        eq = zero(JuMP.AffExpr)
        _pq_coupling_expr_socwr!(eq, ws.zones, dy_all, ws.edge_pairs[k], :q)
        cons_q[k] = JuMP.@constraint(model, eq == s_q[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-hybrid/SOCWR coordinator QP")

    _extract_delta_socwr!(ws, dy_all)

    nu_v2 = Vector{Float64}(undef, n_bus)
    nu_r = Vector{Float64}(undef, n_edge)
    nu_t = Vector{Float64}(undef, n_edge)
    nu_p = Vector{Float64}(undef, n_edge)
    nu_q = Vector{Float64}(undef, n_edge)
    @inbounds for k in 1:n_bus
        nu_v2[k] = -JuMP.dual(cons_v2[k])
    end
    @inbounds for k in 1:n_edge
        nu_r[k] = -JuMP.dual(cons_r[k])
        nu_t[k] = -JuMP.dual(cons_t[k])
        nu_p[k] = -JuMP.dual(cons_p[k])
        nu_q[k] = -JuMP.dual(cons_q[k])
    end
    return (
        nu_v2 = nu_v2,
        nu_r = nu_r,
        nu_t = nu_t,
        nu_p = nu_p,
        nu_q = nu_q,
    )
end
