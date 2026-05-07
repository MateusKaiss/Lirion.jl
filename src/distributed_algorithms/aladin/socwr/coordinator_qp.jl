# ============================================================
# ALADIN-SOCWR — Phase 2 Coordinator QP
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
    pair::InterzonalLinePairSOCWR,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    ri_col = zi.edge_r_cols[pair.line_local_i]
    rj_col = zj.edge_r_cols[pair.line_local_j]
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
    pair::InterzonalLinePairSOCWR,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    ti_col = zi.edge_t_cols[pair.line_local_i]
    tj_col = zj.edge_t_cols[pair.line_local_j]
    yi = zi.y[ti_col]
    yj = zj.y[tj_col]
    JuMP.add_to_expression!(expr, yi - yj)
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_i][ti_col])
    JuMP.add_to_expression!(expr, -1.0, dy_all[pair.zone_j][tj_col])
    return expr
end

@inline function _pinj_coupling_expr_socwr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneDataSOCWR},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SplitBusPairSOCWR,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    yi = zi.y[pair.p_inj_col_i]
    yj = zj.y[pair.p_inj_col_j]
    JuMP.add_to_expression!(expr, yi + yj)
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_i][pair.p_inj_col_i])
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_j][pair.p_inj_col_j])
    return expr
end

@inline function _qinj_coupling_expr_socwr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneDataSOCWR},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SplitBusPairSOCWR,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    yi = zi.y[pair.q_inj_col_i]
    yj = zj.y[pair.q_inj_col_j]
    JuMP.add_to_expression!(expr, yi + yj)
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_i][pair.q_inj_col_i])
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_j][pair.q_inj_col_j])
    return expr
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
    pair::InterzonalLinePairSOCWR,
    which::Symbol,   # :p or :q
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]

    fi = _socwr_flow_value_at_y(zi, pair.line_local_i, which)
    fj = _socwr_flow_value_at_y(zj, pair.line_local_j, which)
    JuMP.add_to_expression!(expr, fi - fj)

    _add_socwr_flow_grad!(expr, dy_all[pair.zone_i], zi, pair.line_local_i, +1.0, which)
    _add_socwr_flow_grad!(expr, dy_all[pair.zone_j], zj, pair.line_local_j, -1.0, which)
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
    n_edge = length(ws.interzonal_line_pairs)

    s_v2 = JuMP.@variable(model, [1:n_bus], base_name = "sv2")
    s_r = JuMP.@variable(model, [1:n_edge], base_name = "sr")
    s_t = JuMP.@variable(model, [1:n_edge], base_name = "st")

    mu = ws.mu[]
    @inbounds for k in 1:n_bus
        JuMP.add_to_expression!(qp_obj, ws.vars.v2.lambda[k], s_v2[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_v2[k], s_v2[k])
    end
    @inbounds for k in 1:n_edge
        JuMP.add_to_expression!(qp_obj, ws.vars.r.lambda[k], s_r[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_r[k], s_r[k])
        JuMP.add_to_expression!(qp_obj, ws.vars.t.lambda[k], s_t[k])
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
        _edge_r_expr_socwr!(er, ws.zones, dy_all, ws.interzonal_line_pairs[k])
        cons_r[k] = JuMP.@constraint(model, er == s_r[k])

        et = zero(JuMP.AffExpr)
        _edge_t_expr_socwr!(et, ws.zones, dy_all, ws.interzonal_line_pairs[k])
        cons_t[k] = JuMP.@constraint(model, et == s_t[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-state/SOCWR coordinator QP")

    _extract_delta_socwr!(ws, dy_all)

    nu_v2 = ws.vars.v2.nu
    nu_r = ws.vars.r.nu
    nu_t = ws.vars.t.nu
    @inbounds for k in 1:n_bus
        nu_v2[k] = -JuMP.dual(cons_v2[k])
    end
    @inbounds for k in 1:n_edge
        nu_r[k] = -JuMP.dual(cons_r[k])
        nu_t[k] = -JuMP.dual(cons_t[k])
    end
    return nothing
end

# ------------------------------------------------------------
# POWER — p_flow / q_flow consensus
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::SOCWRPowerWorkspaceALADIN)
    model, dy_all, qp_obj = _build_socwr_qp_skeleton(ws)

    n_edge = length(ws.interzonal_line_pairs)

    s_p = JuMP.@variable(model, [1:n_edge], base_name = "sp")
    s_q = JuMP.@variable(model, [1:n_edge], base_name = "sq")

    mu = ws.mu[]
    @inbounds for k in 1:n_edge
        JuMP.add_to_expression!(qp_obj, ws.vars.p_flow.lambda[k], s_p[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_p[k], s_p[k])
        JuMP.add_to_expression!(qp_obj, ws.vars.q_flow.lambda[k], s_q[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_q[k], s_q[k])
    end
    JuMP.@objective(model, Min, qp_obj)

    cons_p = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_q = Vector{JuMP.ConstraintRef}(undef, n_edge)
    @inbounds for k in 1:n_edge
        ep = zero(JuMP.AffExpr)
        _pq_coupling_expr_socwr!(ep, ws.zones, dy_all, ws.interzonal_line_pairs[k], :p)
        cons_p[k] = JuMP.@constraint(model, ep == s_p[k])

        eq = zero(JuMP.AffExpr)
        _pq_coupling_expr_socwr!(eq, ws.zones, dy_all, ws.interzonal_line_pairs[k], :q)
        cons_q[k] = JuMP.@constraint(model, eq == s_q[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-power/SOCWR coordinator QP")

    _extract_delta_socwr!(ws, dy_all)

    nu_p = ws.vars.p_flow.nu
    nu_q = ws.vars.q_flow.nu
    @inbounds for k in 1:n_edge
        nu_p[k] = -JuMP.dual(cons_p[k])
        nu_q[k] = -JuMP.dual(cons_q[k])
    end
    return nothing
end

# ------------------------------------------------------------
# HYBRID — v²/R/T + p_flow/q_flow
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::SOCWRHybridWorkspaceALADIN)
    model, dy_all, qp_obj = _build_socwr_qp_skeleton(ws)

    n_bus = length(ws.bus_pairs)
    n_edge = length(ws.interzonal_line_pairs)
    n_split = length(ws.split_bus_pairs)

    s_v2 = JuMP.@variable(model, [1:n_bus], base_name = "sv2")
    s_r = JuMP.@variable(model, [1:n_edge], base_name = "sr")
    s_t = JuMP.@variable(model, [1:n_edge], base_name = "st")
    s_p = JuMP.@variable(model, [1:n_edge], base_name = "sp")
    s_q = JuMP.@variable(model, [1:n_edge], base_name = "sq")
    s_pinj = JuMP.@variable(model, [1:n_split], base_name = "spi")
    s_qinj = JuMP.@variable(model, [1:n_split], base_name = "sqi")

    mu = ws.mu[]
    @inbounds for k in 1:n_bus
        JuMP.add_to_expression!(qp_obj, ws.vars.v2.lambda[k], s_v2[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_v2[k], s_v2[k])
    end
    @inbounds for k in 1:n_edge
        JuMP.add_to_expression!(qp_obj, ws.vars.r.lambda[k], s_r[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_r[k], s_r[k])
        JuMP.add_to_expression!(qp_obj, ws.vars.t.lambda[k], s_t[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_t[k], s_t[k])
        JuMP.add_to_expression!(qp_obj, ws.vars.p_flow.lambda[k], s_p[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_p[k], s_p[k])
        JuMP.add_to_expression!(qp_obj, ws.vars.q_flow.lambda[k], s_q[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_q[k], s_q[k])
    end
    @inbounds for k in 1:n_split
        JuMP.add_to_expression!(qp_obj, ws.vars.p_inj.lambda[k], s_pinj[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_pinj[k], s_pinj[k])
        JuMP.add_to_expression!(qp_obj, ws.vars.q_inj.lambda[k], s_qinj[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_qinj[k], s_qinj[k])
    end
    JuMP.@objective(model, Min, qp_obj)

    cons_v2 = Vector{JuMP.ConstraintRef}(undef, n_bus)
    cons_r = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_t = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_p = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_q = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_pinj = Vector{JuMP.ConstraintRef}(undef, n_split)
    cons_qinj = Vector{JuMP.ConstraintRef}(undef, n_split)

    @inbounds for k in 1:n_bus
        ev = zero(JuMP.AffExpr)
        _bus_v2_expr_socwr!(ev, ws.zones, dy_all, ws.bus_pairs[k])
        cons_v2[k] = JuMP.@constraint(model, ev == s_v2[k])
    end
    @inbounds for k in 1:n_edge
        er = zero(JuMP.AffExpr)
        _edge_r_expr_socwr!(er, ws.zones, dy_all, ws.interzonal_line_pairs[k])
        cons_r[k] = JuMP.@constraint(model, er == s_r[k])

        et = zero(JuMP.AffExpr)
        _edge_t_expr_socwr!(et, ws.zones, dy_all, ws.interzonal_line_pairs[k])
        cons_t[k] = JuMP.@constraint(model, et == s_t[k])

        ep = zero(JuMP.AffExpr)
        _pq_coupling_expr_socwr!(ep, ws.zones, dy_all, ws.interzonal_line_pairs[k], :p)
        cons_p[k] = JuMP.@constraint(model, ep == s_p[k])

        eq = zero(JuMP.AffExpr)
        _pq_coupling_expr_socwr!(eq, ws.zones, dy_all, ws.interzonal_line_pairs[k], :q)
        cons_q[k] = JuMP.@constraint(model, eq == s_q[k])
    end
    @inbounds for k in 1:n_split
        epi = zero(JuMP.AffExpr)
        _pinj_coupling_expr_socwr!(epi, ws.zones, dy_all, ws.split_bus_pairs[k])
        cons_pinj[k] = JuMP.@constraint(model, epi == s_pinj[k])

        eqi = zero(JuMP.AffExpr)
        _qinj_coupling_expr_socwr!(eqi, ws.zones, dy_all, ws.split_bus_pairs[k])
        cons_qinj[k] = JuMP.@constraint(model, eqi == s_qinj[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-hybrid/SOCWR coordinator QP")

    _extract_delta_socwr!(ws, dy_all)

    nu_v2 = ws.vars.v2.nu
    nu_r = ws.vars.r.nu
    nu_t = ws.vars.t.nu
    nu_p = ws.vars.p_flow.nu
    nu_q = ws.vars.q_flow.nu
    nu_pinj = ws.vars.p_inj.nu
    nu_qinj = ws.vars.q_inj.nu
    @inbounds for k in 1:n_bus
        nu_v2[k] = -JuMP.dual(cons_v2[k])
    end
    @inbounds for k in 1:n_edge
        nu_r[k] = -JuMP.dual(cons_r[k])
        nu_t[k] = -JuMP.dual(cons_t[k])
        nu_p[k] = -JuMP.dual(cons_p[k])
        nu_q[k] = -JuMP.dual(cons_q[k])
    end
    @inbounds for k in 1:n_split
        nu_pinj[k] = -JuMP.dual(cons_pinj[k])
        nu_qinj[k] = -JuMP.dual(cons_qinj[k])
    end
    return nothing
end
