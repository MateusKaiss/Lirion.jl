# ============================================================
# ALADIN-DC — Phase 2 coordinator QP
# ============================================================

# ============================================================
# Pair coupling expressions
# ============================================================

# ------------------------------------------------------------
# Bus consensus
# ------------------------------------------------------------

@inline function _bus_coupling_expr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneData},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SharedBusPair,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    yi = zi.y[pair.col_i]
    yj = zj.y[pair.col_j]
    JuMP.add_to_expression!(expr, yi - yj)
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_i][pair.col_i])
    JuMP.add_to_expression!(expr, -1.0, dy_all[pair.zone_j][pair.col_j])
    return expr
end

# ------------------------------------------------------------
# Split-bus consensus
# ------------------------------------------------------------

@inline function _split_bus_coupling_expr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneData},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SplitBusPair,
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

# ------------------------------------------------------------
# Dummy-edge consensus
# ------------------------------------------------------------

@inline function _dummy_edge_coupling_expr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneData},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::DummyLinePairDC,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]

    pflow_yi = pair.b_i * (zi.y[pair.from_col_i] - zi.y[pair.to_col_i])
    pflow_yj = pair.b_j * (zj.y[pair.from_col_j] - zj.y[pair.to_col_j])

    JuMP.add_to_expression!(expr, pflow_yi - pflow_yj)

    JuMP.add_to_expression!(expr, pair.b_i, dy_all[pair.zone_i][pair.from_col_i])
    JuMP.add_to_expression!(expr, -pair.b_i, dy_all[pair.zone_i][pair.to_col_i])
    JuMP.add_to_expression!(expr, -pair.b_j, dy_all[pair.zone_j][pair.from_col_j])
    JuMP.add_to_expression!(expr, pair.b_j, dy_all[pair.zone_j][pair.to_col_j])
    return expr
end

# ------------------------------------------------------------
# Interzonal edge consensus
# ------------------------------------------------------------

@inline function _edge_coupling_expr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneData},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::InterzonalLinePair,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    li = pair.line_local_i
    lj = pair.line_local_j

    fcol_i, tcol_i, bi = zi.edge_from_cols[li], zi.edge_to_cols[li], zi.edge_b[li]
    fcol_j, tcol_j, bj = zj.edge_from_cols[lj], zj.edge_to_cols[lj], zj.edge_b[lj]

    pflow_yi = bi * (zi.y[fcol_i] - zi.y[tcol_i])
    pflow_yj = bj * (zj.y[fcol_j] - zj.y[tcol_j])

    JuMP.add_to_expression!(expr, pflow_yi - pflow_yj)

    JuMP.add_to_expression!(expr, bi, dy_all[pair.zone_i][fcol_i])
    JuMP.add_to_expression!(expr, -bi, dy_all[pair.zone_i][tcol_i])
    JuMP.add_to_expression!(expr, -bj, dy_all[pair.zone_j][fcol_j])
    JuMP.add_to_expression!(expr, bj, dy_all[pair.zone_j][tcol_j])
    return expr
end

# ============================================================
# Dispatch
# ============================================================

# ------------------------------------------------------------
# STATE
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::DCStateWorkspaceALADIN)
    model, dy_all, qp_obj = _build_qp_skeleton(ws)

    n_pairs = length(ws.bus_pairs)
    s = JuMP.@variable(model, [1:n_pairs], base_name = "s")

    mu = ws.mu[]
    @inbounds for k in 1:n_pairs
        JuMP.add_to_expression!(qp_obj, ws.vars.theta.lambda[k], s[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s[k], s[k])
    end
    JuMP.@objective(model, Min, qp_obj)

    cons_refs = Vector{JuMP.ConstraintRef}(undef, n_pairs)
    @inbounds for k in 1:n_pairs
        expr = zero(JuMP.AffExpr)
        _bus_coupling_expr!(expr, ws.zones, dy_all, ws.bus_pairs[k])
        cons_refs[k] = JuMP.@constraint(model, expr == s[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-state/DC coordinator QP")

    _extract_delta!(ws, dy_all)

    nu = ws.vars.theta.nu
    @inbounds for k in 1:n_pairs
        nu[k] = -JuMP.dual(cons_refs[k])
    end
    return nothing
end

# ------------------------------------------------------------
# POWER
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::DCPowerWorkspaceALADIN)
    model, dy_all, qp_obj = _build_qp_skeleton(ws)

    n_pairs = length(ws.interzonal_line_pairs)
    n_dummy = length(ws.dummy_line_pairs)
    s = JuMP.@variable(model, [1:n_pairs], base_name = "s")
    s_dummy = JuMP.@variable(model, [1:n_dummy], base_name = "sd")

    mu = ws.mu[]
    @inbounds for k in 1:n_pairs
        JuMP.add_to_expression!(qp_obj, ws.vars.p_flow.lambda[k], s[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s[k], s[k])
    end
    @inbounds for k in 1:n_dummy
        JuMP.add_to_expression!(qp_obj, ws.vars.p_dummy.lambda[k], s_dummy[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_dummy[k], s_dummy[k])
    end
    JuMP.@objective(model, Min, qp_obj)

    cons_refs = Vector{JuMP.ConstraintRef}(undef, n_pairs)
    @inbounds for k in 1:n_pairs
        expr = zero(JuMP.AffExpr)
        _edge_coupling_expr!(expr, ws.zones, dy_all, ws.interzonal_line_pairs[k])
        cons_refs[k] = JuMP.@constraint(model, expr == s[k])
    end

    cons_dummy_refs = Vector{JuMP.ConstraintRef}(undef, n_dummy)
    @inbounds for k in 1:n_dummy
        expr = zero(JuMP.AffExpr)
        _dummy_edge_coupling_expr!(expr, ws.zones, dy_all, ws.dummy_line_pairs[k])
        cons_dummy_refs[k] = JuMP.@constraint(model, expr == s_dummy[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-power/DC coordinator QP")

    _extract_delta!(ws, dy_all)

    nu = ws.vars.p_flow.nu
    @inbounds for k in 1:n_pairs
        nu[k] = -JuMP.dual(cons_refs[k])
    end
    nu_dummy = ws.vars.p_dummy.nu
    @inbounds for k in 1:n_dummy
        nu_dummy[k] = -JuMP.dual(cons_dummy_refs[k])
    end
    return nothing
end

# ------------------------------------------------------------
# HYBRID
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::DCHybridWorkspaceALADIN)
    model, dy_all, qp_obj = _build_qp_skeleton(ws)

    n_bus = length(ws.bus_pairs)
    n_edge = length(ws.interzonal_line_pairs)
    n_split = length(ws.split_bus_pairs)
    n_dummy = length(ws.dummy_line_pairs)
    s_bus = JuMP.@variable(model, [1:n_bus], base_name = "sb")
    s_edge = JuMP.@variable(model, [1:n_edge], base_name = "se")
    s_split = JuMP.@variable(model, [1:n_split], base_name = "sp")
    s_dummy = JuMP.@variable(model, [1:n_dummy], base_name = "sd")

    mu = ws.mu[]
    @inbounds for k in 1:n_bus
        JuMP.add_to_expression!(qp_obj, ws.vars.theta.lambda[k], s_bus[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_bus[k], s_bus[k])
    end
    @inbounds for k in 1:n_edge
        JuMP.add_to_expression!(qp_obj, ws.vars.p_flow.lambda[k], s_edge[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_edge[k], s_edge[k])
    end
    @inbounds for k in 1:n_split
        JuMP.add_to_expression!(qp_obj, ws.vars.p_inj.lambda[k], s_split[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_split[k], s_split[k])
    end
    @inbounds for k in 1:n_dummy
        JuMP.add_to_expression!(qp_obj, ws.vars.p_dummy.lambda[k], s_dummy[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_dummy[k], s_dummy[k])
    end
    JuMP.@objective(model, Min, qp_obj)

    cons_bus_refs = Vector{JuMP.ConstraintRef}(undef, n_bus)
    @inbounds for k in 1:n_bus
        expr = zero(JuMP.AffExpr)
        _bus_coupling_expr!(expr, ws.zones, dy_all, ws.bus_pairs[k])
        cons_bus_refs[k] = JuMP.@constraint(model, expr == s_bus[k])
    end

    cons_edge_refs = Vector{JuMP.ConstraintRef}(undef, n_edge)
    @inbounds for k in 1:n_edge
        expr = zero(JuMP.AffExpr)
        _edge_coupling_expr!(expr, ws.zones, dy_all, ws.interzonal_line_pairs[k])
        cons_edge_refs[k] = JuMP.@constraint(model, expr == s_edge[k])
    end

    cons_split_refs = Vector{JuMP.ConstraintRef}(undef, n_split)
    @inbounds for k in 1:n_split
        expr = zero(JuMP.AffExpr)
        _split_bus_coupling_expr!(expr, ws.zones, dy_all, ws.split_bus_pairs[k])
        cons_split_refs[k] = JuMP.@constraint(model, expr == s_split[k])
    end

    cons_dummy_refs = Vector{JuMP.ConstraintRef}(undef, n_dummy)
    @inbounds for k in 1:n_dummy
        expr = zero(JuMP.AffExpr)
        _dummy_edge_coupling_expr!(expr, ws.zones, dy_all, ws.dummy_line_pairs[k])
        cons_dummy_refs[k] = JuMP.@constraint(model, expr == s_dummy[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-hybrid/DC coordinator QP")

    _extract_delta!(ws, dy_all)

    nu_bus = ws.vars.theta.nu
    @inbounds for k in 1:n_bus
        nu_bus[k] = -JuMP.dual(cons_bus_refs[k])
    end
    nu_edge = ws.vars.p_flow.nu
    @inbounds for k in 1:n_edge
        nu_edge[k] = -JuMP.dual(cons_edge_refs[k])
    end
    nu_split = ws.vars.p_inj.nu
    @inbounds for k in 1:n_split
        nu_split[k] = -JuMP.dual(cons_split_refs[k])
    end
    nu_dummy = ws.vars.p_dummy.nu
    @inbounds for k in 1:n_dummy
        nu_dummy[k] = -JuMP.dual(cons_dummy_refs[k])
    end
    return nothing
end
