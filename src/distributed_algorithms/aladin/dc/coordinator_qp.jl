@inline function _add_zone_cost!(
    qp_obj::JuMP.QuadExpr,
    zone::AladinZoneData,
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

function _add_active_constraints!(
    model::JuMP.Model,
    zone,
    dy::Vector{JuMP.VariableRef},
)
    C = build_active_jacobian(zone)
    nrows = size(C, 1)
    if nrows == 0
        return nothing
    end
    rows = SparseArrays.rowvals(C)
    vals = SparseArrays.nonzeros(C)

    row_exprs = [zero(JuMP.AffExpr) for _ in 1:nrows]
    @inbounds for col in 1:size(C, 2)
        for idx in SparseArrays.nzrange(C, col)
            r = rows[idx]
            JuMP.add_to_expression!(row_exprs[r], vals[idx], dy[col])
        end
    end
    @inbounds for r in 1:nrows
        JuMP.@constraint(model, row_exprs[r] == 0)
    end
    return nothing
end

# ------------------------------------------------------------
# Helper: fix Δy_i,k = 0 for each active box bound.
# ------------------------------------------------------------
function _add_active_bound_constraints!(
    model::JuMP.Model,
    zone,
    dy::Vector{JuMP.VariableRef},
)
    @inbounds for k in zone.active_bound_cols
        JuMP.@constraint(model, dy[k] == 0)
    end
    return nothing
end

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

@inline function _split_bus_coupling_expr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneData},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SplitBusPair,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    yi = zi.y[pair.pinj_col_i]
    yj = zj.y[pair.pinj_col_j]
    JuMP.add_to_expression!(expr, yi + yj)
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_i][pair.pinj_col_i])
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_j][pair.pinj_col_j])
    return expr
end

@inline function _dummy_edge_coupling_expr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneData},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::DummyBranchPairDC,
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

@inline function _edge_coupling_expr!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneData},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SharedEdgePair,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    li = pair.edge_local_i
    lj = pair.edge_local_j

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

function solve_coordinator_qp!(ws::DCStateWorkspaceALADIN)
    model = JuMP.Model(ws.qp_optimizer)
    n_zones = length(ws.zones)

    # Δy variables per zone.
    dy_all = Vector{Vector{JuMP.VariableRef}}(undef, n_zones)
    @inbounds for i in 1:n_zones
        dy_all[i] = JuMP.@variable(model, [1:ws.zones[i].n_vars], base_name = "dy$(i)")
    end

    # Slack s per bus pair.
    n_pairs = length(ws.bus_pairs)
    s = JuMP.@variable(model, [1:n_pairs], base_name = "s")

    # Objective
    qp_obj = zero(JuMP.QuadExpr)
    @inbounds for i in 1:n_zones
        _add_zone_cost!(qp_obj, ws.zones[i], dy_all[i])
    end
    mu = ws.mu[]
    @inbounds for k in 1:n_pairs
        JuMP.add_to_expression!(qp_obj, ws.vars.theta.lambda[k], s[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s[k], s[k])
    end
    JuMP.@objective(model, Min, qp_obj)

    # Equality constraints per zone
    @inbounds for i in 1:n_zones
        _add_active_constraints!(model, ws.zones[i], dy_all[i])
        _add_active_bound_constraints!(model, ws.zones[i], dy_all[i])
    end

    # Consensus constraints per pair
    cons_refs = Vector{JuMP.ConstraintRef}(undef, n_pairs)
    @inbounds for k in 1:n_pairs
        expr = zero(JuMP.AffExpr)
        _bus_coupling_expr!(expr, ws.zones, dy_all, ws.bus_pairs[k])
        cons_refs[k] = JuMP.@constraint(model, expr == s[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-state/DC coordinator QP")

    # Extract Δy into each zone's delta buffer.
    @inbounds for i in 1:n_zones
        delta = ws.zones[i].delta
        v = dy_all[i]
        for k in 1:ws.zones[i].n_vars
            delta[k] = JuMP.value(v[k])
        end
    end

    # Store QP duals of the consensus rows back into ws.vars.theta.nu.
    nu = ws.vars.theta.nu
    @inbounds for k in 1:n_pairs
        nu[k] = -JuMP.dual(cons_refs[k])
    end
    return nothing
end

function solve_coordinator_qp!(ws::DCPowerWorkspaceALADIN)
    model = JuMP.Model(ws.qp_optimizer)
    n_zones = length(ws.zones)

    dy_all = Vector{Vector{JuMP.VariableRef}}(undef, n_zones)
    @inbounds for i in 1:n_zones
        dy_all[i] = JuMP.@variable(model, [1:ws.zones[i].n_vars], base_name = "dy$(i)")
    end

    n_pairs = length(ws.edge_pairs)
    n_dummy = length(ws.dummy_branch_pairs)
    s = JuMP.@variable(model, [1:n_pairs], base_name = "s")
    s_dummy = JuMP.@variable(model, [1:n_dummy], base_name = "sd")

    qp_obj = zero(JuMP.QuadExpr)
    @inbounds for i in 1:n_zones
        _add_zone_cost!(qp_obj, ws.zones[i], dy_all[i])
    end
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

    @inbounds for i in 1:n_zones
        _add_active_constraints!(model, ws.zones[i], dy_all[i])
        _add_active_bound_constraints!(model, ws.zones[i], dy_all[i])
    end

    cons_refs = Vector{JuMP.ConstraintRef}(undef, n_pairs)
    @inbounds for k in 1:n_pairs
        expr = zero(JuMP.AffExpr)
        _edge_coupling_expr!(expr, ws.zones, dy_all, ws.edge_pairs[k])
        cons_refs[k] = JuMP.@constraint(model, expr == s[k])
    end

    cons_dummy_refs = Vector{JuMP.ConstraintRef}(undef, n_dummy)
    @inbounds for k in 1:n_dummy
        expr = zero(JuMP.AffExpr)
        _dummy_edge_coupling_expr!(expr, ws.zones, dy_all, ws.dummy_branch_pairs[k])
        cons_dummy_refs[k] = JuMP.@constraint(model, expr == s_dummy[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-power/DC coordinator QP")

    @inbounds for i in 1:n_zones
        delta = ws.zones[i].delta
        v = dy_all[i]
        for k in 1:ws.zones[i].n_vars
            delta[k] = JuMP.value(v[k])
        end
    end

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

function solve_coordinator_qp!(ws::DCHybridWorkspaceALADIN)
    model = JuMP.Model(ws.qp_optimizer)
    n_zones = length(ws.zones)

    dy_all = Vector{Vector{JuMP.VariableRef}}(undef, n_zones)
    @inbounds for i in 1:n_zones
        dy_all[i] = JuMP.@variable(model, [1:ws.zones[i].n_vars], base_name = "dy$(i)")
    end

    n_bus = length(ws.bus_pairs)
    n_edge = length(ws.edge_pairs)
    n_split = length(ws.split_bus_pairs)
    n_dummy = length(ws.dummy_branch_pairs)
    s_bus = JuMP.@variable(model, [1:n_bus], base_name = "sb")
    s_edge = JuMP.@variable(model, [1:n_edge], base_name = "se")
    s_split = JuMP.@variable(model, [1:n_split], base_name = "sp")
    s_dummy = JuMP.@variable(model, [1:n_dummy], base_name = "sd")

    qp_obj = zero(JuMP.QuadExpr)
    @inbounds for i in 1:n_zones
        _add_zone_cost!(qp_obj, ws.zones[i], dy_all[i])
    end
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

    @inbounds for i in 1:n_zones
        _add_active_constraints!(model, ws.zones[i], dy_all[i])
        _add_active_bound_constraints!(model, ws.zones[i], dy_all[i])
    end

    cons_bus_refs = Vector{JuMP.ConstraintRef}(undef, n_bus)
    @inbounds for k in 1:n_bus
        expr = zero(JuMP.AffExpr)
        _bus_coupling_expr!(expr, ws.zones, dy_all, ws.bus_pairs[k])
        cons_bus_refs[k] = JuMP.@constraint(model, expr == s_bus[k])
    end

    cons_edge_refs = Vector{JuMP.ConstraintRef}(undef, n_edge)
    @inbounds for k in 1:n_edge
        expr = zero(JuMP.AffExpr)
        _edge_coupling_expr!(expr, ws.zones, dy_all, ws.edge_pairs[k])
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
        _dummy_edge_coupling_expr!(expr, ws.zones, dy_all, ws.dummy_branch_pairs[k])
        cons_dummy_refs[k] = JuMP.@constraint(model, expr == s_dummy[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-hybrid/DC coordinator QP")

    @inbounds for i in 1:n_zones
        delta = ws.zones[i].delta
        v = dy_all[i]
        for k in 1:ws.zones[i].n_vars
            delta[k] = JuMP.value(v[k])
        end
    end

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

# ------------------------------------------------------------
# Termination-status guard
# ------------------------------------------------------------

function _check_qp_status!(model::JuMP.Model, tag::AbstractString)
    ts = JuMP.termination_status(model)
    acceptable = ts == MOI.OPTIMAL || ts == MOI.LOCALLY_SOLVED
    if !acceptable
        raw = try
            JuMP.raw_status(model)
        catch
            "raw_status unavailable"
        end
        error("""
        $tag failed
          termination_status = $ts
          raw_status         = $raw
        """)
    end
    return nothing
end
