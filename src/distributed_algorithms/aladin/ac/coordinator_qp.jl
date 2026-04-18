# ============================================================
# ALADIN-AC — Phase 2 Coordinator QP (canonical HFD 2016)
# ------------------------------------------------------------

@inline function _add_zone_cost!(
    qp_obj::JuMP.QuadExpr,
    zone::AladinZoneDataAC,
    dy::Vector{JuMP.VariableRef},
)
    @inbounds for k in 1:zone.n_vars
        gk = zone.g[k]
        if gk != 0.0
            JuMP.add_to_expression!(qp_obj, gk, dy[k])
        end
    end

    H = zone.hess
    row_idx = SparseArrays.rowvals(H)
    nz_vals = SparseArrays.nonzeros(H)
    @inbounds for col in 1:size(H, 2)
        for k in SparseArrays.nzrange(H, col)
            r = row_idx[k]
            r > col && continue
            v = nz_vals[k]
            v == 0.0 && continue
            if r == col
                JuMP.add_to_expression!(qp_obj, 0.5 * v, dy[r], dy[col])
            else
                JuMP.add_to_expression!(qp_obj, v, dy[r], dy[col])
            end
        end
    end
    return nothing
end

struct _ACFlowLin
    pf::Float64
    qf::Float64
    dpf_dvf::Float64
    dpf_dvt::Float64
    dpf_dθf::Float64
    dpf_dθt::Float64
    dqf_dvf::Float64
    dqf_dvt::Float64
    dqf_dθf::Float64
    dqf_dθt::Float64
end

@inline function _ac_flow_linearisation(
    zone::AladinZoneDataAC,
    le::Int,
)::_ACFlowLin
    # Pull current variable values at y
    θf = zone.y[zone.edge_from_theta_col[le]]
    θt = zone.y[zone.edge_to_theta_col[le]]
    vf = zone.y[zone.edge_from_vm_col[le]]
    vt = zone.y[zone.edge_to_vm_col[le]]

    g = zone.edge_g[le]
    b = zone.edge_b[le]
    b_sh = zone.edge_b_sh[le]
    tau = zone.edge_inv_tap[le]
    tap_angle = zone.edge_tap_angle[le]

    Δθ = θf - θt - tap_angle
    c = cos(Δθ)
    s = sin(Δθ)

    # Flow values at y
    pf = g * (tau * vf)^2 - tau * vf * vt * (g * c + b * s)
    qf = -(b + b_sh / 2) * (tau * vf)^2 - tau * vf * vt * (-b * c + g * s)

    # Gradients:  d(pf)/d(vf, vt, θf, θt)
    dpf_dvf = 2.0 * g * tau^2 * vf - tau * vt * (g * c + b * s)
    dpf_dvt = -tau * vf * (g * c + b * s)
    dpf_dθf = -tau * vf * vt * (-g * s + b * c)
    dpf_dθt = -dpf_dθf   # symmetric: d/dθt = -d/dθf

    # Gradients:  d(qf)/d(vf, vt, θf, θt)
    dqf_dvf = -2.0 * (b + b_sh / 2) * tau^2 * vf - tau * vt * (-b * c + g * s)
    dqf_dvt = -tau * vf * (-b * c + g * s)
    dqf_dθf = -tau * vf * vt * (b * s + g * c)
    dqf_dθt = -dqf_dθf

    return _ACFlowLin(
        pf, qf,
        dpf_dvf, dpf_dvt, dpf_dθf, dpf_dθt,
        dqf_dvf, dqf_dvt, dqf_dθf, dqf_dθt,
    )
end

# ------------------------------------------------------------
# Coupling expression helpers
# ------------------------------------------------------------

@inline function _bus_theta_expr_ac!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneDataAC},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SharedBusPairAC,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    yi = zi.y[pair.theta_col_i]
    yj = zj.y[pair.theta_col_j]
    JuMP.add_to_expression!(expr, yi - yj)
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_i][pair.theta_col_i])
    JuMP.add_to_expression!(expr, -1.0, dy_all[pair.zone_j][pair.theta_col_j])
    return expr
end

@inline function _bus_vm_expr_ac!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneDataAC},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SharedBusPairAC,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    yi = zi.y[pair.vm_col_i]
    yj = zj.y[pair.vm_col_j]
    JuMP.add_to_expression!(expr, yi - yj)
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_i][pair.vm_col_i])
    JuMP.add_to_expression!(expr, -1.0, dy_all[pair.zone_j][pair.vm_col_j])
    return expr
end

@inline function _pq_coupling_expr_ac!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneDataAC},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SharedEdgePairAC,
    flow::Symbol,   # :p or :q
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]

    lin_i = _ac_flow_linearisation(zi, pair.edge_local_i)
    lin_j = _ac_flow_linearisation(zj, pair.edge_local_j)

    if flow === :p
        JuMP.add_to_expression!(expr, lin_i.pf - lin_j.pf)
        # Zone i side: +∇pf_i ⋅ Δy_i
        _add_flow_grad!(expr, dy_all[pair.zone_i], zi, pair.edge_local_i,
            +1.0, lin_i.dpf_dvf, lin_i.dpf_dvt, lin_i.dpf_dθf, lin_i.dpf_dθt)
        # Zone j side: -∇pf_j ⋅ Δy_j
        _add_flow_grad!(expr, dy_all[pair.zone_j], zj, pair.edge_local_j,
            -1.0, lin_j.dpf_dvf, lin_j.dpf_dvt, lin_j.dpf_dθf, lin_j.dpf_dθt)
    else
        JuMP.add_to_expression!(expr, lin_i.qf - lin_j.qf)
        _add_flow_grad!(expr, dy_all[pair.zone_i], zi, pair.edge_local_i,
            +1.0, lin_i.dqf_dvf, lin_i.dqf_dvt, lin_i.dqf_dθf, lin_i.dqf_dθt)
        _add_flow_grad!(expr, dy_all[pair.zone_j], zj, pair.edge_local_j,
            -1.0, lin_j.dqf_dvf, lin_j.dqf_dvt, lin_j.dqf_dθf, lin_j.dqf_dθt)
    end
    return expr
end

@inline function _dummy_ac_flow_linearisation(
    y::Vector{Float64},
    from_theta_col::Int,
    to_theta_col::Int,
    from_vm_col::Int,
    to_vm_col::Int,
    g::Float64,
    b::Float64,
    b_sh::Float64,
    tau::Float64,
    tap_angle::Float64,
)::_ACFlowLin
    θf = y[from_theta_col]
    θt = y[to_theta_col]
    vf = y[from_vm_col]
    vt = y[to_vm_col]

    Δθ = θf - θt - tap_angle
    c = cos(Δθ)
    s = sin(Δθ)

    pf = g * (tau * vf)^2 - tau * vf * vt * (g * c + b * s)
    qf = -(b + b_sh / 2) * (tau * vf)^2 - tau * vf * vt * (-b * c + g * s)

    dpf_dvf = 2.0 * g * tau^2 * vf - tau * vt * (g * c + b * s)
    dpf_dvt = -tau * vf * (g * c + b * s)
    dpf_dθf = -tau * vf * vt * (-g * s + b * c)
    dpf_dθt = -dpf_dθf

    dqf_dvf = -2.0 * (b + b_sh / 2) * tau^2 * vf - tau * vt * (-b * c + g * s)
    dqf_dvt = -tau * vf * (-b * c + g * s)
    dqf_dθf = -tau * vf * vt * (b * s + g * c)
    dqf_dθt = -dqf_dθf

    return _ACFlowLin(
        pf, qf,
        dpf_dvf, dpf_dvt, dpf_dθf, dpf_dθt,
        dqf_dvf, dqf_dvt, dqf_dθf, dqf_dθt,
    )
end

@inline function _dummy_pq_coupling_expr_ac!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneDataAC},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::DummyBranchPairAC,
    flow::Symbol,  # :p or :q
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]

    lin_i = _dummy_ac_flow_linearisation(
        zi.y,
        pair.from_theta_col_i, pair.to_theta_col_i,
        pair.from_vm_col_i, pair.to_vm_col_i,
        pair.g_i, pair.b_i, pair.b_sh_i,
        pair.inv_tap_i, pair.tap_angle_i,
    )
    lin_j = _dummy_ac_flow_linearisation(
        zj.y,
        pair.from_theta_col_j, pair.to_theta_col_j,
        pair.from_vm_col_j, pair.to_vm_col_j,
        pair.g_j, pair.b_j, pair.b_sh_j,
        pair.inv_tap_j, pair.tap_angle_j,
    )

    if flow === :p
        JuMP.add_to_expression!(expr, lin_i.pf - lin_j.pf)
        dy_i = dy_all[pair.zone_i]
        dy_j = dy_all[pair.zone_j]
        JuMP.add_to_expression!(expr, lin_i.dpf_dvf, dy_i[pair.from_vm_col_i])
        JuMP.add_to_expression!(expr, lin_i.dpf_dvt, dy_i[pair.to_vm_col_i])
        JuMP.add_to_expression!(expr, lin_i.dpf_dθf, dy_i[pair.from_theta_col_i])
        JuMP.add_to_expression!(expr, lin_i.dpf_dθt, dy_i[pair.to_theta_col_i])
        JuMP.add_to_expression!(expr, -lin_j.dpf_dvf, dy_j[pair.from_vm_col_j])
        JuMP.add_to_expression!(expr, -lin_j.dpf_dvt, dy_j[pair.to_vm_col_j])
        JuMP.add_to_expression!(expr, -lin_j.dpf_dθf, dy_j[pair.from_theta_col_j])
        JuMP.add_to_expression!(expr, -lin_j.dpf_dθt, dy_j[pair.to_theta_col_j])
    else
        JuMP.add_to_expression!(expr, lin_i.qf - lin_j.qf)
        dy_i = dy_all[pair.zone_i]
        dy_j = dy_all[pair.zone_j]
        JuMP.add_to_expression!(expr, lin_i.dqf_dvf, dy_i[pair.from_vm_col_i])
        JuMP.add_to_expression!(expr, lin_i.dqf_dvt, dy_i[pair.to_vm_col_i])
        JuMP.add_to_expression!(expr, lin_i.dqf_dθf, dy_i[pair.from_theta_col_i])
        JuMP.add_to_expression!(expr, lin_i.dqf_dθt, dy_i[pair.to_theta_col_i])
        JuMP.add_to_expression!(expr, -lin_j.dqf_dvf, dy_j[pair.from_vm_col_j])
        JuMP.add_to_expression!(expr, -lin_j.dqf_dvt, dy_j[pair.to_vm_col_j])
        JuMP.add_to_expression!(expr, -lin_j.dqf_dθf, dy_j[pair.from_theta_col_j])
        JuMP.add_to_expression!(expr, -lin_j.dqf_dθt, dy_j[pair.to_theta_col_j])
    end
    return expr
end

@inline function _pinj_coupling_expr_ac!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneDataAC},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SplitBusPairAC,
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

@inline function _qinj_coupling_expr_ac!(
    expr::JuMP.AffExpr,
    zones::Vector{AladinZoneDataAC},
    dy_all::Vector{Vector{JuMP.VariableRef}},
    pair::SplitBusPairAC,
)
    zi = zones[pair.zone_i]
    zj = zones[pair.zone_j]
    yi = zi.y[pair.qinj_col_i]
    yj = zj.y[pair.qinj_col_j]
    JuMP.add_to_expression!(expr, yi + yj)
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_i][pair.qinj_col_i])
    JuMP.add_to_expression!(expr, 1.0, dy_all[pair.zone_j][pair.qinj_col_j])
    return expr
end

@inline function _add_flow_grad!(
    expr::JuMP.AffExpr,
    dy::Vector{JuMP.VariableRef},
    zone::AladinZoneDataAC,
    le::Int,
    sign::Float64,
    dvf::Float64, dvt::Float64, dθf::Float64, dθt::Float64,
)
    JuMP.add_to_expression!(expr, sign * dvf, dy[zone.edge_from_vm_col[le]])
    JuMP.add_to_expression!(expr, sign * dvt, dy[zone.edge_to_vm_col[le]])
    JuMP.add_to_expression!(expr, sign * dθf, dy[zone.edge_from_theta_col[le]])
    JuMP.add_to_expression!(expr, sign * dθt, dy[zone.edge_to_theta_col[le]])
    return nothing
end

# ------------------------------------------------------------
# Shared QP scaffolding for AC
# ------------------------------------------------------------

function _build_ac_qp_skeleton(ws::AbstractALADINWorkspace)
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

function _extract_delta!(ws::AbstractALADINWorkspace, dy_all)
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
# STATE
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::ACStateWorkspaceALADIN)
    model, dy_all, qp_obj = _build_ac_qp_skeleton(ws)

    n_bus = length(ws.bus_pairs)

    s_theta = JuMP.@variable(model, [1:n_bus], base_name = "sθ")
    s_vm = JuMP.@variable(model, [1:n_bus], base_name = "svm")

    mu = ws.mu[]
    @inbounds for k in 1:n_bus
        JuMP.add_to_expression!(qp_obj, ws.vars.theta.lambda[k], s_theta[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_theta[k], s_theta[k])
        JuMP.add_to_expression!(qp_obj, ws.vars.v_mag.lambda[k], s_vm[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_vm[k], s_vm[k])
    end
    JuMP.@objective(model, Min, qp_obj)

    cons_theta = Vector{JuMP.ConstraintRef}(undef, n_bus)
    cons_vm = Vector{JuMP.ConstraintRef}(undef, n_bus)
    @inbounds for k in 1:n_bus
        et = zero(JuMP.AffExpr)
        _bus_theta_expr_ac!(et, ws.zones, dy_all, ws.bus_pairs[k])
        cons_theta[k] = JuMP.@constraint(model, et == s_theta[k])

        ev = zero(JuMP.AffExpr)
        _bus_vm_expr_ac!(ev, ws.zones, dy_all, ws.bus_pairs[k])
        cons_vm[k] = JuMP.@constraint(model, ev == s_vm[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-state/AC coordinator QP")

    _extract_delta!(ws, dy_all)

    nu_theta = ws.vars.theta.nu
    nu_vm = ws.vars.v_mag.nu
    @inbounds for k in 1:n_bus
        nu_theta[k] = -JuMP.dual(cons_theta[k])
        nu_vm[k] = -JuMP.dual(cons_vm[k])
    end
    return nothing
end

# ------------------------------------------------------------
# POWER
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::ACPowerWorkspaceALADIN)
    model, dy_all, qp_obj = _build_ac_qp_skeleton(ws)

    n_edge = length(ws.edge_pairs)
    n_dummy = length(ws.dummy_branch_pairs)
    s_p = JuMP.@variable(model, [1:n_edge], base_name = "sp")
    s_q = JuMP.@variable(model, [1:n_edge], base_name = "sq")
    s_pd = JuMP.@variable(model, [1:n_dummy], base_name = "spd")
    s_qd = JuMP.@variable(model, [1:n_dummy], base_name = "sqd")

    mu = ws.mu[]
    @inbounds for k in 1:n_edge
        JuMP.add_to_expression!(qp_obj, ws.vars.p_flow.lambda[k], s_p[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_p[k], s_p[k])
        JuMP.add_to_expression!(qp_obj, ws.vars.q_flow.lambda[k], s_q[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_q[k], s_q[k])
    end
    @inbounds for k in 1:n_dummy
        JuMP.add_to_expression!(qp_obj, ws.vars.p_dummy.lambda[k], s_pd[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_pd[k], s_pd[k])
        JuMP.add_to_expression!(qp_obj, ws.vars.q_dummy.lambda[k], s_qd[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_qd[k], s_qd[k])
    end
    JuMP.@objective(model, Min, qp_obj)

    cons_p = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_q = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_pd = Vector{JuMP.ConstraintRef}(undef, n_dummy)
    cons_qd = Vector{JuMP.ConstraintRef}(undef, n_dummy)
    @inbounds for k in 1:n_edge
        ep = zero(JuMP.AffExpr)
        _pq_coupling_expr_ac!(ep, ws.zones, dy_all, ws.edge_pairs[k], :p)
        cons_p[k] = JuMP.@constraint(model, ep == s_p[k])

        eq = zero(JuMP.AffExpr)
        _pq_coupling_expr_ac!(eq, ws.zones, dy_all, ws.edge_pairs[k], :q)
        cons_q[k] = JuMP.@constraint(model, eq == s_q[k])
    end
    @inbounds for k in 1:n_dummy
        epd = zero(JuMP.AffExpr)
        _dummy_pq_coupling_expr_ac!(epd, ws.zones, dy_all, ws.dummy_branch_pairs[k], :p)
        cons_pd[k] = JuMP.@constraint(model, epd == s_pd[k])

        eqd = zero(JuMP.AffExpr)
        _dummy_pq_coupling_expr_ac!(eqd, ws.zones, dy_all, ws.dummy_branch_pairs[k], :q)
        cons_qd[k] = JuMP.@constraint(model, eqd == s_qd[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-power/AC coordinator QP")

    _extract_delta!(ws, dy_all)

    nu_p = ws.vars.p_flow.nu
    nu_q = ws.vars.q_flow.nu
    @inbounds for k in 1:n_edge
        nu_p[k] = -JuMP.dual(cons_p[k])
        nu_q[k] = -JuMP.dual(cons_q[k])
    end
    nu_pd = ws.vars.p_dummy.nu
    nu_qd = ws.vars.q_dummy.nu
    @inbounds for k in 1:n_dummy
        nu_pd[k] = -JuMP.dual(cons_pd[k])
        nu_qd[k] = -JuMP.dual(cons_qd[k])
    end
    return nothing
end

# ------------------------------------------------------------
# HYBRID
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::ACHybridWorkspaceALADIN)
    model, dy_all, qp_obj = _build_ac_qp_skeleton(ws)

    n_bus = length(ws.bus_pairs)
    n_edge = length(ws.edge_pairs)
    n_split = length(ws.split_bus_pairs)
    n_dummy = length(ws.dummy_branch_pairs)

    s_theta = JuMP.@variable(model, [1:n_bus], base_name = "sθ")
    s_vm = JuMP.@variable(model, [1:n_bus], base_name = "svm")
    s_p = JuMP.@variable(model, [1:n_edge], base_name = "sp")
    s_q = JuMP.@variable(model, [1:n_edge], base_name = "sq")
    s_pinj = JuMP.@variable(model, [1:n_split], base_name = "spi")
    s_qinj = JuMP.@variable(model, [1:n_split], base_name = "sqi")
    s_pd = JuMP.@variable(model, [1:n_dummy], base_name = "spd")
    s_qd = JuMP.@variable(model, [1:n_dummy], base_name = "sqd")

    mu = ws.mu[]
    @inbounds for k in 1:n_bus
        JuMP.add_to_expression!(qp_obj, ws.vars.theta.lambda[k], s_theta[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_theta[k], s_theta[k])
        JuMP.add_to_expression!(qp_obj, ws.vars.v_mag.lambda[k], s_vm[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_vm[k], s_vm[k])
    end
    @inbounds for k in 1:n_edge
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
    @inbounds for k in 1:n_dummy
        JuMP.add_to_expression!(qp_obj, ws.vars.p_dummy.lambda[k], s_pd[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_pd[k], s_pd[k])
        JuMP.add_to_expression!(qp_obj, ws.vars.q_dummy.lambda[k], s_qd[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_qd[k], s_qd[k])
    end
    JuMP.@objective(model, Min, qp_obj)

    cons_theta = Vector{JuMP.ConstraintRef}(undef, n_bus)
    cons_vm = Vector{JuMP.ConstraintRef}(undef, n_bus)
    @inbounds for k in 1:n_bus
        et = zero(JuMP.AffExpr)
        _bus_theta_expr_ac!(et, ws.zones, dy_all, ws.bus_pairs[k])
        cons_theta[k] = JuMP.@constraint(model, et == s_theta[k])

        ev = zero(JuMP.AffExpr)
        _bus_vm_expr_ac!(ev, ws.zones, dy_all, ws.bus_pairs[k])
        cons_vm[k] = JuMP.@constraint(model, ev == s_vm[k])
    end

    cons_p = Vector{JuMP.ConstraintRef}(undef, n_edge)
    cons_q = Vector{JuMP.ConstraintRef}(undef, n_edge)
    @inbounds for k in 1:n_edge
        ep = zero(JuMP.AffExpr)
        _pq_coupling_expr_ac!(ep, ws.zones, dy_all, ws.edge_pairs[k], :p)
        cons_p[k] = JuMP.@constraint(model, ep == s_p[k])

        eq = zero(JuMP.AffExpr)
        _pq_coupling_expr_ac!(eq, ws.zones, dy_all, ws.edge_pairs[k], :q)
        cons_q[k] = JuMP.@constraint(model, eq == s_q[k])
    end

    cons_pinj = Vector{JuMP.ConstraintRef}(undef, n_split)
    cons_qinj = Vector{JuMP.ConstraintRef}(undef, n_split)
    @inbounds for k in 1:n_split
        epi = zero(JuMP.AffExpr)
        _pinj_coupling_expr_ac!(epi, ws.zones, dy_all, ws.split_bus_pairs[k])
        cons_pinj[k] = JuMP.@constraint(model, epi == s_pinj[k])

        eqi = zero(JuMP.AffExpr)
        _qinj_coupling_expr_ac!(eqi, ws.zones, dy_all, ws.split_bus_pairs[k])
        cons_qinj[k] = JuMP.@constraint(model, eqi == s_qinj[k])
    end

    cons_pd = Vector{JuMP.ConstraintRef}(undef, n_dummy)
    cons_qd = Vector{JuMP.ConstraintRef}(undef, n_dummy)
    @inbounds for k in 1:n_dummy
        epd = zero(JuMP.AffExpr)
        _dummy_pq_coupling_expr_ac!(epd, ws.zones, dy_all, ws.dummy_branch_pairs[k], :p)
        cons_pd[k] = JuMP.@constraint(model, epd == s_pd[k])

        eqd = zero(JuMP.AffExpr)
        _dummy_pq_coupling_expr_ac!(eqd, ws.zones, dy_all, ws.dummy_branch_pairs[k], :q)
        cons_qd[k] = JuMP.@constraint(model, eqd == s_qd[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-hybrid/AC coordinator QP")

    _extract_delta!(ws, dy_all)

    nu_theta = ws.vars.theta.nu
    nu_vm = ws.vars.v_mag.nu
    nu_p = ws.vars.p_flow.nu
    nu_q = ws.vars.q_flow.nu
    nu_pinj = ws.vars.p_inj.nu
    nu_qinj = ws.vars.q_inj.nu
    @inbounds for k in 1:n_bus
        nu_theta[k] = -JuMP.dual(cons_theta[k])
        nu_vm[k] = -JuMP.dual(cons_vm[k])
    end
    @inbounds for k in 1:n_edge
        nu_p[k] = -JuMP.dual(cons_p[k])
        nu_q[k] = -JuMP.dual(cons_q[k])
    end
    @inbounds for k in 1:n_split
        nu_pinj[k] = -JuMP.dual(cons_pinj[k])
        nu_qinj[k] = -JuMP.dual(cons_qinj[k])
    end
    nu_pd = ws.vars.p_dummy.nu
    nu_qd = ws.vars.q_dummy.nu
    @inbounds for k in 1:n_dummy
        nu_pd[k] = -JuMP.dual(cons_pd[k])
        nu_qd[k] = -JuMP.dual(cons_qd[k])
    end
    return nothing
end
