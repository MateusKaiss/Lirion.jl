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

# ------------------------------------------------------------
# Closed-form branch flow values + gradients for the coupling rows
# ------------------------------------------------------------
#
# Lirion's AC forward flows (see src/modeling/branch/ac.jl):
#
#   p_f = g · (inv_tap · v_f)²
#        − inv_tap · v_f · v_t · ( g·cos(Δθ) + b·sin(Δθ) )
#
#   q_f = −(b + b_sh/2) · (inv_tap · v_f)²
#        − inv_tap · v_f · v_t · (−b·cos(Δθ) + g·sin(Δθ) )
#
# where Δθ = θ_f − θ_t − tap_angle.
#
# We need p_f(y), q_f(y) and their gradients w.r.t. the four
# variables (v_f, v_t, θ_f, θ_t) to assemble the Phase-2 consensus
# rows. All derivatives are obtained by hand from the formulas above
# — no autodiff is needed for this linearisation step.

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
        JuMP.add_to_expression!(qp_obj, ws.lambda_theta[k], s_theta[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_theta[k], s_theta[k])
        JuMP.add_to_expression!(qp_obj, ws.lambda_vm[k], s_vm[k])
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

    nu_theta = Vector{Float64}(undef, n_bus)
    nu_vm = Vector{Float64}(undef, n_bus)
    @inbounds for k in 1:n_bus
        nu_theta[k] = -JuMP.dual(cons_theta[k])
        nu_vm[k] = -JuMP.dual(cons_vm[k])
    end
    return (nu_theta = nu_theta, nu_vm = nu_vm)
end

# ------------------------------------------------------------
# POWER
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::ACPowerWorkspaceALADIN)
    model, dy_all, qp_obj = _build_ac_qp_skeleton(ws)

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
        _pq_coupling_expr_ac!(ep, ws.zones, dy_all, ws.edge_pairs[k], :p)
        cons_p[k] = JuMP.@constraint(model, ep == s_p[k])

        eq = zero(JuMP.AffExpr)
        _pq_coupling_expr_ac!(eq, ws.zones, dy_all, ws.edge_pairs[k], :q)
        cons_q[k] = JuMP.@constraint(model, eq == s_q[k])
    end

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-power/AC coordinator QP")

    _extract_delta!(ws, dy_all)

    nu_p = Vector{Float64}(undef, n_edge)
    nu_q = Vector{Float64}(undef, n_edge)
    @inbounds for k in 1:n_edge
        nu_p[k] = -JuMP.dual(cons_p[k])
        nu_q[k] = -JuMP.dual(cons_q[k])
    end
    return (nu_p = nu_p, nu_q = nu_q)
end

# ------------------------------------------------------------
# HYBRID
# ------------------------------------------------------------

function solve_coordinator_qp!(ws::ACHybridWorkspaceALADIN)
    model, dy_all, qp_obj = _build_ac_qp_skeleton(ws)

    n_bus = length(ws.bus_pairs)
    n_edge = length(ws.edge_pairs)

    s_theta = JuMP.@variable(model, [1:n_bus], base_name = "sθ")
    s_vm = JuMP.@variable(model, [1:n_bus], base_name = "svm")
    s_p = JuMP.@variable(model, [1:n_edge], base_name = "sp")
    s_q = JuMP.@variable(model, [1:n_edge], base_name = "sq")

    mu = ws.mu[]
    @inbounds for k in 1:n_bus
        JuMP.add_to_expression!(qp_obj, ws.lambda_theta[k], s_theta[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_theta[k], s_theta[k])
        JuMP.add_to_expression!(qp_obj, ws.lambda_vm[k], s_vm[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_vm[k], s_vm[k])
    end
    @inbounds for k in 1:n_edge
        JuMP.add_to_expression!(qp_obj, ws.lambda_p[k], s_p[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_p[k], s_p[k])
        JuMP.add_to_expression!(qp_obj, ws.lambda_q[k], s_q[k])
        JuMP.add_to_expression!(qp_obj, 0.5 * mu, s_q[k], s_q[k])
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

    JuMP.optimize!(model)
    _check_qp_status!(model, "ALADIN-hybrid/AC coordinator QP")

    _extract_delta!(ws, dy_all)

    nu_theta = Vector{Float64}(undef, n_bus)
    nu_vm = Vector{Float64}(undef, n_bus)
    nu_p = Vector{Float64}(undef, n_edge)
    nu_q = Vector{Float64}(undef, n_edge)
    @inbounds for k in 1:n_bus
        nu_theta[k] = -JuMP.dual(cons_theta[k])
        nu_vm[k] = -JuMP.dual(cons_vm[k])
    end
    @inbounds for k in 1:n_edge
        nu_p[k] = -JuMP.dual(cons_p[k])
        nu_q[k] = -JuMP.dual(cons_q[k])
    end
    return (
        nu_theta = nu_theta,
        nu_vm = nu_vm,
        nu_p = nu_p,
        nu_q = nu_q,
    )
end
