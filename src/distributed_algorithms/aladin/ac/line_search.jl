# ============================================================
# ALADIN-AC — Armijo line search on the HFD 2016 merit function
# ============================================================
#
# AC OPF is non-convex, so Theorem 4 of Houska-Frasch-Diehl 2016
# requires Algorithm 3 globalisation for local convergence guarantees.
# This file implements Armijo backtracking on the merit function:
#
#   Φ(α) = Σ cost_i(y + αδ) + λᵀ·r(α) + (μ/2)·‖r(α)‖²
#
# The flow residuals r(α) are non-linear in α (p_flow involves
# cos/sin of θ_f − θ_t), so Φ is evaluated point-wise via the analytic
# flow expressions in `flow_accessors.jl`. The directional derivative
# dΦ/dα|_{α=0} is computed by central finite difference with a small
# step h = 1e-6 — sufficient for the Armijo sufficient-decrease test.

const _ACPhase2Workspace = Union{
    ACStateWorkspaceALADIN,
    ACPowerWorkspaceALADIN,
    ACHybridWorkspaceALADIN,
}

# ------------------------------------------------------------
# Cost evaluation at trial point
# ------------------------------------------------------------

@inline function _ac_cost_at_trial(zone::AladinZoneDataAC, α::Float64)
    acc = 0.0
    y = zone.y
    delta = zone.delta
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        xk = y[col] + α * delta[col]
        acc += 0.5 * zone.pg_c2_sbase2[k] * xk * xk + zone.pg_c1_sbase[k] * xk
    end
    return acc
end

# ------------------------------------------------------------
# Trial-point flow evaluation (analytic, exact in α)
# ------------------------------------------------------------

@inline function _ac_pf_at_trial(zone::AladinZoneDataAC, le::Int, α::Float64)
    fcol_θ = zone.edge_from_theta_col[le]
    tcol_θ = zone.edge_to_theta_col[le]
    fcol_v = zone.edge_from_vm_col[le]
    tcol_v = zone.edge_to_vm_col[le]
    θf = zone.y[fcol_θ] + α * zone.delta[fcol_θ]
    θt = zone.y[tcol_θ] + α * zone.delta[tcol_θ]
    vf = zone.y[fcol_v] + α * zone.delta[fcol_v]
    vt = zone.y[tcol_v] + α * zone.delta[tcol_v]
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    tau = zone.edge_inv_tap[le]
    Δθ = θf - θt - zone.edge_tap_angle[le]
    return g * (tau * vf)^2 - tau * vf * vt * (g * cos(Δθ) + b * sin(Δθ))
end

@inline function _ac_qf_at_trial(zone::AladinZoneDataAC, le::Int, α::Float64)
    fcol_θ = zone.edge_from_theta_col[le]
    tcol_θ = zone.edge_to_theta_col[le]
    fcol_v = zone.edge_from_vm_col[le]
    tcol_v = zone.edge_to_vm_col[le]
    θf = zone.y[fcol_θ] + α * zone.delta[fcol_θ]
    θt = zone.y[tcol_θ] + α * zone.delta[tcol_θ]
    vf = zone.y[fcol_v] + α * zone.delta[fcol_v]
    vt = zone.y[tcol_v] + α * zone.delta[tcol_v]
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    b_sh = zone.edge_b_sh[le]
    tau = zone.edge_inv_tap[le]
    Δθ = θf - θt - zone.edge_tap_angle[le]
    return -(b + b_sh / 2) * (tau * vf)^2 -
           tau * vf * vt * (-b * cos(Δθ) + g * sin(Δθ))
end

@inline function _ac_dummy_pf_at_trial(
    pair::DummyBranchPairAC, zone::AladinZoneDataAC, side::Symbol, α::Float64,
)
    if side === :i
        fc_θ, tc_θ = pair.from_theta_col_i, pair.to_theta_col_i
        fc_v, tc_v = pair.from_vm_col_i, pair.to_vm_col_i
        g, b, tau, ta = pair.g_i, pair.b_i, pair.inv_tap_i, pair.tap_angle_i
    else
        fc_θ, tc_θ = pair.from_theta_col_j, pair.to_theta_col_j
        fc_v, tc_v = pair.from_vm_col_j, pair.to_vm_col_j
        g, b, tau, ta = pair.g_j, pair.b_j, pair.inv_tap_j, pair.tap_angle_j
    end
    θf = zone.y[fc_θ] + α * zone.delta[fc_θ]
    θt = zone.y[tc_θ] + α * zone.delta[tc_θ]
    vf = zone.y[fc_v] + α * zone.delta[fc_v]
    vt = zone.y[tc_v] + α * zone.delta[tc_v]
    Δθ = θf - θt - ta
    return g * (tau * vf)^2 - tau * vf * vt * (g * cos(Δθ) + b * sin(Δθ))
end

@inline function _ac_dummy_qf_at_trial(
    pair::DummyBranchPairAC, zone::AladinZoneDataAC, side::Symbol, α::Float64,
)
    if side === :i
        fc_θ, tc_θ = pair.from_theta_col_i, pair.to_theta_col_i
        fc_v, tc_v = pair.from_vm_col_i, pair.to_vm_col_i
        g, b, b_sh = pair.g_i, pair.b_i, pair.b_sh_i
        tau, ta = pair.inv_tap_i, pair.tap_angle_i
    else
        fc_θ, tc_θ = pair.from_theta_col_j, pair.to_theta_col_j
        fc_v, tc_v = pair.from_vm_col_j, pair.to_vm_col_j
        g, b, b_sh = pair.g_j, pair.b_j, pair.b_sh_j
        tau, ta = pair.inv_tap_j, pair.tap_angle_j
    end
    θf = zone.y[fc_θ] + α * zone.delta[fc_θ]
    θt = zone.y[tc_θ] + α * zone.delta[tc_θ]
    vf = zone.y[fc_v] + α * zone.delta[fc_v]
    vt = zone.y[tc_v] + α * zone.delta[tc_v]
    Δθ = θf - θt - ta
    return -(b + b_sh / 2) * (tau * vf)^2 -
           tau * vf * vt * (-b * cos(Δθ) + g * sin(Δθ))
end

# ------------------------------------------------------------
# Merit function Φ(α)
# ------------------------------------------------------------

function _ac_merit(ws::ACStateWorkspaceALADIN, α::Float64)
    Φ = 0.0
    for zone in ws.zones
        Φ += _ac_cost_at_trial(zone, α)
    end
    r_sq = 0.0
    λr = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        θ_i = zi.y[pair.theta_col_i] + α * zi.delta[pair.theta_col_i]
        θ_j = zj.y[pair.theta_col_j] + α * zj.delta[pair.theta_col_j]
        v_i = zi.y[pair.vm_col_i] + α * zi.delta[pair.vm_col_i]
        v_j = zj.y[pair.vm_col_j] + α * zj.delta[pair.vm_col_j]
        dθ = θ_i - θ_j
        dv = v_i - v_j
        r_sq += dθ * dθ + dv * dv
        λr += ws.vars.theta.lambda[k] * dθ + ws.vars.v_mag.lambda[k] * dv
    end
    return Φ + λr + 0.5 * ws.mu[] * r_sq
end

function _ac_merit(ws::ACPowerWorkspaceALADIN, α::Float64)
    Φ = 0.0
    for zone in ws.zones
        Φ += _ac_cost_at_trial(zone, α)
    end
    r_sq = 0.0
    λr = 0.0
    @inbounds for (k, pair) in enumerate(ws.edge_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _ac_pf_at_trial(zi, pair.edge_local_i, α)
        pj_val = _ac_pf_at_trial(zj, pair.edge_local_j, α)
        qi_val = _ac_qf_at_trial(zi, pair.edge_local_i, α)
        qj_val = _ac_qf_at_trial(zj, pair.edge_local_j, α)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        λr += ws.vars.p_flow.lambda[k] * dp + ws.vars.q_flow.lambda[k] * dq
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_branch_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _ac_dummy_pf_at_trial(pair, zi, :i, α)
        pj_val = _ac_dummy_pf_at_trial(pair, zj, :j, α)
        qi_val = _ac_dummy_qf_at_trial(pair, zi, :i, α)
        qj_val = _ac_dummy_qf_at_trial(pair, zj, :j, α)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        λr += ws.vars.p_dummy.lambda[k] * dp + ws.vars.q_dummy.lambda[k] * dq
    end
    return Φ + λr + 0.5 * ws.mu[] * r_sq
end

function _ac_merit(ws::ACHybridWorkspaceALADIN, α::Float64)
    Φ = 0.0
    for zone in ws.zones
        Φ += _ac_cost_at_trial(zone, α)
    end
    r_sq = 0.0
    λr = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        θ_i = zi.y[pair.theta_col_i] + α * zi.delta[pair.theta_col_i]
        θ_j = zj.y[pair.theta_col_j] + α * zj.delta[pair.theta_col_j]
        v_i = zi.y[pair.vm_col_i] + α * zi.delta[pair.vm_col_i]
        v_j = zj.y[pair.vm_col_j] + α * zj.delta[pair.vm_col_j]
        dθ = θ_i - θ_j
        dv = v_i - v_j
        r_sq += dθ * dθ + dv * dv
        λr += ws.vars.theta.lambda[k] * dθ + ws.vars.v_mag.lambda[k] * dv
    end
    @inbounds for (k, pair) in enumerate(ws.edge_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _ac_pf_at_trial(zi, pair.edge_local_i, α)
        pj_val = _ac_pf_at_trial(zj, pair.edge_local_j, α)
        qi_val = _ac_qf_at_trial(zi, pair.edge_local_i, α)
        qj_val = _ac_qf_at_trial(zj, pair.edge_local_j, α)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        λr += ws.vars.p_flow.lambda[k] * dp + ws.vars.q_flow.lambda[k] * dq
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_branch_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _ac_dummy_pf_at_trial(pair, zi, :i, α)
        pj_val = _ac_dummy_pf_at_trial(pair, zj, :j, α)
        qi_val = _ac_dummy_qf_at_trial(pair, zi, :i, α)
        qj_val = _ac_dummy_qf_at_trial(pair, zj, :j, α)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        λr += ws.vars.p_dummy.lambda[k] * dp + ws.vars.q_dummy.lambda[k] * dq
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi = zi.y[pair.pinj_col_i] + α * zi.delta[pair.pinj_col_i]
        pj = zj.y[pair.pinj_col_j] + α * zj.delta[pair.pinj_col_j]
        qi = zi.y[pair.qinj_col_i] + α * zi.delta[pair.qinj_col_i]
        qj = zj.y[pair.qinj_col_j] + α * zj.delta[pair.qinj_col_j]
        # split-bus consensus is x_i + x_j = 0 (sign-flipped sum convention).
        dp = pi + pj
        dq = qi + qj
        r_sq += dp * dp + dq * dq
        λr += ws.vars.p_inj.lambda[k] * dp + ws.vars.q_inj.lambda[k] * dq
    end
    return Φ + λr + 0.5 * ws.mu[] * r_sq
end

# ------------------------------------------------------------
# Armijo backtracking with finite-difference dΦ/dα|_{α=0}
# ------------------------------------------------------------

function _line_search_step!(
    ws::_ACPhase2Workspace;
    c1::Float64 = 1e-4,
    max_backtrack::Int = 25,
    α_min::Float64 = 1e-3,
    fd_step::Float64 = 1e-6,
)
    ws.line_search || return nothing

    Φ_0 = _ac_merit(ws, 0.0)
    # Central finite difference for dΦ/dα at α=0.
    Φ_fd = _ac_merit(ws, fd_step)
    dΦ_0 = (Φ_fd - Φ_0) / fd_step

    if !isfinite(dΦ_0) || dΦ_0 >= 0.0
        ws.alpha2 = α_min
        return nothing
    end

    α = ws.alpha2_init
    @inbounds for _ in 1:max_backtrack
        Φ_α = _ac_merit(ws, α)
        if isfinite(Φ_α) && Φ_α <= Φ_0 + c1 * α * dΦ_0
            ws.alpha2 = α
            return nothing
        end
        α *= 0.5
        if α < α_min
            break
        end
    end
    ws.alpha2 = α_min
    return nothing
end
