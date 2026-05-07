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
# cos/sin of θ_f − θ_t), so Φ is evaluated point-wise from the analytic
# flow expressions. The directional derivative dΦ/dα|_{α=0} is computed
# **analytically** via the chain rule, matching the convention used in
# the SOCWR (`socwr/line_search.jl`) and DC (`dc/line_search.jl`)
# implementations — full canonical parity across the three models.

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

@inline function _ac_cost_dir_deriv_0(zone::AladinZoneDataAC)
    acc = 0.0
    y = zone.y
    delta = zone.delta
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        acc += (zone.pg_c2_sbase2[k] * y[col] + zone.pg_c1_sbase[k]) * delta[col]
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

# ------------------------------------------------------------
# Analytic directional derivatives dpf/dα and dqf/dα at α = 0
# ------------------------------------------------------------
#
# pf = g·τ²·vf² − τ·vf·vt·(g·cos(Δθ) + b·sin(Δθ))
# qf = −(b + b_sh/2)·τ²·vf² − τ·vf·vt·(−b·cos(Δθ) + g·sin(Δθ))
#
# Partials at the current iterate (α=0):
#   ∂pf/∂vf = 2g·τ²·vf − τ·vt·(g·cos(Δθ) + b·sin(Δθ))
#   ∂pf/∂vt = − τ·vf·(g·cos(Δθ) + b·sin(Δθ))
#   ∂pf/∂θf = τ·vf·vt·(g·sin(Δθ) − b·cos(Δθ))
#   ∂pf/∂θt = − ∂pf/∂θf
#   ∂qf/∂vf = −2(b + b_sh/2)·τ²·vf − τ·vt·(−b·cos(Δθ) + g·sin(Δθ))
#   ∂qf/∂vt = − τ·vf·(−b·cos(Δθ) + g·sin(Δθ))
#   ∂qf/∂θf = − τ·vf·vt·(b·sin(Δθ) + g·cos(Δθ))
#   ∂qf/∂θt = − ∂qf/∂θf
#
# Chain rule: dpf/dα|_{α=0} = Σ ∂pf/∂x · δx evaluated at y.

@inline function _ac_pf_dir_deriv(zone::AladinZoneDataAC, le::Int)
    fcol_θ = zone.edge_from_theta_col[le]
    tcol_θ = zone.edge_to_theta_col[le]
    fcol_v = zone.edge_from_vm_col[le]
    tcol_v = zone.edge_to_vm_col[le]
    θf = zone.y[fcol_θ]
    θt = zone.y[tcol_θ]
    vf = zone.y[fcol_v]
    vt = zone.y[tcol_v]
    δθf = zone.delta[fcol_θ]
    δθt = zone.delta[tcol_θ]
    δvf = zone.delta[fcol_v]
    δvt = zone.delta[tcol_v]
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    tau = zone.edge_inv_tap[le]
    Δθ = θf - θt - zone.edge_tap_angle[le]
    cΔ = cos(Δθ)
    sΔ = sin(Δθ)
    dpf_dvf = 2 * g * tau * tau * vf - tau * vt * (g * cΔ + b * sΔ)
    dpf_dvt = -tau * vf * (g * cΔ + b * sΔ)
    dpf_dθf = tau * vf * vt * (g * sΔ - b * cΔ)
    return dpf_dvf * δvf + dpf_dvt * δvt + dpf_dθf * (δθf - δθt)
end

@inline function _ac_qf_dir_deriv(zone::AladinZoneDataAC, le::Int)
    fcol_θ = zone.edge_from_theta_col[le]
    tcol_θ = zone.edge_to_theta_col[le]
    fcol_v = zone.edge_from_vm_col[le]
    tcol_v = zone.edge_to_vm_col[le]
    θf = zone.y[fcol_θ]
    θt = zone.y[tcol_θ]
    vf = zone.y[fcol_v]
    vt = zone.y[tcol_v]
    δθf = zone.delta[fcol_θ]
    δθt = zone.delta[tcol_θ]
    δvf = zone.delta[fcol_v]
    δvt = zone.delta[tcol_v]
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    b_sh = zone.edge_b_sh[le]
    tau = zone.edge_inv_tap[le]
    Δθ = θf - θt - zone.edge_tap_angle[le]
    cΔ = cos(Δθ)
    sΔ = sin(Δθ)
    dqf_dvf = -2 * (b + b_sh / 2) * tau * tau * vf - tau * vt * (-b * cΔ + g * sΔ)
    dqf_dvt = -tau * vf * (-b * cΔ + g * sΔ)
    dqf_dθf = -tau * vf * vt * (b * sΔ + g * cΔ)
    return dqf_dvf * δvf + dqf_dvt * δvt + dqf_dθf * (δθf - δθt)
end

@inline function _ac_dummy_pf_at_trial(
    pair::DummyLinePairAC, zone::AladinZoneDataAC, side::Symbol, α::Float64,
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
    pair::DummyLinePairAC, zone::AladinZoneDataAC, side::Symbol, α::Float64,
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

@inline function _ac_dummy_pf_dir_deriv(
    pair::DummyLinePairAC, zone::AladinZoneDataAC, side::Symbol,
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
    θf = zone.y[fc_θ]
    θt = zone.y[tc_θ]
    vf = zone.y[fc_v]
    vt = zone.y[tc_v]
    δθf = zone.delta[fc_θ]
    δθt = zone.delta[tc_θ]
    δvf = zone.delta[fc_v]
    δvt = zone.delta[tc_v]
    Δθ = θf - θt - ta
    cΔ = cos(Δθ)
    sΔ = sin(Δθ)
    dpf_dvf = 2 * g * tau * tau * vf - tau * vt * (g * cΔ + b * sΔ)
    dpf_dvt = -tau * vf * (g * cΔ + b * sΔ)
    dpf_dθf = tau * vf * vt * (g * sΔ - b * cΔ)
    return dpf_dvf * δvf + dpf_dvt * δvt + dpf_dθf * (δθf - δθt)
end

@inline function _ac_dummy_qf_dir_deriv(
    pair::DummyLinePairAC, zone::AladinZoneDataAC, side::Symbol,
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
    θf = zone.y[fc_θ]
    θt = zone.y[tc_θ]
    vf = zone.y[fc_v]
    vt = zone.y[tc_v]
    δθf = zone.delta[fc_θ]
    δθt = zone.delta[tc_θ]
    δvf = zone.delta[fc_v]
    δvt = zone.delta[tc_v]
    Δθ = θf - θt - ta
    cΔ = cos(Δθ)
    sΔ = sin(Δθ)
    dqf_dvf = -2 * (b + b_sh / 2) * tau * tau * vf - tau * vt * (-b * cΔ + g * sΔ)
    dqf_dvt = -tau * vf * (-b * cΔ + g * sΔ)
    dqf_dθf = -tau * vf * vt * (b * sΔ + g * cΔ)
    return dqf_dvf * δvf + dqf_dvt * δvt + dqf_dθf * (δθf - δθt)
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
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _ac_pf_at_trial(zi, pair.line_local_i, α)
        pj_val = _ac_pf_at_trial(zj, pair.line_local_j, α)
        qi_val = _ac_qf_at_trial(zi, pair.line_local_i, α)
        qj_val = _ac_qf_at_trial(zj, pair.line_local_j, α)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        λr += ws.vars.p_flow.lambda[k] * dp + ws.vars.q_flow.lambda[k] * dq
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
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
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _ac_pf_at_trial(zi, pair.line_local_i, α)
        pj_val = _ac_pf_at_trial(zj, pair.line_local_j, α)
        qi_val = _ac_qf_at_trial(zi, pair.line_local_i, α)
        qj_val = _ac_qf_at_trial(zj, pair.line_local_j, α)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        λr += ws.vars.p_flow.lambda[k] * dp + ws.vars.q_flow.lambda[k] * dq
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
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
# Analytical directional derivative dΦ/dα at α = 0
# ------------------------------------------------------------
#
# Φ(α) = Σ cost_i(y + αδ) + λᵀ·r(α) + (μ/2)·‖r(α)‖²
# dΦ/dα|_{α=0} = Σ ∇cost_i(y)·δ + (λ + μ·r₀)ᵀ·Δr
#
# where r₀ = r(0) is the current consensus residual and Δr = dr/dα|_0.
# For state pairs r is affine in α; for AC flows Δr comes from the
# analytic chain-rule helpers `_ac_pf_dir_deriv` / `_ac_qf_dir_deriv`.

function _ac_dΦ0(ws::ACStateWorkspaceALADIN)
    dΦ = 0.0
    for zone in ws.zones
        dΦ += _ac_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0_θ = zi.y[pair.theta_col_i] - zj.y[pair.theta_col_j]
        Δr_θ = zi.delta[pair.theta_col_i] - zj.delta[pair.theta_col_j]
        r0_v = zi.y[pair.vm_col_i] - zj.y[pair.vm_col_j]
        Δr_v = zi.delta[pair.vm_col_i] - zj.delta[pair.vm_col_j]
        dΦ += (ws.vars.theta.lambda[k] + mu * r0_θ) * Δr_θ
        dΦ += (ws.vars.v_mag.lambda[k] + mu * r0_v) * Δr_v
    end
    return dΦ
end

function _ac_dΦ0(ws::ACPowerWorkspaceALADIN)
    dΦ = 0.0
    for zone in ws.zones
        dΦ += _ac_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0_p =
            _ac_pf_at_trial(zi, pair.line_local_i, 0.0) -
            _ac_pf_at_trial(zj, pair.line_local_j, 0.0)
        Δr_p =
            _ac_pf_dir_deriv(zi, pair.line_local_i) -
            _ac_pf_dir_deriv(zj, pair.line_local_j)
        r0_q =
            _ac_qf_at_trial(zi, pair.line_local_i, 0.0) -
            _ac_qf_at_trial(zj, pair.line_local_j, 0.0)
        Δr_q =
            _ac_qf_dir_deriv(zi, pair.line_local_i) -
            _ac_qf_dir_deriv(zj, pair.line_local_j)
        dΦ += (ws.vars.p_flow.lambda[k] + mu * r0_p) * Δr_p
        dΦ += (ws.vars.q_flow.lambda[k] + mu * r0_q) * Δr_q
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0_p =
            _ac_dummy_pf_at_trial(pair, zi, :i, 0.0) -
            _ac_dummy_pf_at_trial(pair, zj, :j, 0.0)
        Δr_p = _ac_dummy_pf_dir_deriv(pair, zi, :i) -
               _ac_dummy_pf_dir_deriv(pair, zj, :j)
        r0_q =
            _ac_dummy_qf_at_trial(pair, zi, :i, 0.0) -
            _ac_dummy_qf_at_trial(pair, zj, :j, 0.0)
        Δr_q = _ac_dummy_qf_dir_deriv(pair, zi, :i) -
               _ac_dummy_qf_dir_deriv(pair, zj, :j)
        dΦ += (ws.vars.p_dummy.lambda[k] + mu * r0_p) * Δr_p
        dΦ += (ws.vars.q_dummy.lambda[k] + mu * r0_q) * Δr_q
    end
    return dΦ
end

function _ac_dΦ0(ws::ACHybridWorkspaceALADIN)
    dΦ = 0.0
    for zone in ws.zones
        dΦ += _ac_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0_θ = zi.y[pair.theta_col_i] - zj.y[pair.theta_col_j]
        Δr_θ = zi.delta[pair.theta_col_i] - zj.delta[pair.theta_col_j]
        r0_v = zi.y[pair.vm_col_i] - zj.y[pair.vm_col_j]
        Δr_v = zi.delta[pair.vm_col_i] - zj.delta[pair.vm_col_j]
        dΦ += (ws.vars.theta.lambda[k] + mu * r0_θ) * Δr_θ
        dΦ += (ws.vars.v_mag.lambda[k] + mu * r0_v) * Δr_v
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0_p =
            _ac_pf_at_trial(zi, pair.line_local_i, 0.0) -
            _ac_pf_at_trial(zj, pair.line_local_j, 0.0)
        Δr_p =
            _ac_pf_dir_deriv(zi, pair.line_local_i) -
            _ac_pf_dir_deriv(zj, pair.line_local_j)
        r0_q =
            _ac_qf_at_trial(zi, pair.line_local_i, 0.0) -
            _ac_qf_at_trial(zj, pair.line_local_j, 0.0)
        Δr_q =
            _ac_qf_dir_deriv(zi, pair.line_local_i) -
            _ac_qf_dir_deriv(zj, pair.line_local_j)
        dΦ += (ws.vars.p_flow.lambda[k] + mu * r0_p) * Δr_p
        dΦ += (ws.vars.q_flow.lambda[k] + mu * r0_q) * Δr_q
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0_p =
            _ac_dummy_pf_at_trial(pair, zi, :i, 0.0) -
            _ac_dummy_pf_at_trial(pair, zj, :j, 0.0)
        Δr_p = _ac_dummy_pf_dir_deriv(pair, zi, :i) -
               _ac_dummy_pf_dir_deriv(pair, zj, :j)
        r0_q =
            _ac_dummy_qf_at_trial(pair, zi, :i, 0.0) -
            _ac_dummy_qf_at_trial(pair, zj, :j, 0.0)
        Δr_q = _ac_dummy_qf_dir_deriv(pair, zi, :i) -
               _ac_dummy_qf_dir_deriv(pair, zj, :j)
        dΦ += (ws.vars.p_dummy.lambda[k] + mu * r0_p) * Δr_p
        dΦ += (ws.vars.q_dummy.lambda[k] + mu * r0_q) * Δr_q
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        # split-bus consensus is x_i + x_j = 0 (sign-flipped sum convention).
        r0_p = zi.y[pair.pinj_col_i] + zj.y[pair.pinj_col_j]
        Δr_p = zi.delta[pair.pinj_col_i] + zj.delta[pair.pinj_col_j]
        r0_q = zi.y[pair.qinj_col_i] + zj.y[pair.qinj_col_j]
        Δr_q = zi.delta[pair.qinj_col_i] + zj.delta[pair.qinj_col_j]
        dΦ += (ws.vars.p_inj.lambda[k] + mu * r0_p) * Δr_p
        dΦ += (ws.vars.q_inj.lambda[k] + mu * r0_q) * Δr_q
    end
    return dΦ
end

# ------------------------------------------------------------
# Armijo backtracking with analytic dΦ/dα|_{α=0}
# ------------------------------------------------------------

function _line_search_step!(
    ws::_ACPhase2Workspace;
    c1::Float64 = 1e-4,
    max_backtrack::Int = 25,
    α_min::Float64 = 1e-3,
)
    ws.line_search || return nothing

    Φ_0 = _ac_merit(ws, 0.0)
    dΦ_0 = _ac_dΦ0(ws)

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
