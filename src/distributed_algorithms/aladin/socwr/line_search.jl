# ============================================================
# ALADIN-SOCWR — Armijo line search on the HFD 2016 merit function
# ============================================================

const _SOCWRPhase2Workspace = Union{
    SOCWRStateWorkspaceALADIN,
    SOCWRPowerWorkspaceALADIN,
    SOCWRHybridWorkspaceALADIN,
}

@inline function _socwr_cost_at_trial(zone::AladinZoneDataSOCWR, α::Float64)
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

@inline function _socwr_flow_at_trial(
    zone::AladinZoneDataSOCWR,
    le::Int,
    α::Float64,
    which::Symbol,
)::Float64
    cols = zone.edge_flow_cols[le]
    coefs = which === :p ? zone.edge_pflow_coefs[le] : zone.edge_qflow_coefs[le]
    const_val = which === :p ? zone.edge_pflow_const[le] : zone.edge_qflow_const[le]
    acc = const_val
    y = zone.y
    delta = zone.delta
    @inbounds for i in eachindex(cols)
        c = cols[i]
        acc += coefs[i] * (y[c] + α * delta[c])
    end
    return acc
end

function _socwr_merit(ws::SOCWRStateWorkspaceALADIN, α::Float64)
    Φ = 0.0
    for zone in ws.zones
        Φ += _socwr_cost_at_trial(zone, α)
    end

    r_sq = 0.0
    λr = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.v2_col_i] + α * zi.delta[pair.v2_col_i]
        xj = zj.y[pair.v2_col_j] + α * zj.delta[pair.v2_col_j]
        d = xi - xj
        r_sq += d * d
        λr += ws.vars.v2.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.edge_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri_col = zi.edge_r_cols[pair.edge_local_i]
        rj_col = zj.edge_r_cols[pair.edge_local_j]
        ti_col = zi.edge_t_cols[pair.edge_local_i]
        tj_col = zj.edge_t_cols[pair.edge_local_j]
        dr = (zi.y[ri_col] + α * zi.delta[ri_col]) -
             (zj.y[rj_col] + α * zj.delta[rj_col])
        dt = (zi.y[ti_col] + α * zi.delta[ti_col]) -
             (zj.y[tj_col] + α * zj.delta[tj_col])
        r_sq += dr * dr + dt * dt
        λr += ws.vars.r.lambda[k] * dr + ws.vars.t.lambda[k] * dt
    end

    mu = ws.mu[]
    return Φ + λr + 0.5 * mu * r_sq
end

function _socwr_merit(ws::SOCWRPowerWorkspaceALADIN, α::Float64)
    Φ = 0.0
    for zone in ws.zones
        Φ += _socwr_cost_at_trial(zone, α)
    end

    r_sq = 0.0
    λr = 0.0
    @inbounds for (k, pair) in enumerate(ws.edge_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _socwr_flow_at_trial(zi, pair.edge_local_i, α, :p)
        pj_val = _socwr_flow_at_trial(zj, pair.edge_local_j, α, :p)
        qi_val = _socwr_flow_at_trial(zi, pair.edge_local_i, α, :q)
        qj_val = _socwr_flow_at_trial(zj, pair.edge_local_j, α, :q)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        λr += ws.vars.p_flow.lambda[k] * dp + ws.vars.q_flow.lambda[k] * dq
    end

    mu = ws.mu[]
    return Φ + λr + 0.5 * mu * r_sq
end

function _socwr_merit(ws::SOCWRHybridWorkspaceALADIN, α::Float64)
    Φ = 0.0
    for zone in ws.zones
        Φ += _socwr_cost_at_trial(zone, α)
    end

    r_sq = 0.0
    λr = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.v2_col_i] + α * zi.delta[pair.v2_col_i]
        xj = zj.y[pair.v2_col_j] + α * zj.delta[pair.v2_col_j]
        d = xi - xj
        r_sq += d * d
        λr += ws.vars.v2.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.edge_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri_col = zi.edge_r_cols[pair.edge_local_i]
        rj_col = zj.edge_r_cols[pair.edge_local_j]
        ti_col = zi.edge_t_cols[pair.edge_local_i]
        tj_col = zj.edge_t_cols[pair.edge_local_j]
        dr = (zi.y[ri_col] + α * zi.delta[ri_col]) -
             (zj.y[rj_col] + α * zj.delta[rj_col])
        dt = (zi.y[ti_col] + α * zi.delta[ti_col]) -
             (zj.y[tj_col] + α * zj.delta[tj_col])

        pi_val = _socwr_flow_at_trial(zi, pair.edge_local_i, α, :p)
        pj_val = _socwr_flow_at_trial(zj, pair.edge_local_j, α, :p)
        qi_val = _socwr_flow_at_trial(zi, pair.edge_local_i, α, :q)
        qj_val = _socwr_flow_at_trial(zj, pair.edge_local_j, α, :q)
        dp = pi_val - pj_val
        dq = qi_val - qj_val

        r_sq += dr * dr + dt * dt + dp * dp + dq * dq
        λr +=
            ws.vars.r.lambda[k] * dr + ws.vars.t.lambda[k] * dt +
            ws.vars.p_flow.lambda[k] * dp + ws.vars.q_flow.lambda[k] * dq
    end

    mu = ws.mu[]
    return Φ + λr + 0.5 * mu * r_sq
end

function _line_search_step!(
    ws::_SOCWRPhase2Workspace;
    c1::Float64 = 1e-4,
    ε_fd::Float64 = 1e-6,
    max_backtrack::Int = 25,
    α_min::Float64 = 1e-3,
)
    ws.line_search || return nothing

    Φ_0 = _socwr_merit(ws, 0.0)
    Φ_ε = _socwr_merit(ws, ε_fd)
    dΦ_0 = (Φ_ε - Φ_0) / ε_fd

    # If the QP direction is not descent, fall back to α_min to keep
    # the iteration making progress via the strong prox pull.
    if !isfinite(dΦ_0) || dΦ_0 >= 0.0
        ws.alpha2 = α_min
        return nothing
    end

    α = ws.alpha2_init
    @inbounds for _ in 1:max_backtrack
        Φ_α = _socwr_merit(ws, α)
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
