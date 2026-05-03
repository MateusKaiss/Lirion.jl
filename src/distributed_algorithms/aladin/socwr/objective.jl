# ============================================================
# ALADIN-SOCWR — Local (Phase 1) objective builder
# ============================================================
#
# The Phase-1 prox objective has three pieces:
#   1. Static quadratic    `(ρ/2)·Σ σ_k·v_k²`        (cached, refresh on ρ change)
#   2. Linear prox          `−ρ·Σ σ_k·x_prev_k·v_k`  (rebuilt every iter)
#   3. Coupling λᵀ·A_i·x   (rebuilt every iter)
#
# This file holds all three plus the per-variant `update_objective!`
# dispatch — the ALADIN analogue of `admm/socwr/objective.jl`.

const _SOCWRLocalWorkspace = Union{
    SOCWRStateWorkspaceALADIN,
    SOCWRPowerWorkspaceALADIN,
    SOCWRHybridWorkspaceALADIN,
}

# ============================================================
# Quadratic (static, cached) part
# ============================================================
#
# `(ρ/2)·Σ σ_k·v_k²` is constant in (ρ, σ); cached on `zone.quad_static`
# at workspace build time and refreshed only when ρ changes
# (`adapt_rho_mu!`). Two variants:
#   - `:primal` (canonical): sigmas by class — θ / v² / r / t / pg / qg
#   - `:pflow` (experimental): sigmas (pg, qg, other) on every var +
#     extra σ_pflow / σ_qflow on each interzonal edge's flow expression

# ------------------------------------------------------------
# Primal
# ------------------------------------------------------------

function _build_quad_static_primal_socwr!(
    zone::AladinZoneDataSOCWR,
    rho::Float64;
    sigma_theta::Float64,
    sigma_v2::Float64,
    sigma_r::Float64,
    sigma_t::Float64,
    sigma_pg::Float64,
    sigma_qg::Float64,
    sigma_other::Float64,
)
    class = zeros(UInt8, zone.n_vars)
    @inbounds for c in zone.theta_cols
        class[c] = 1
    end
    @inbounds for c in zone.v2_cols
        class[c] = 2
    end
    @inbounds for c in zone.r_cols
        class[c] = 3
    end
    @inbounds for c in zone.t_cols
        class[c] = 4
    end
    @inbounds for c in zone.pg_cols
        class[c] = 5
    end
    @inbounds for c in zone.qg_cols
        class[c] = 6
    end

    quad = zero(JuMP.QuadExpr)
    @inbounds for k in 1:zone.n_vars
        sig =
            class[k] == 1 ? sigma_theta :
            class[k] == 2 ? sigma_v2 :
            class[k] == 3 ? sigma_r :
            class[k] == 4 ? sigma_t :
            class[k] == 5 ? sigma_pg :
            class[k] == 6 ? sigma_qg :
            sigma_other
        sig == 0.0 && continue
        v = zone.all_vars[k]
        JuMP.add_to_expression!(quad, 0.5 * rho * sig, v, v)
    end
    zone.quad_static = quad
    return nothing
end

# ------------------------------------------------------------
# Pflow
# ------------------------------------------------------------

function _build_quad_static_pflow_socwr!(
    zone::AladinZoneDataSOCWR,
    edge_pairs::Vector{SharedEdgePairSOCWR},
    zone_pos::Int,
    rho::Float64;
    sigma_pflow::Float64,
    sigma_qflow::Float64,
    sigma_pg::Float64,
    sigma_qg::Float64,
    sigma_other::Float64,
)
    is_pg = falses(zone.n_vars)
    @inbounds for c in zone.pg_cols
        is_pg[c] = true
    end
    is_qg = falses(zone.n_vars)
    @inbounds for c in zone.qg_cols
        is_qg[c] = true
    end

    quad = zero(JuMP.QuadExpr)
    @inbounds for k in 1:zone.n_vars
        sig = is_pg[k] ? sigma_pg : (is_qg[k] ? sigma_qg : sigma_other)
        sig == 0.0 && continue
        v = zone.all_vars[k]
        JuMP.add_to_expression!(quad, 0.5 * rho * sig, v, v)
    end

    @inbounds for pair in edge_pairs
        if pair.zone_i == zone_pos
            le = pair.edge_local_i
        elseif pair.zone_j == zone_pos
            le = pair.edge_local_j
        else
            continue
        end
        cols = zone.edge_flow_cols[le]
        cp = zone.edge_pflow_coefs[le]
        cq = zone.edge_qflow_coefs[le]

        for i in eachindex(cols)
            ci_p = cp[i]
            ci_q = cq[i]
            v_i = zone.all_vars[cols[i]]
            if ci_p != 0.0
                JuMP.add_to_expression!(
                    quad, 0.5 * rho * sigma_pflow * ci_p * ci_p, v_i, v_i,
                )
            end
            if ci_q != 0.0
                JuMP.add_to_expression!(
                    quad, 0.5 * rho * sigma_qflow * ci_q * ci_q, v_i, v_i,
                )
            end
            for j in (i+1):length(cols)
                cj_p = cp[j]
                cj_q = cq[j]
                v_j = zone.all_vars[cols[j]]
                if ci_p != 0.0 && cj_p != 0.0
                    JuMP.add_to_expression!(
                        quad, rho * sigma_pflow * ci_p * cj_p, v_i, v_j,
                    )
                end
                if ci_q != 0.0 && cj_q != 0.0
                    JuMP.add_to_expression!(
                        quad, rho * sigma_qflow * ci_q * cj_q, v_i, v_j,
                    )
                end
            end
        end
    end
    zone.quad_static = quad
    return nothing
end

# ------------------------------------------------------------
# Refresh helpers (all zones)
# ------------------------------------------------------------

@inline function _refresh_quad_static_socwr_primal_all!(zones, rho, ws_or_cfg)
    @inbounds for zone in zones
        _build_quad_static_primal_socwr!(
            zone, rho;
            sigma_theta = ws_or_cfg.sigma_theta,
            sigma_v2 = ws_or_cfg.sigma_v2,
            sigma_r = ws_or_cfg.sigma_r,
            sigma_t = ws_or_cfg.sigma_t,
            sigma_pg = ws_or_cfg.sigma_pg,
            sigma_qg = ws_or_cfg.sigma_qg,
            sigma_other = ws_or_cfg.sigma_other,
        )
    end
    return nothing
end

@inline function _refresh_quad_static_socwr_pflow_all!(zones, edge_pairs, rho, ws_or_cfg)
    @inbounds for (zone_pos, zone) in enumerate(zones)
        _build_quad_static_pflow_socwr!(
            zone, edge_pairs, zone_pos, rho;
            sigma_pflow = ws_or_cfg.sigma_pflow,
            sigma_qflow = ws_or_cfg.sigma_qflow,
            sigma_pg = ws_or_cfg.sigma_pg,
            sigma_qg = ws_or_cfg.sigma_qg,
            sigma_other = ws_or_cfg.sigma_other,
        )
    end
    return nothing
end

# ============================================================
# Linear (per-iteration) prox part
# ============================================================

# ------------------------------------------------------------
# Primal
# ------------------------------------------------------------

@inline function _socwr_prox_lin!(
    lin::JuMP.AffExpr,
    zone::AladinZoneDataSOCWR,
    rho::Float64;
    sigma_theta::Float64,
    sigma_v2::Float64,
    sigma_r::Float64,
    sigma_t::Float64,
    sigma_pg::Float64,
    sigma_qg::Float64,
    sigma_other::Float64,
)
    # 0 = other, 1 = θ, 2 = v², 3 = R, 4 = T, 5 = pg, 6 = qg
    class = zeros(UInt8, zone.n_vars)
    @inbounds for c in zone.theta_cols
        class[c] = 1
    end
    @inbounds for c in zone.v2_cols
        class[c] = 2
    end
    @inbounds for c in zone.r_cols
        class[c] = 3
    end
    @inbounds for c in zone.t_cols
        class[c] = 4
    end
    @inbounds for c in zone.pg_cols
        class[c] = 5
    end
    @inbounds for c in zone.qg_cols
        class[c] = 6
    end

    @inbounds for k in 1:zone.n_vars
        sig =
            class[k] == 1 ? sigma_theta :
            class[k] == 2 ? sigma_v2 :
            class[k] == 3 ? sigma_r :
            class[k] == 4 ? sigma_t :
            class[k] == 5 ? sigma_pg :
            class[k] == 6 ? sigma_qg :
            sigma_other
        sig == 0.0 && continue
        JuMP.add_to_expression!(lin, -rho * sig * zone.x_prev[k], zone.all_vars[k])
    end
    return lin
end

# ------------------------------------------------------------
# Pflow
# ------------------------------------------------------------

@inline function _prox_pflow_lin_socwr!(
    lin::JuMP.AffExpr,
    zone::AladinZoneDataSOCWR,
    edge_pairs::Vector{SharedEdgePairSOCWR},
    zone_pos::Int,
    rho::Float64;
    sigma_pflow::Float64,
    sigma_qflow::Float64,
    sigma_pg::Float64,
    sigma_qg::Float64,
    sigma_other::Float64,
)
    is_pg = falses(zone.n_vars)
    @inbounds for c in zone.pg_cols
        is_pg[c] = true
    end
    is_qg = falses(zone.n_vars)
    @inbounds for c in zone.qg_cols
        is_qg[c] = true
    end

    @inbounds for k in 1:zone.n_vars
        sig = is_pg[k] ? sigma_pg : (is_qg[k] ? sigma_qg : sigma_other)
        sig == 0.0 && continue
        JuMP.add_to_expression!(lin, -rho * sig * zone.x_prev[k], zone.all_vars[k])
    end

    @inbounds for pair in edge_pairs
        if pair.zone_i == zone_pos
            le = pair.edge_local_i
        elseif pair.zone_j == zone_pos
            le = pair.edge_local_j
        else
            continue
        end
        cols = zone.edge_flow_cols[le]
        cp = zone.edge_pflow_coefs[le]
        cq = zone.edge_qflow_coefs[le]

        D_p = 0.0
        D_q = 0.0
        for i in eachindex(cols)
            D_p += cp[i] * zone.x_prev[cols[i]]
            D_q += cq[i] * zone.x_prev[cols[i]]
        end

        for i in eachindex(cols)
            ci_p = cp[i]
            ci_q = cq[i]
            v_i = zone.all_vars[cols[i]]
            if ci_p != 0.0
                JuMP.add_to_expression!(lin, -rho * sigma_pflow * D_p * ci_p, v_i)
            end
            if ci_q != 0.0
                JuMP.add_to_expression!(lin, -rho * sigma_qflow * D_q * ci_q, v_i)
            end
        end
    end
    return lin
end

# ------------------------------------------------------------
# Primal — full-arg helper
# ------------------------------------------------------------

@inline function _socwr_prox_lin_full!(lin::JuMP.AffExpr, zdata::AladinZoneDataSOCWR, ws::_SOCWRLocalWorkspace)
    _socwr_prox_lin!(
        lin, zdata, ws.rho[];
        sigma_theta = ws.sigma_theta,
        sigma_v2 = ws.sigma_v2,
        sigma_r = ws.sigma_r,
        sigma_t = ws.sigma_t,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )
    return lin
end

# ============================================================
# Linear coupling terms (λᵀ · A_i · x)
# ============================================================

# ------------------------------------------------------------
# Bus (v²)
# ------------------------------------------------------------

@inline function _lambda_bus_v2_socwr!(
    lin::JuMP.AffExpr,
    ws::Union{SOCWRStateWorkspaceALADIN, SOCWRHybridWorkspaceALADIN},
    zone_pos::Int,
)
    zdata = ws.zones[zone_pos]
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        if pair.zone_i == zone_pos
            col = pair.v2_col_i
            sign = 1.0
        elseif pair.zone_j == zone_pos
            col = pair.v2_col_j
            sign = -1.0
        else
            continue
        end
        JuMP.add_to_expression!(lin, sign * ws.vars.v2.lambda[k], zdata.all_vars[col])
    end
    return lin
end

# ------------------------------------------------------------
# Edge (R, T)
# ------------------------------------------------------------

@inline function _lambda_edge_rt_socwr!(
    lin::JuMP.AffExpr,
    ws::Union{SOCWRStateWorkspaceALADIN, SOCWRHybridWorkspaceALADIN},
    zone_pos::Int,
)
    zdata = ws.zones[zone_pos]
    @inbounds for (k, pair) in enumerate(ws.edge_pairs)
        if pair.zone_i == zone_pos
            le = pair.edge_local_i
            sign = 1.0
        elseif pair.zone_j == zone_pos
            le = pair.edge_local_j
            sign = -1.0
        else
            continue
        end
        r_col = zdata.edge_r_cols[le]
        t_col = zdata.edge_t_cols[le]
        JuMP.add_to_expression!(lin, sign * ws.vars.r.lambda[k], zdata.all_vars[r_col])
        JuMP.add_to_expression!(lin, sign * ws.vars.t.lambda[k], zdata.all_vars[t_col])
    end
    return lin
end

# ------------------------------------------------------------
# Edge (P / Q flow)
# ------------------------------------------------------------

# Coupling term on p_flow / q_flow: the flows are affine expressions
# in the zone variables (SOCWR Jabr), so the coupling collapses to a
# plain AffExpr via the cached expansion.
@inline function _lambda_edge_pq_socwr!(
    lin::JuMP.AffExpr,
    ws::Union{SOCWRPowerWorkspaceALADIN, SOCWRHybridWorkspaceALADIN},
    zone_pos::Int,
)
    zdata = ws.zones[zone_pos]
    @inbounds for (k, pair) in enumerate(ws.edge_pairs)
        if pair.zone_i == zone_pos
            le = pair.edge_local_i
            sign = 1.0
        elseif pair.zone_j == zone_pos
            le = pair.edge_local_j
            sign = -1.0
        else
            continue
        end
        lp = sign * ws.vars.p_flow.lambda[k]
        lq = sign * ws.vars.q_flow.lambda[k]
        cols = zdata.edge_flow_cols[le]
        cp = zdata.edge_pflow_coefs[le]
        cq = zdata.edge_qflow_coefs[le]
        @inbounds for i in eachindex(cols)
            JuMP.add_to_expression!(lin, lp * cp[i], zdata.all_vars[cols[i]])
            JuMP.add_to_expression!(lin, lq * cq[i], zdata.all_vars[cols[i]])
        end
        # The constant offset in the flow expansion does not affect
        # the gradient in x but shifts the Lagrangian value by a
        # constant; we add it to keep the objective consistent when
        # users read part.generation_cost.
        lin.constant += lp * zdata.edge_pflow_const[le]
        lin.constant += lq * zdata.edge_qflow_const[le]
    end
    return lin
end

# ------------------------------------------------------------
# Split bus
# ------------------------------------------------------------

@inline function _lambda_split_bus_socwr!(
    lin::JuMP.AffExpr,
    ws::SOCWRHybridWorkspaceALADIN,
    zone_pos::Int,
)
    zdata = ws.zones[zone_pos]
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        if pair.zone_i == zone_pos
            pcol = pair.pinj_col_i
            qcol = pair.qinj_col_i
        elseif pair.zone_j == zone_pos
            pcol = pair.pinj_col_j
            qcol = pair.qinj_col_j
        else
            continue
        end
        # Sign-flipped consensus: both zones contribute +λ (no sign switch).
        JuMP.add_to_expression!(lin, ws.vars.p_inj.lambda[k], zdata.all_vars[pcol])
        JuMP.add_to_expression!(lin, ws.vars.q_inj.lambda[k], zdata.all_vars[qcol])
    end
    return lin
end

# ============================================================
# Update objective
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function update_objective!(part::PartitionSOCWR, ws::SOCWRStateWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_bus_v2_socwr!(lin, ws, zone_pos)
    _lambda_edge_rt_socwr!(lin, ws, zone_pos)
    _socwr_prox_lin_full!(lin, zdata, ws)

    @objective(part.model, Min, part.obj_base + zdata.quad_static + lin)
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function update_objective!(part::PartitionSOCWR, ws::SOCWRPowerWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_edge_pq_socwr!(lin, ws, zone_pos)

    if ws.proximal_mode == :pflow
        _prox_pflow_lin_socwr!(
            lin, zdata, ws.edge_pairs, zone_pos, ws.rho[];
            sigma_pflow = ws.sigma_pflow,
            sigma_qflow = ws.sigma_qflow,
            sigma_pg = ws.sigma_pg,
            sigma_qg = ws.sigma_qg,
            sigma_other = ws.sigma_other,
        )
    else
        _socwr_prox_lin_full!(lin, zdata, ws)
    end

    @objective(part.model, Min, part.obj_base + zdata.quad_static + lin)
    return nothing
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function update_objective!(part::PartitionSOCWR, ws::SOCWRHybridWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_bus_v2_socwr!(lin, ws, zone_pos)
    _lambda_edge_rt_socwr!(lin, ws, zone_pos)
    _lambda_edge_pq_socwr!(lin, ws, zone_pos)
    _lambda_split_bus_socwr!(lin, ws, zone_pos)
    _socwr_prox_lin_full!(lin, zdata, ws)

    @objective(part.model, Min, part.obj_base + zdata.quad_static + lin)
    return nothing
end
