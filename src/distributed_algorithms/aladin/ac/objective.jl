# ============================================================
# ALADIN-AC — Local (Phase 1) objective builder
# ============================================================
#
# The Phase-1 prox objective has three pieces:
#   1. Static quadratic    `(ρ/2)·Σ σ_k·v_k²`        (cached, refresh on ρ change)
#   2. Linear prox          `−ρ·Σ σ_k·x_prev_k·v_k`  (rebuilt every iter)
#   3. Coupling λᵀ·A_i·x   (rebuilt every iter)
#
# This file holds all three plus the per-variant `update_objective!`
# dispatch — the ALADIN analogue of `admm/ac/objective.jl`.

# ============================================================
# Quadratic (static, cached) part
# ============================================================
#
# `(ρ/2)·Σ σ_k·v_k²` is constant in (ρ, σ); cached on `zone.quad_static`
# at workspace build time and refreshed only when ρ changes
# (`adapt_rho_mu!`). AC has a single sigma-by-class variant:
# θ / vm / pg / qg / other.

# ------------------------------------------------------------
# Sigma-by-class
# ------------------------------------------------------------

function _build_quad_static_ac!(
    zone::AladinZoneDataAC,
    rho::Float64;
    sigma_theta::Float64,
    sigma_vm::Float64,
    sigma_pg::Float64,
    sigma_qg::Float64,
    sigma_other::Float64,
)
    class = zeros(UInt8, zone.n_vars)
    @inbounds for c in zone.theta_cols
        class[c] = 1
    end
    @inbounds for c in zone.vm_cols
        class[c] = 2
    end
    @inbounds for c in zone.pg_cols
        class[c] = 3
    end
    @inbounds for c in zone.qg_cols
        class[c] = 4
    end

    quad = zero(JuMP.QuadExpr)
    @inbounds for k in 1:zone.n_vars
        sig =
            class[k] == 1 ? sigma_theta :
            class[k] == 2 ? sigma_vm :
            class[k] == 3 ? sigma_pg :
            class[k] == 4 ? sigma_qg :
            sigma_other
        sig == 0.0 && continue
        v = zone.all_vars[k]
        JuMP.add_to_expression!(quad, 0.5 * rho * sig, v, v)
    end
    zone.quad_static = quad
    return nothing
end

# ------------------------------------------------------------
# Refresh helper (all zones)
# ------------------------------------------------------------

@inline function _refresh_quad_static_ac_all!(zones, rho, ws_or_cfg)
    @inbounds for zone in zones
        _build_quad_static_ac!(
            zone, rho;
            sigma_theta = ws_or_cfg.sigma_theta,
            sigma_vm = ws_or_cfg.sigma_vm,
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

@inline function _ac_prox_lin!(
    lin::JuMP.AffExpr,
    zone::AladinZoneDataAC,
    rho::Float64;
    sigma_theta::Float64,
    sigma_vm::Float64,
    sigma_pg::Float64,
    sigma_qg::Float64,
    sigma_other::Float64,
)
    class = zeros(UInt8, zone.n_vars)
    @inbounds for c in zone.theta_cols
        class[c] = 1
    end
    @inbounds for c in zone.vm_cols
        class[c] = 2
    end
    @inbounds for c in zone.pg_cols
        class[c] = 3
    end
    @inbounds for c in zone.qg_cols
        class[c] = 4
    end

    @inbounds for k in 1:zone.n_vars
        sig =
            class[k] == 1 ? sigma_theta :
            class[k] == 2 ? sigma_vm :
            class[k] == 3 ? sigma_pg :
            class[k] == 4 ? sigma_qg :
            sigma_other
        sig == 0.0 && continue
        JuMP.add_to_expression!(lin, -rho * sig * zone.x_prev[k], zone.all_vars[k])
    end
    return lin
end

# ============================================================
# Linear coupling terms (λᵀ · A_i · x)
# ============================================================

# ------------------------------------------------------------
# Bus
# ------------------------------------------------------------

@inline function _lambda_bus_ac!(
    lin::JuMP.AffExpr,
    ws::Union{ACStateWorkspaceALADIN, ACHybridWorkspaceALADIN},
    zone_pos::Int,
)
    zdata = ws.zones[zone_pos]
    for (k, pair) in enumerate(ws.bus_pairs)
        if pair.zone_i == zone_pos
            θ_col = pair.theta_col_i
            vm_col = pair.vm_col_i
            sign = 1.0
        elseif pair.zone_j == zone_pos
            θ_col = pair.theta_col_j
            vm_col = pair.vm_col_j
            sign = -1.0
        else
            continue
        end
        JuMP.add_to_expression!(lin, sign * ws.vars.theta.lambda[k], zdata.all_vars[θ_col])
        JuMP.add_to_expression!(lin, sign * ws.vars.v_mag.lambda[k], zdata.all_vars[vm_col])
    end
    return lin
end

# ------------------------------------------------------------
# Edge (P / Q flow)
# ------------------------------------------------------------

function _lambda_edge_pq_ac!(
    model::JuMP.Model,
    obj_expr,
    ws::Union{ACPowerWorkspaceALADIN, ACHybridWorkspaceALADIN},
    zone_pos::Int,
)
    part_pf = model[:p_flow_branches]
    part_qf = model[:q_flow_branches]
    zdata = ws.zones[zone_pos]
    branch_map =
        ws.partitions.parts[zone_pos].branch_idx_full_to_zone

    for (k, pair) in enumerate(ws.interzonal_line_pairs)
        if pair.zone_i == zone_pos
            le = pair.line_local_i
            sign = 1.0
        elseif pair.zone_j == zone_pos
            le = pair.line_local_j
            sign = -1.0
        else
            continue
        end
        line_global = pair.line_global
        l = branch_map[line_global]
        l == 0 && continue
        lp = sign * ws.vars.p_flow.lambda[k]
        lq = sign * ws.vars.q_flow.lambda[k]
        obj_expr += lp * part_pf[l, 1] + lq * part_qf[l, 1]
    end
    return obj_expr
end

# ------------------------------------------------------------
# Dummy branch (P / Q flow)
# ------------------------------------------------------------

function _lambda_dummy_pq_ac!(
    model::JuMP.Model,
    obj_expr,
    ws::Union{ACPowerWorkspaceALADIN, ACHybridWorkspaceALADIN},
    zone_pos::Int,
)
    part_pf = model[:p_flow_branches]
    part_qf = model[:q_flow_branches]
    for (k, pair) in enumerate(ws.dummy_line_pairs)
        if pair.zone_i == zone_pos
            l = pair.local_i
            sign = 1.0
        elseif pair.zone_j == zone_pos
            l = pair.local_j
            sign = -1.0
        else
            continue
        end
        lp = sign * ws.vars.p_dummy.lambda[k]
        lq = sign * ws.vars.q_dummy.lambda[k]
        obj_expr += lp * part_pf[l, 1] + lq * part_qf[l, 1]
    end
    return obj_expr
end

# ------------------------------------------------------------
# Split bus
# ------------------------------------------------------------

@inline function _lambda_split_bus_ac!(
    lin::JuMP.AffExpr,
    ws::ACHybridWorkspaceALADIN,
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
        # Sign-flipped consensus: both zones contribute +λ.
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

function update_objective!(part::PartitionAC, ws::ACStateWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_bus_ac!(lin, ws, zone_pos)
    _ac_prox_lin!(
        lin, zdata, ws.rho[];
        sigma_theta = ws.sigma_theta,
        sigma_vm = ws.sigma_vm,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )

    @objective(part.model, Min, part.obj_base + zdata.quad_static + lin)
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function update_objective!(part::PartitionAC, ws::ACPowerWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]
    model = part.model

    # AC only supports proximal_mode = :primal. The :pflow mode is rejected
    # at ALADINConfig construction time because p_flow_branches is non-linear
    # in (V, θ), so (p_flow − p_prev)² is not quadratic in primal variables.
    lin = zero(JuMP.AffExpr)
    _ac_prox_lin!(
        lin, zdata, ws.rho[];
        sigma_theta = ws.sigma_theta,
        sigma_vm = ws.sigma_vm,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )

    obj_expr = part.obj_base + zdata.quad_static + lin
    obj_expr = _lambda_edge_pq_ac!(model, obj_expr, ws, zone_pos)
    obj_expr = _lambda_dummy_pq_ac!(model, obj_expr, ws, zone_pos)
    @objective(model, Min, obj_expr)
    return nothing
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function update_objective!(part::PartitionAC, ws::ACHybridWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]
    model = part.model

    lin = zero(JuMP.AffExpr)
    _lambda_bus_ac!(lin, ws, zone_pos)
    _lambda_split_bus_ac!(lin, ws, zone_pos)
    _ac_prox_lin!(
        lin, zdata, ws.rho[];
        sigma_theta = ws.sigma_theta,
        sigma_vm = ws.sigma_vm,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )

    obj_expr = part.obj_base + zdata.quad_static + lin
    obj_expr = _lambda_edge_pq_ac!(model, obj_expr, ws, zone_pos)
    obj_expr = _lambda_dummy_pq_ac!(model, obj_expr, ws, zone_pos)
    @objective(model, Min, obj_expr)
    return nothing
end
