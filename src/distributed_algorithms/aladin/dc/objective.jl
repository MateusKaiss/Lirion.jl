# ============================================================
# ALADIN-DC — Local (Phase 1) objective builder
# ============================================================

# ------------------------------------------------------------
# Primal
# ------------------------------------------------------------

function _build_quad_static_primal_dc!(
    zone::AladinZoneData,
    rho::Float64;
    sigma_theta_shared::Float64,
    sigma_pg::Float64,
    sigma_other::Float64,
)
    is_shared_theta = falses(zone.n_vars)
    @inbounds for c in zone.shared_theta_cols
        is_shared_theta[c] = true
    end
    is_pg = falses(zone.n_vars)
    @inbounds for c in zone.pg_cols
        is_pg[c] = true
    end

    quad = zero(JuMP.QuadExpr)
    @inbounds for k in 1:zone.n_vars
        sig =
            is_shared_theta[k] ? sigma_theta_shared :
            is_pg[k] ? sigma_pg :
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

function _build_quad_static_pflow_dc!(
    zone::AladinZoneData,
    rho::Float64;
    sigma_pflow::Float64,
    sigma_pg::Float64,
    sigma_other::Float64,
)
    is_pg = falses(zone.n_vars)
    @inbounds for c in zone.pg_cols
        is_pg[c] = true
    end

    quad = zero(JuMP.QuadExpr)
    @inbounds for k in 1:zone.n_vars
        sig = is_pg[k] ? sigma_pg : sigma_other
        sig == 0.0 && continue
        v = zone.all_vars[k]
        JuMP.add_to_expression!(quad, 0.5 * rho * sig, v, v)
    end
    @inbounds for k in eachindex(zone.edge_from_cols)
        fcol = zone.edge_from_cols[k]
        tcol = zone.edge_to_cols[k]
        bval = zone.edge_b[k]
        coef_sq = 0.5 * rho * sigma_pflow * bval * bval
        v_f = zone.all_vars[fcol]
        v_t = zone.all_vars[tcol]
        JuMP.add_to_expression!(quad, coef_sq, v_f, v_f)
        JuMP.add_to_expression!(quad, -2 * coef_sq, v_f, v_t)
        JuMP.add_to_expression!(quad, coef_sq, v_t, v_t)
    end
    zone.quad_static = quad
    return nothing
end

# ------------------------------------------------------------
# Hybrid + pflow
# ------------------------------------------------------------
#
# Hybrid TieLine has consensus on BOTH shared θ AND p_flow on tie-lines.
# The `:pflow` proximal here keeps the primal prox on shared θ / pg /
# other AND adds the σ_pflow·b²·(θ_f − θ_t)² term so each consensus
# direction has its own prox-aligned regulariser.

function _build_quad_static_hybrid_pflow_dc!(
    zone::AladinZoneData,
    rho::Float64;
    sigma_theta_shared::Float64,
    sigma_pg::Float64,
    sigma_other::Float64,
    sigma_pflow::Float64,
)
    _build_quad_static_primal_dc!(
        zone, rho;
        sigma_theta_shared = sigma_theta_shared,
        sigma_pg = sigma_pg,
        sigma_other = sigma_other,
    )
    quad = zone.quad_static
    @inbounds for k in eachindex(zone.edge_from_cols)
        fcol = zone.edge_from_cols[k]
        tcol = zone.edge_to_cols[k]
        bval = zone.edge_b[k]
        coef_sq = 0.5 * rho * sigma_pflow * bval * bval
        v_f = zone.all_vars[fcol]
        v_t = zone.all_vars[tcol]
        JuMP.add_to_expression!(quad, coef_sq, v_f, v_f)
        JuMP.add_to_expression!(quad, -2 * coef_sq, v_f, v_t)
        JuMP.add_to_expression!(quad, coef_sq, v_t, v_t)
    end
    zone.quad_static = quad
    return nothing
end

# ------------------------------------------------------------
# Cache obj_static = obj_base + quad_static
# ------------------------------------------------------------

@inline function _set_obj_static_dc!(part::PartitionDC, zone::AladinZoneData)
    part.obj_static = part.obj_base + zone.quad_static
    return nothing
end

# ------------------------------------------------------------
# Refresh helpers (all zones)
# ------------------------------------------------------------

@inline function _refresh_quad_dc_primal!(parts, zones, rho, ws_or_cfg)
    @inbounds for (i, zone) in enumerate(zones)
        _build_quad_static_primal_dc!(
            zone, rho;
            sigma_theta_shared = ws_or_cfg.sigma_theta_shared,
            sigma_pg = ws_or_cfg.sigma_pg,
            sigma_other = ws_or_cfg.sigma_other,
        )
        _set_obj_static_dc!(parts[i], zone)
    end
    return nothing
end

@inline function _refresh_quad_dc_pflow!(parts, zones, rho, ws_or_cfg)
    @inbounds for (i, zone) in enumerate(zones)
        _build_quad_static_pflow_dc!(
            zone, rho;
            sigma_pflow = ws_or_cfg.sigma_pflow,
            sigma_pg = ws_or_cfg.sigma_pg,
            sigma_other = ws_or_cfg.sigma_other,
        )
        _set_obj_static_dc!(parts[i], zone)
    end
    return nothing
end

@inline function _refresh_quad_dc_hybrid_pflow!(parts, zones, rho, ws_or_cfg)
    @inbounds for (i, zone) in enumerate(zones)
        _build_quad_static_hybrid_pflow_dc!(
            zone, rho;
            sigma_theta_shared = ws_or_cfg.sigma_theta_shared,
            sigma_pg = ws_or_cfg.sigma_pg,
            sigma_other = ws_or_cfg.sigma_other,
            sigma_pflow = ws_or_cfg.sigma_pflow,
        )
        _set_obj_static_dc!(parts[i], zone)
    end
    return nothing
end

# ------------------------------------------------------------
# Public entry point: refresh_rho! (analogue of admm/refresh_rho!)
# ------------------------------------------------------------

refresh_rho!(ws::DCStateWorkspaceALADIN) =
    _refresh_quad_dc_primal!(ws.partitions.parts, ws.zones, ws.rho[], ws)

refresh_rho!(ws::DCPowerWorkspaceALADIN) =
    ws.proximal_mode === :pflow ?
        _refresh_quad_dc_pflow!(ws.partitions.parts, ws.zones, ws.rho[], ws) :
        _refresh_quad_dc_primal!(ws.partitions.parts, ws.zones, ws.rho[], ws)

refresh_rho!(ws::DCHybridWorkspaceALADIN) =
    ws.proximal_mode === :pflow ?
        _refresh_quad_dc_hybrid_pflow!(
            ws.partitions.parts, ws.zones, ws.rho_bus[], ws,
        ) :
        _refresh_quad_dc_primal!(ws.partitions.parts, ws.zones, ws.rho_bus[], ws)

# ============================================================
# Linear (per-iteration) prox part
# ============================================================

# ------------------------------------------------------------
# Primal
# ------------------------------------------------------------

@inline function _prox_lin_primal_dc!(
    lin::JuMP.AffExpr,
    zone::AladinZoneData,
    rho::Float64;
    sigma_theta_shared::Float64,
    sigma_pg::Float64,
    sigma_other::Float64,
)
    is_shared_theta = falses(zone.n_vars)
    @inbounds for c in zone.shared_theta_cols
        is_shared_theta[c] = true
    end
    is_pg = falses(zone.n_vars)
    @inbounds for c in zone.pg_cols
        is_pg[c] = true
    end

    @inbounds for k in 1:zone.n_vars
        sig =
            is_shared_theta[k] ? sigma_theta_shared :
            is_pg[k] ? sigma_pg :
            sigma_other
        sig == 0.0 && continue
        JuMP.add_to_expression!(lin, -rho * sig * zone.x_prev[k], zone.all_vars[k])
    end
    return lin
end

# ------------------------------------------------------------
# Pflow
# ------------------------------------------------------------

@inline function _prox_lin_pflow_dc!(
    lin::JuMP.AffExpr,
    zone::AladinZoneData,
    rho::Float64;
    sigma_pflow::Float64,
    sigma_pg::Float64,
    sigma_other::Float64,
)
    is_pg = falses(zone.n_vars)
    @inbounds for c in zone.pg_cols
        is_pg[c] = true
    end

    @inbounds for k in 1:zone.n_vars
        sig = is_pg[k] ? sigma_pg : sigma_other
        sig == 0.0 && continue
        JuMP.add_to_expression!(lin, -rho * sig * zone.x_prev[k], zone.all_vars[k])
    end

    @inbounds for k in eachindex(zone.edge_from_cols)
        fcol = zone.edge_from_cols[k]
        tcol = zone.edge_to_cols[k]
        bval = zone.edge_b[k]
        p_prev = bval * (zone.x_prev[fcol] - zone.x_prev[tcol])
        coef_lin = -rho * sigma_pflow * bval * p_prev
        JuMP.add_to_expression!(lin, coef_lin, zone.all_vars[fcol])
        JuMP.add_to_expression!(lin, -coef_lin, zone.all_vars[tcol])
    end
    return lin
end

# ------------------------------------------------------------
# Hybrid + pflow
# ------------------------------------------------------------
#
# Linear counterpart to `_build_quad_static_hybrid_pflow_dc!`: keeps
# the primal linear prox (shared θ / pg / other) and adds the edge
# pflow term so the per-iteration linearisation matches the cached
# quadratic.

@inline function _prox_lin_hybrid_pflow_dc!(
    lin::JuMP.AffExpr,
    zone::AladinZoneData,
    rho::Float64;
    sigma_theta_shared::Float64,
    sigma_pg::Float64,
    sigma_other::Float64,
    sigma_pflow::Float64,
)
    _prox_lin_primal_dc!(
        lin, zone, rho;
        sigma_theta_shared = sigma_theta_shared,
        sigma_pg = sigma_pg,
        sigma_other = sigma_other,
    )
    @inbounds for k in eachindex(zone.edge_from_cols)
        fcol = zone.edge_from_cols[k]
        tcol = zone.edge_to_cols[k]
        bval = zone.edge_b[k]
        p_prev = bval * (zone.x_prev[fcol] - zone.x_prev[tcol])
        coef_lin = -rho * sigma_pflow * bval * p_prev
        JuMP.add_to_expression!(lin, coef_lin, zone.all_vars[fcol])
        JuMP.add_to_expression!(lin, -coef_lin, zone.all_vars[tcol])
    end
    return lin
end

# ============================================================
# Linear coupling terms (λᵀ · A_i · x)
# ============================================================

# ------------------------------------------------------------
# Bus
# ------------------------------------------------------------

@inline function _lambda_bus_affexpr!(
    lin::JuMP.AffExpr,
    ws::DCStateWorkspaceALADIN,
    zone_pos::Int,
)
    zdata = ws.zones[zone_pos]
    for (k, pair) in enumerate(ws.bus_pairs)
        if pair.zone_i == zone_pos
            col = pair.col_i
            JuMP.add_to_expression!(lin, ws.vars.theta.lambda[k], zdata.all_vars[col])
        elseif pair.zone_j == zone_pos
            col = pair.col_j
            JuMP.add_to_expression!(lin, -ws.vars.theta.lambda[k], zdata.all_vars[col])
        end
    end
    return lin
end

# ------------------------------------------------------------
# Edge
# ------------------------------------------------------------

@inline function _lambda_edge_affexpr!(
    lin::JuMP.AffExpr,
    ws_zones::Vector{AladinZoneData},
    interzonal_line_pairs::Vector{InterzonalLinePair},
    lambda_edge::Vector{Float64},
    zone_pos::Int,
)
    zdata = ws_zones[zone_pos]
    for (k, pair) in enumerate(interzonal_line_pairs)
        if pair.zone_i == zone_pos
            le = pair.line_local_i
            sign = 1.0
        elseif pair.zone_j == zone_pos
            le = pair.line_local_j
            sign = -1.0
        else
            continue
        end
        fcol = zdata.edge_from_cols[le]
        tcol = zdata.edge_to_cols[le]
        bval = zdata.edge_b[le]
        coef = sign * lambda_edge[k] * bval
        JuMP.add_to_expression!(lin, coef, zdata.all_vars[fcol])
        JuMP.add_to_expression!(lin, -coef, zdata.all_vars[tcol])
    end
    return lin
end

# ------------------------------------------------------------
# Dummy branch
# ------------------------------------------------------------

@inline function _lambda_dummy_affexpr!(
    lin::JuMP.AffExpr,
    ws_zones::Vector{AladinZoneData},
    dummy_pairs::Vector{DummyLinePairDC},
    lambda_dummy::Vector{Float64},
    zone_pos::Int,
)
    for (k, pair) in enumerate(dummy_pairs)
        if pair.zone_i == zone_pos
            fcol, tcol, bval = pair.from_col_i, pair.to_col_i, pair.b_i
            sign = 1.0
        elseif pair.zone_j == zone_pos
            fcol, tcol, bval = pair.from_col_j, pair.to_col_j, pair.b_j
            sign = -1.0
        else
            continue
        end
        zdata = ws_zones[zone_pos]
        coef = sign * lambda_dummy[k] * bval
        JuMP.add_to_expression!(lin, coef, zdata.all_vars[fcol])
        JuMP.add_to_expression!(lin, -coef, zdata.all_vars[tcol])
    end
    return lin
end

# ============================================================
# Update objective
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function update_objective!(part::PartitionDC, ws::DCStateWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_bus_affexpr!(lin, ws, zone_pos)
    _prox_lin_primal_dc!(
        lin, zdata, ws.rho[];
        sigma_theta_shared = ws.sigma_theta_shared,
        sigma_pg = ws.sigma_pg,
        sigma_other = ws.sigma_other,
    )

    @objective(part.model, Min, part.obj_static + lin)
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function update_objective!(part::PartitionDC, ws::DCPowerWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_edge_affexpr!(
        lin,
        ws.zones,
        ws.interzonal_line_pairs,
        ws.vars.p_flow.lambda,
        zone_pos,
    )
    _lambda_dummy_affexpr!(
        lin, ws.zones, ws.dummy_line_pairs, ws.vars.p_dummy.lambda, zone_pos,
    )

    if ws.proximal_mode == :pflow
        _prox_lin_pflow_dc!(
            lin, zdata, ws.rho[];
            sigma_pflow = ws.sigma_pflow,
            sigma_pg = ws.sigma_pg,
            sigma_other = ws.sigma_other,
        )
    else
        _prox_lin_primal_dc!(
            lin, zdata, ws.rho[];
            sigma_theta_shared = ws.sigma_theta_shared,
            sigma_pg = ws.sigma_pg,
            sigma_other = ws.sigma_other,
        )
    end

    @objective(part.model, Min, part.obj_static + lin)
    return nothing
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function update_objective!(part::PartitionDC, ws::DCHybridWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    for (k, pair) in enumerate(ws.bus_pairs)
        if pair.zone_i == zone_pos
            col = pair.col_i
            JuMP.add_to_expression!(lin, ws.vars.theta.lambda[k], zdata.all_vars[col])
        elseif pair.zone_j == zone_pos
            col = pair.col_j
            JuMP.add_to_expression!(lin, -ws.vars.theta.lambda[k], zdata.all_vars[col])
        end
    end
    _lambda_edge_affexpr!(
        lin,
        ws.zones,
        ws.interzonal_line_pairs,
        ws.vars.p_flow.lambda,
        zone_pos,
    )
    _lambda_dummy_affexpr!(
        lin, ws.zones, ws.dummy_line_pairs, ws.vars.p_dummy.lambda, zone_pos,
    )

    for (k, pair) in enumerate(ws.split_bus_pairs)
        if pair.zone_i == zone_pos
            col = pair.p_inj_col_i
            JuMP.add_to_expression!(lin, ws.vars.p_inj.lambda[k], zdata.all_vars[col])
        elseif pair.zone_j == zone_pos
            col = pair.p_inj_col_j
            JuMP.add_to_expression!(lin, ws.vars.p_inj.lambda[k], zdata.all_vars[col])
        end
    end

    if ws.proximal_mode === :pflow
        _prox_lin_hybrid_pflow_dc!(
            lin, zdata, ws.rho_bus[];
            sigma_theta_shared = ws.sigma_theta_shared,
            sigma_pg = ws.sigma_pg,
            sigma_other = ws.sigma_other,
            sigma_pflow = ws.sigma_pflow,
        )
    else
        _prox_lin_primal_dc!(
            lin, zdata, ws.rho_bus[];
            sigma_theta_shared = ws.sigma_theta_shared,
            sigma_pg = ws.sigma_pg,
            sigma_other = ws.sigma_other,
        )
    end

    @objective(part.model, Min, part.obj_static + lin)
    return nothing
end
