# ============================================================
# ALADIN-DC — Local (Phase 1) objective builder
# ------------------------------------------------------------
# Canonical HFD 2016 local subproblem:
#
#   min_{x_i}  f_i(x_i) + λᵀ (A_i x_i) + (ρ/2) ‖x_i - x_iᵏ‖²_{Σ}
#   s.t.      local equality + inequality constraints
#
# where:
#   * A_i x_i extracts the coupling variables (θ at shared buses for State,
#                                               p_flow_branches on interzonal
#                                               edges for Power, both for
#                                               Hybrid)
#   * x_iᵏ  is this zone's previous iterate (NOT a global consensus z)
#   * Σ    is the per-variable proximal weight (identity in our implementation)
#
# Expanded algebra:
#   (ρ/2) ‖x - x_prev‖² = (ρ/2) ‖x‖² - ρ x_prevᵀ x + const
#
# so the full objective becomes a base cost plus a per-variable affine/quad
# penalty. We rebuild this from scratch every iteration because x_prev and λ
# change.
# ============================================================

@inline function _prox_quadexpr_cols!(
    quad::JuMP.QuadExpr,
    zone::AladinZoneData,
    cols::AbstractVector{Int},
    rho::Float64,
)
    @inbounds for k in cols
        v = zone.all_vars[k]
        x_prev = zone.x_prev[k]
        JuMP.add_to_expression!(quad, 0.5 * rho, v, v)
        JuMP.add_to_expression!(quad, -rho * x_prev, v)
    end
    return quad
end

# Canonical HFD proximal penalty: weighted quadratic around x_prev on ALL
# local variables of the zone. Weights are differentiated per variable class
@inline function _prox_all_quadexpr!(
    quad::JuMP.QuadExpr,
    zone::AladinZoneData,
    rho::Float64;
    sigma_theta_shared::Float64 = 10.0,
    sigma_theta_other::Float64 = 10.0,
    sigma_pg::Float64 = 1.0,
    sigma_other::Float64 = 1.0,
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
        v = zone.all_vars[k]
        x_prev = zone.x_prev[k]
        JuMP.add_to_expression!(quad, 0.5 * rho * sig, v, v)
        JuMP.add_to_expression!(quad, -rho * sig * x_prev, v)
    end
    return quad
end

# Canonical HFD proximal is a weighted ‖x_i - x_iᵏ‖²_{Σ_i} on ALL local
# variables, independent of which sub-space is being coupled via A_i. So
# the proximal term is the same across State, Power and Hybrid — only the
# linear coupling term λᵀ(A_i x) differs. There is no extra edge-flow
# proximal; any penalisation on the branch-flow quantities is expressed
# through the linear λ term and through the consensus constraint in the
# Phase-2 QP.

@inline function _lambda_bus_affexpr!(
    lin::JuMP.AffExpr,
    ws::DCStateWorkspaceALADIN,
    zone_pos::Int,
)
    # ∑ λ_k · (sign_k · θ_i[bus]) for each shared-bus pair touching zone `zone_pos`.
    # We adopt the sign convention: + for zone_i, - for zone_j (zone_i < zone_j
    # because pairs are emitted in position order, see _build_shared_bus_pairs).
    zdata = ws.zones[zone_pos]
    for (k, pair) in enumerate(ws.bus_pairs)
        if pair.zone_i == zone_pos
            col = pair.col_i
            JuMP.add_to_expression!(lin, ws.lambda_bus[k], zdata.all_vars[col])
        elseif pair.zone_j == zone_pos
            col = pair.col_j
            JuMP.add_to_expression!(lin, -ws.lambda_bus[k], zdata.all_vars[col])
        end
    end
    return lin
end

@inline function _lambda_edge_affexpr!(
    lin::JuMP.AffExpr,
    ws_zones::Vector{AladinZoneData},
    edge_pairs::Vector{SharedEdgePair},
    lambda_edge::Vector{Float64},
    zone_pos::Int,
)
    # For each interzonal edge e touching this zone, add λ_e · sign · (b·(θ_from - θ_to))
    zdata = ws_zones[zone_pos]
    for (k, pair) in enumerate(edge_pairs)
        if pair.zone_i == zone_pos
            le = pair.edge_local_i
            sign = 1.0
        elseif pair.zone_j == zone_pos
            le = pair.edge_local_j
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

# ============================================================
# DISPATCH — update_objective!(part, ws) rebuilds the JuMP objective.
# ============================================================

function update_objective!(part::PartitionDC, ws::DCStateWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_bus_affexpr!(lin, ws, zone_pos)

    quad = zero(JuMP.QuadExpr)
    _prox_all_quadexpr!(
        quad, zdata, ws.rho[];
        sigma_theta_shared = ws.sigma_theta_shared,
        sigma_theta_other = ws.sigma_theta_other,
        sigma_pg = ws.sigma_pg,
        sigma_other = ws.sigma_other,
    )

    @objective(part.model, Min, part.obj_base + lin + quad)
    return nothing
end

function update_objective!(part::PartitionDC, ws::DCPowerWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_edge_affexpr!(lin, ws.zones, ws.edge_pairs, ws.lambda_edge, zone_pos)

    quad = zero(JuMP.QuadExpr)
    _prox_all_quadexpr!(
        quad, zdata, ws.rho[];
        sigma_theta_shared = ws.sigma_theta_shared,
        sigma_theta_other = ws.sigma_theta_other,
        sigma_pg = ws.sigma_pg,
        sigma_other = ws.sigma_other,
    )

    @objective(part.model, Min, part.obj_base + lin + quad)
    return nothing
end

function update_objective!(part::PartitionDC, ws::DCHybridWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    for (k, pair) in enumerate(ws.bus_pairs)
        if pair.zone_i == zone_pos
            col = pair.col_i
            JuMP.add_to_expression!(lin, ws.lambda_bus[k], zdata.all_vars[col])
        elseif pair.zone_j == zone_pos
            col = pair.col_j
            JuMP.add_to_expression!(lin, -ws.lambda_bus[k], zdata.all_vars[col])
        end
    end
    _lambda_edge_affexpr!(lin, ws.zones, ws.edge_pairs, ws.lambda_edge, zone_pos)

    # Single proximal on all local vars. For Hybrid we use rho_bus as the
    # scalar — the Σ weights handle variable-class scaling.
    quad = zero(JuMP.QuadExpr)
    _prox_all_quadexpr!(
        quad, zdata, ws.rho_bus[];
        sigma_theta_shared = ws.sigma_theta_shared,
        sigma_theta_other = ws.sigma_theta_other,
        sigma_pg = ws.sigma_pg,
        sigma_other = ws.sigma_other,
    )

    @objective(part.model, Min, part.obj_base + lin + quad)
    return nothing
end
