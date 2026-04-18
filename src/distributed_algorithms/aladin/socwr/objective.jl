# ============================================================
# ALADIN-SOCWR — Local (Phase 1) objective builder

const _SOCWRLocalWorkspace = Union{
    SOCWRStateWorkspaceALADIN,
    SOCWRPowerWorkspaceALADIN,
    SOCWRHybridWorkspaceALADIN,
}

@inline function _socwr_prox_quadexpr!(
    quad::JuMP.QuadExpr,
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
        v = zone.all_vars[k]
        x_prev = zone.x_prev[k]
        JuMP.add_to_expression!(quad, 0.5 * rho * sig, v, v)
        JuMP.add_to_expression!(quad, -rho * sig * x_prev, v)
    end
    return quad
end

# ------------------------------------------------------------
# Linear coupling terms (λᵀ · A_i · x)
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
# DISPATCH
# ------------------------------------------------------------

@inline function _socwr_prox_quad(zdata::AladinZoneDataSOCWR, ws::_SOCWRLocalWorkspace)
    quad = zero(JuMP.QuadExpr)
    _socwr_prox_quadexpr!(
        quad, zdata, ws.rho[];
        sigma_theta = ws.sigma_theta,
        sigma_v2 = ws.sigma_v2,
        sigma_r = ws.sigma_r,
        sigma_t = ws.sigma_t,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )
    return quad
end

function update_objective!(part::PartitionSOCWR, ws::SOCWRStateWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_bus_v2_socwr!(lin, ws, zone_pos)
    _lambda_edge_rt_socwr!(lin, ws, zone_pos)

    quad = _socwr_prox_quad(zdata, ws)

    @objective(part.model, Min, part.obj_base + lin + quad)
    return nothing
end

function update_objective!(part::PartitionSOCWR, ws::SOCWRPowerWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_edge_pq_socwr!(lin, ws, zone_pos)

    quad = _socwr_prox_quad(zdata, ws)

    @objective(part.model, Min, part.obj_base + lin + quad)
    return nothing
end

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

function update_objective!(part::PartitionSOCWR, ws::SOCWRHybridWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_bus_v2_socwr!(lin, ws, zone_pos)
    _lambda_edge_rt_socwr!(lin, ws, zone_pos)
    _lambda_edge_pq_socwr!(lin, ws, zone_pos)
    _lambda_split_bus_socwr!(lin, ws, zone_pos)

    quad = _socwr_prox_quad(zdata, ws)

    @objective(part.model, Min, part.obj_base + lin + quad)
    return nothing
end
