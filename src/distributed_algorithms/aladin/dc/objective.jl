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

@inline function _prox_all_quadexpr!(
    quad::JuMP.QuadExpr,
    zone::AladinZoneData,
    rho::Float64;
    sigma_theta_shared::Float64 = 10.0,
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

@inline function _prox_pflow_dc!(
    quad::JuMP.QuadExpr,
    zone::AladinZoneData,
    rho::Float64;
    sigma_pflow::Float64 = 10.0,
    sigma_pg::Float64 = 1.0,
    sigma_other::Float64 = 1.0,
)
    is_pg = falses(zone.n_vars)
    @inbounds for c in zone.pg_cols
        is_pg[c] = true
    end

    @inbounds for k in 1:zone.n_vars
        sig = is_pg[k] ? sigma_pg : sigma_other
        v = zone.all_vars[k]
        x_prev = zone.x_prev[k]
        JuMP.add_to_expression!(quad, 0.5 * rho * sig, v, v)
        JuMP.add_to_expression!(quad, -rho * sig * x_prev, v)
    end

    @inbounds for k in eachindex(zone.edge_from_cols)
        fcol = zone.edge_from_cols[k]
        tcol = zone.edge_to_cols[k]
        bval = zone.edge_b[k]
        p_prev = bval * (zone.x_prev[fcol] - zone.x_prev[tcol])

        coef_sq = 0.5 * rho * sigma_pflow * bval * bval
        v_f = zone.all_vars[fcol]
        v_t = zone.all_vars[tcol]
        JuMP.add_to_expression!(quad, coef_sq, v_f, v_f)
        JuMP.add_to_expression!(quad, -2 * coef_sq, v_f, v_t)
        JuMP.add_to_expression!(quad, coef_sq, v_t, v_t)

        coef_lin = -rho * sigma_pflow * bval * p_prev
        JuMP.add_to_expression!(quad, coef_lin, v_f)
        JuMP.add_to_expression!(quad, -coef_lin, v_t)
    end
    return quad
end

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

@inline function _lambda_dummy_affexpr!(
    lin::JuMP.AffExpr,
    ws_zones::Vector{AladinZoneData},
    dummy_pairs::Vector{DummyBranchPairDC},
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

@inline function _lambda_edge_affexpr!(
    lin::JuMP.AffExpr,
    ws_zones::Vector{AladinZoneData},
    edge_pairs::Vector{SharedEdgePair},
    lambda_edge::Vector{Float64},
    zone_pos::Int,
)
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

function update_objective!(part::PartitionDC, ws::DCStateWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_bus_affexpr!(lin, ws, zone_pos)

    quad = zero(JuMP.QuadExpr)
    _prox_all_quadexpr!(
        quad, zdata, ws.rho[];
        sigma_theta_shared = ws.sigma_theta_shared,
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
    _lambda_edge_affexpr!(lin, ws.zones, ws.edge_pairs, ws.vars.p_flow.lambda, zone_pos)
    _lambda_dummy_affexpr!(
        lin, ws.zones, ws.dummy_branch_pairs, ws.vars.p_dummy.lambda, zone_pos,
    )

    quad = zero(JuMP.QuadExpr)
    if ws.proximal_mode == :pflow
        _prox_pflow_dc!(
            quad, zdata, ws.rho[];
            sigma_pflow = ws.sigma_pflow,
            sigma_pg = ws.sigma_pg,
            sigma_other = ws.sigma_other,
        )
    else
        _prox_all_quadexpr!(
            quad, zdata, ws.rho[];
            sigma_theta_shared = ws.sigma_theta_shared,
            sigma_pg = ws.sigma_pg,
            sigma_other = ws.sigma_other,
        )
    end

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
            JuMP.add_to_expression!(lin, ws.vars.theta.lambda[k], zdata.all_vars[col])
        elseif pair.zone_j == zone_pos
            col = pair.col_j
            JuMP.add_to_expression!(lin, -ws.vars.theta.lambda[k], zdata.all_vars[col])
        end
    end
    _lambda_edge_affexpr!(lin, ws.zones, ws.edge_pairs, ws.vars.p_flow.lambda, zone_pos)
    _lambda_dummy_affexpr!(
        lin, ws.zones, ws.dummy_branch_pairs, ws.vars.p_dummy.lambda, zone_pos,
    )

    for (k, pair) in enumerate(ws.split_bus_pairs)
        if pair.zone_i == zone_pos
            col = pair.pinj_col_i
            JuMP.add_to_expression!(lin, ws.vars.p_inj.lambda[k], zdata.all_vars[col])
        elseif pair.zone_j == zone_pos
            col = pair.pinj_col_j
            JuMP.add_to_expression!(lin, ws.vars.p_inj.lambda[k], zdata.all_vars[col])
        end
    end

    quad = zero(JuMP.QuadExpr)
    _prox_all_quadexpr!(
        quad, zdata, ws.rho_bus[];
        sigma_theta_shared = ws.sigma_theta_shared,
        sigma_pg = ws.sigma_pg,
        sigma_other = ws.sigma_other,
    )

    @objective(part.model, Min, part.obj_base + lin + quad)
    return nothing
end
