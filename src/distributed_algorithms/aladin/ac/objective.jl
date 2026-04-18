@inline function _ac_prox_quadexpr!(
    quad::JuMP.QuadExpr,
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

    for (k, pair) in enumerate(ws.edge_pairs)
        if pair.zone_i == zone_pos
            le = pair.edge_local_i
            sign = 1.0
        elseif pair.zone_j == zone_pos
            le = pair.edge_local_j
            sign = -1.0
        else
            continue
        end
        edge_global = pair.edge_global
        l = branch_map[edge_global]
        l == 0 && continue
        lp = sign * ws.vars.p_flow.lambda[k]
        lq = sign * ws.vars.q_flow.lambda[k]
        obj_expr += lp * part_pf[l, 1] + lq * part_qf[l, 1]
    end
    return obj_expr
end

# ------------------------------------------------------------
# DISPATCH
# ------------------------------------------------------------

function update_objective!(part::PartitionAC, ws::ACStateWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]

    lin = zero(JuMP.AffExpr)
    _lambda_bus_ac!(lin, ws, zone_pos)

    quad = zero(JuMP.QuadExpr)
    _ac_prox_quadexpr!(
        quad, zdata, ws.rho[];
        sigma_theta = ws.sigma_theta,
        sigma_vm = ws.sigma_vm,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )

    @objective(part.model, Min, part.obj_base + lin + quad)
    return nothing
end

function update_objective!(part::PartitionAC, ws::ACPowerWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]
    model = part.model

    quad = zero(JuMP.QuadExpr)
    _ac_prox_quadexpr!(
        quad, zdata, ws.rho[];
        sigma_theta = ws.sigma_theta,
        sigma_vm = ws.sigma_vm,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )

    obj_expr = part.obj_base + quad
    obj_expr = _lambda_edge_pq_ac!(model, obj_expr, ws, zone_pos)
    @objective(model, Min, obj_expr)
    return nothing
end

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

function update_objective!(part::PartitionAC, ws::ACHybridWorkspaceALADIN)
    zone_pos = part.position
    zdata = ws.zones[zone_pos]
    model = part.model

    lin = zero(JuMP.AffExpr)
    _lambda_bus_ac!(lin, ws, zone_pos)
    _lambda_split_bus_ac!(lin, ws, zone_pos)

    quad = zero(JuMP.QuadExpr)
    _ac_prox_quadexpr!(
        quad, zdata, ws.rho[];
        sigma_theta = ws.sigma_theta,
        sigma_vm = ws.sigma_vm,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )

    obj_expr = part.obj_base + lin + quad
    obj_expr = _lambda_edge_pq_ac!(model, obj_expr, ws, zone_pos)
    @objective(model, Min, obj_expr)
    return nothing
end
