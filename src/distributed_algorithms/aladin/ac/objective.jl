# ============================================================
# ALADIN-AC — Local (Phase 1) objective builder
# ============================================================

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
# Pflow — raw-var part only
# ------------------------------------------------------------

function _build_quad_static_ac_pflow!(
    zone::AladinZoneDataAC,
    rho::Float64;
    sigma_pg::Float64,
    sigma_qg::Float64,
    sigma_other::Float64,
)
    class = zeros(UInt8, zone.n_vars)
    @inbounds for c in zone.pg_cols
        class[c] = 3
    end
    @inbounds for c in zone.qg_cols
        class[c] = 4
    end

    quad = zero(JuMP.QuadExpr)
    @inbounds for k in 1:zone.n_vars
        sig =
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
# Refresh helpers (all zones)
# ------------------------------------------------------------

@inline function _refresh_quad_ac!(zones, rho, ws_or_cfg)
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

@inline function _refresh_quad_ac_pflow!(zones, rho, ws_or_cfg)
    @inbounds for zone in zones
        _build_quad_static_ac_pflow!(
            zone, rho;
            sigma_pg = ws_or_cfg.sigma_pg,
            sigma_qg = ws_or_cfg.sigma_qg,
            sigma_other = ws_or_cfg.sigma_other,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# Public entry point: refresh_rho!
# ------------------------------------------------------------

refresh_rho!(ws::ACStateWorkspaceALADIN) = _refresh_quad_ac!(ws.zones, ws.rho[], ws)

refresh_rho!(ws::ACPowerWorkspaceALADIN) =
    ws.proximal_mode === :pflow ?
    _refresh_quad_ac_pflow!(ws.zones, ws.rho[], ws) :
    _refresh_quad_ac!(ws.zones, ws.rho[], ws)

refresh_rho!(ws::ACHybridWorkspaceALADIN) =
    ws.proximal_mode === :pflow ?
    _refresh_quad_ac_pflow!(ws.zones, ws.rho[], ws) :
    _refresh_quad_ac!(ws.zones, ws.rho[], ws)

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
    return nothing
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
    zone_data = ws.zones[zone_pos]
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
        JuMP.add_to_expression!(lin, sign * ws.vars.theta.lambda[k], zone_data.all_vars[θ_col])
        JuMP.add_to_expression!(lin, sign * ws.vars.v_mag.lambda[k], zone_data.all_vars[vm_col])
    end
    return nothing
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
    zone_data = ws.zones[zone_pos]
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
    part_prev = model[:p_reverse_flow_branches]
    part_qrev = model[:q_reverse_flow_branches]
    for (k, pair) in enumerate(ws.dummy_line_pairs)
        if pair.zone_i == zone_pos
            l = pair.local_i
            lp = ws.vars.p_dummy.lambda[k]
            lq = ws.vars.q_dummy.lambda[k]
            obj_expr += lp * (-part_prev[l, 1]) + lq * (-part_qrev[l, 1])
        elseif pair.zone_j == zone_pos
            l = pair.local_j
            lp = -ws.vars.p_dummy.lambda[k]
            lq = -ws.vars.q_dummy.lambda[k]
            obj_expr += lp * part_pf[l, 1] + lq * part_qf[l, 1]
        end
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
    zone_data = ws.zones[zone_pos]
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        c = findfirst(==(zone_pos), pair.zones)
        c === nothing && continue
        pcol = pair.p_inj_cols[c]
        qcol = pair.q_inj_cols[c]
        JuMP.add_to_expression!(lin, ws.vars.p_inj.lambda[k], zone_data.all_vars[pcol])
        JuMP.add_to_expression!(lin, ws.vars.q_inj.lambda[k], zone_data.all_vars[qcol])
    end
    return nothing
end

# ============================================================
# Flow proximal (NonlinearExpr, :pflow mode)
# ============================================================


@inline function _ac_flow_at(
    g::Float64, b::Float64, b_sh::Float64,
    inv_tap::Float64, tap_angle::Float64,
    θ_f::Float64, θ_t::Float64, vm_f::Float64, vm_t::Float64,
)
    Δθ = θ_f - θ_t - tap_angle
    cos_Δθ = cos(Δθ)
    sin_Δθ = sin(Δθ)
    inv_tap_vm_f = inv_tap * vm_f
    inv_tap_vm_f_sq = inv_tap_vm_f * inv_tap_vm_f
    cross = inv_tap_vm_f * vm_t

    p = g * inv_tap_vm_f_sq - cross * (g * cos_Δθ + b * sin_Δθ)
    q = -(b + 0.5 * b_sh) * inv_tap_vm_f_sq -
        cross * (-b * cos_Δθ + g * sin_Δθ)
    return (p, q)
end

@inline function _ac_flow_at_to(
    g::Float64, b::Float64, b_sh::Float64,
    inv_tap::Float64, tap_angle::Float64,
    θ_f::Float64, θ_t::Float64, vm_f::Float64, vm_t::Float64,
)
    Δθ = θ_f - θ_t - tap_angle
    cos_Δθ = cos(Δθ)
    sin_Δθ = sin(Δθ)
    cross = inv_tap * vm_f * vm_t
    vm_t_sq = vm_t * vm_t

    p = -g * vm_t_sq + cross * (g * cos_Δθ - b * sin_Δθ)
    q = (b + 0.5 * b_sh) * vm_t_sq + cross * (-b * cos_Δθ - g * sin_Δθ)
    return (p, q)
end

# ------------------------------------------------------------
# Previous-iterate flow on an interzonal edge
# ------------------------------------------------------------

@inline function _ac_flow_prev_interzonal(zone::AladinZoneDataAC, le::Int)
    θ_f = zone.x_prev[zone.edge_from_theta_col[le]]
    θ_t = zone.x_prev[zone.edge_to_theta_col[le]]
    vm_f = zone.x_prev[zone.edge_from_vm_col[le]]
    vm_t = zone.x_prev[zone.edge_to_vm_col[le]]
    return _ac_flow_at(
        zone.edge_g[le], zone.edge_b[le], zone.edge_b_sh[le],
        zone.edge_inv_tap[le], zone.edge_tap_angle[le],
        θ_f, θ_t, vm_f, vm_t,
    )
end

# ------------------------------------------------------------
# Previous-iterate flow on a dummy half (parameters live on the pair)
# ------------------------------------------------------------

@inline function _ac_flow_prev_dummy(
    zone::AladinZoneDataAC, pair::DummyLinePairAC, side::Symbol,
)
    if side === :i
        θ_f = zone.x_prev[pair.from_theta_col_i]
        θ_t = zone.x_prev[pair.to_theta_col_i]
        vm_f = zone.x_prev[pair.from_vm_col_i]
        vm_t = zone.x_prev[pair.to_vm_col_i]
        return _ac_flow_at_to(
            pair.g_i, pair.b_i, pair.b_sh_i,
            pair.inv_tap_i, pair.tap_angle_i,
            θ_f, θ_t, vm_f, vm_t,
        )
    else
        θ_f = zone.x_prev[pair.from_theta_col_j]
        θ_t = zone.x_prev[pair.to_theta_col_j]
        vm_f = zone.x_prev[pair.from_vm_col_j]
        vm_t = zone.x_prev[pair.to_vm_col_j]
        return _ac_flow_at(
            pair.g_j, pair.b_j, pair.b_sh_j,
            pair.inv_tap_j, pair.tap_angle_j,
            θ_f, θ_t, vm_f, vm_t,
        )
    end
end

# ------------------------------------------------------------
# Add (ρ/2)·σ·(p_flow − p_prev)² + (ρ/2)·σ·(q_flow − q_prev)² to obj_expr
# ------------------------------------------------------------

function _ac_prox_flow!(
    obj_expr,
    model::JuMP.Model,
    ws::Union{ACPowerWorkspaceALADIN, ACHybridWorkspaceALADIN},
    zone_pos::Int,
    rho::Float64,
    sigma_pflow::Float64,
    sigma_qflow::Float64,
)
    part_pf = model[:p_flow_branches]
    part_qf = model[:q_flow_branches]
    zone = ws.zones[zone_pos]
    branch_map = ws.partitions.parts[zone_pos].branch_idx_full_to_zone

    half_rho = 0.5 * rho

    for pair in ws.interzonal_line_pairs
        if pair.zone_i == zone_pos
            le = pair.line_local_i
        elseif pair.zone_j == zone_pos
            le = pair.line_local_j
        else
            continue
        end
        line_global = pair.line_global
        l = branch_map[line_global]
        l == 0 && continue

        p_prev, q_prev = _ac_flow_prev_interzonal(zone, le)
        p_expr = part_pf[l, 1]
        q_expr = part_qf[l, 1]
        obj_expr += half_rho * sigma_pflow * (p_expr - p_prev)^2
        obj_expr += half_rho * sigma_qflow * (q_expr - q_prev)^2
    end

    for pair in ws.dummy_line_pairs
        if pair.zone_i == zone_pos
            local_idx = pair.local_i
            side = :i
        elseif pair.zone_j == zone_pos
            local_idx = pair.local_j
            side = :j
        else
            continue
        end

        p_prev, q_prev = _ac_flow_prev_dummy(zone, pair, side)
        if side === :i
            p_expr = -model[:p_reverse_flow_branches][local_idx, 1]
            q_expr = -model[:q_reverse_flow_branches][local_idx, 1]
        else
            p_expr = part_pf[local_idx, 1]
            q_expr = part_qf[local_idx, 1]
        end
        obj_expr += half_rho * sigma_pflow * (p_expr - p_prev)^2
        obj_expr += half_rho * sigma_qflow * (q_expr - q_prev)^2
    end

    return obj_expr
end

# ============================================================
# Update objective
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function update_objective!(part::PartitionAC, ws::ACStateWorkspaceALADIN)
    zone_pos = part.position
    zone_data = ws.zones[zone_pos]

    if haskey(ENV, "LIRION_GAUGEFIX") && !has_reference_bus(part.bus)
        pairs_here = [p for p in ws.bus_pairs
                      if p.zone_i == zone_pos || p.zone_j == zone_pos]
        if !isempty(pairs_here)
            idx_raw = tryparse(Int, ENV["LIRION_GAUGEFIX"])
            idx = clamp(idx_raw === nothing ? 1 : idx_raw, 1, length(pairs_here))
            pair = pairs_here[idx]
            if pair.zone_j == zone_pos
                src = ws.zones[pair.zone_i]
                val = src.x_prev[pair.theta_col_i]
                c = pair.theta_col_j
            else
                src = ws.zones[pair.zone_j]
                val = src.x_prev[pair.theta_col_j]
                c = pair.theta_col_i
            end
            JuMP.fix(zone_data.all_vars[c], val; force = true)
            zone_data.lb[c] = val
            zone_data.ub[c] = val
        end
    end

    lin = zero(JuMP.AffExpr)
    _lambda_bus_ac!(lin, ws, zone_pos)
    _ac_prox_lin!(
        lin, zone_data, ws.rho[];
        sigma_theta = ws.sigma_theta,
        sigma_vm = ws.sigma_vm,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )

    @objective(part.local_model, Min, part.obj_base + zone_data.quad_static + lin)
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function update_objective!(part::PartitionAC, ws::ACPowerWorkspaceALADIN)
    zone_pos = part.position
    zone_data = ws.zones[zone_pos]
    model = part.local_model

    pflow_mode = ws.proximal_mode === :pflow

    sigma_theta_eff = pflow_mode ? ws.sigma_other : ws.sigma_theta
    sigma_vm_eff = pflow_mode ? ws.sigma_other : ws.sigma_vm

    lin = zero(JuMP.AffExpr)
    _ac_prox_lin!(
        lin, zone_data, ws.rho[];
        sigma_theta = sigma_theta_eff,
        sigma_vm = sigma_vm_eff,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )

    obj_expr = part.obj_base + zone_data.quad_static + lin
    obj_expr = _lambda_edge_pq_ac!(model, obj_expr, ws, zone_pos)
    obj_expr = _lambda_dummy_pq_ac!(model, obj_expr, ws, zone_pos)

    if pflow_mode
        obj_expr = _ac_prox_flow!(
            obj_expr, model, ws, zone_pos, ws.rho[],
            ws.sigma_pflow, ws.sigma_qflow,
        )
    end

    @objective(model, Min, obj_expr)
    return nothing
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function update_objective!(part::PartitionAC, ws::ACHybridWorkspaceALADIN)
    zone_pos = part.position
    zone_data = ws.zones[zone_pos]
    model = part.local_model

    pflow_mode = ws.proximal_mode === :pflow

    sigma_theta_eff = pflow_mode ? ws.sigma_other : ws.sigma_theta
    sigma_vm_eff = pflow_mode ? ws.sigma_other : ws.sigma_vm

    lin = zero(JuMP.AffExpr)
    _lambda_bus_ac!(lin, ws, zone_pos)
    _lambda_split_bus_ac!(lin, ws, zone_pos)
    _ac_prox_lin!(
        lin, zone_data, ws.rho[];
        sigma_theta = sigma_theta_eff,
        sigma_vm = sigma_vm_eff,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )

    obj_expr = part.obj_base + zone_data.quad_static + lin
    obj_expr = _lambda_edge_pq_ac!(model, obj_expr, ws, zone_pos)
    obj_expr = _lambda_dummy_pq_ac!(model, obj_expr, ws, zone_pos)

    if pflow_mode
        obj_expr = _ac_prox_flow!(
            obj_expr, model, ws, zone_pos, ws.rho[],
            ws.sigma_pflow, ws.sigma_qflow,
        )
    end

    @objective(model, Min, obj_expr)
    return nothing
end

# ============================================================
# Stationarity-residual contribution of the AC :pflow flow proximal
# ============================================================
function _add_extra_prox_grad!(
    grad::Vector{Float64},
    ws::Union{ACPowerWorkspaceALADIN, ACHybridWorkspaceALADIN},
    zone_pos::Int,
    zone::AladinZoneDataAC,
)
    ws.proximal_mode === :pflow || return nothing
    rho = ws.rho[]
    sp = ws.sigma_pflow
    sq = ws.sigma_qflow

    @inbounds for pair in ws.interzonal_line_pairs
        le = pair.zone_i == zone_pos ? pair.line_local_i :
             pair.zone_j == zone_pos ? pair.line_local_j : 0
        le == 0 && continue
        lin = _ac_flow_linearisation(zone, le)
        p_prev, q_prev = _ac_flow_prev_interzonal(zone, le)
        cp = rho * sp * (lin.pf - p_prev)
        cq = rho * sq * (lin.qf - q_prev)
        grad[zone.edge_from_vm_col[le]] += cp * lin.dpf_dvf + cq * lin.dqf_dvf
        grad[zone.edge_to_vm_col[le]] += cp * lin.dpf_dvt + cq * lin.dqf_dvt
        grad[zone.edge_from_theta_col[le]] +=
            cp * lin.dpf_dtheta_f + cq * lin.dqf_dtheta_f
        grad[zone.edge_to_theta_col[le]] +=
            cp * lin.dpf_dtheta_t + cq * lin.dqf_dtheta_t
    end
    # dummy-bus pairs add the same flow-proximal gradient terms.
    @inbounds for pair in ws.dummy_line_pairs
        if pair.zone_i == zone_pos
            lin = _dummy_ac_flow_linearisation_to(
                zone.y,
                pair.from_theta_col_i, pair.to_theta_col_i,
                pair.from_vm_col_i, pair.to_vm_col_i,
                pair.g_i, pair.b_i, pair.b_sh_i,
                pair.inv_tap_i, pair.tap_angle_i,
            )
            p_prev, q_prev = _ac_flow_prev_dummy(zone, pair, :i)
            fvc, tvc = pair.from_vm_col_i, pair.to_vm_col_i
            ftc, ttc = pair.from_theta_col_i, pair.to_theta_col_i
        elseif pair.zone_j == zone_pos
            lin = _dummy_ac_flow_linearisation(
                zone.y,
                pair.from_theta_col_j, pair.to_theta_col_j,
                pair.from_vm_col_j, pair.to_vm_col_j,
                pair.g_j, pair.b_j, pair.b_sh_j,
                pair.inv_tap_j, pair.tap_angle_j,
            )
            p_prev, q_prev = _ac_flow_prev_dummy(zone, pair, :j)
            fvc, tvc = pair.from_vm_col_j, pair.to_vm_col_j
            ftc, ttc = pair.from_theta_col_j, pair.to_theta_col_j
        else
            continue
        end
        cp = rho * sp * (lin.pf - p_prev)
        cq = rho * sq * (lin.qf - q_prev)
        grad[fvc] += cp * lin.dpf_dvf + cq * lin.dqf_dvf
        grad[tvc] += cp * lin.dpf_dvt + cq * lin.dqf_dvt
        grad[ftc] += cp * lin.dpf_dtheta_f + cq * lin.dqf_dtheta_f
        grad[ttc] += cp * lin.dpf_dtheta_t + cq * lin.dqf_dtheta_t
    end
    return nothing
end

function _extra_prox_support!(
    sup::Set{Int},
    ws::Union{ACPowerWorkspaceALADIN, ACHybridWorkspaceALADIN},
    zone_pos::Int,
    zone::AladinZoneDataAC,
)
    ws.proximal_mode === :pflow || return nothing
    @inbounds for pair in ws.interzonal_line_pairs
        le = pair.zone_i == zone_pos ? pair.line_local_i :
             pair.zone_j == zone_pos ? pair.line_local_j : 0
        le == 0 && continue
        push!(sup, zone.edge_from_vm_col[le])
        push!(sup, zone.edge_to_vm_col[le])
        push!(sup, zone.edge_from_theta_col[le])
        push!(sup, zone.edge_to_theta_col[le])
    end
    @inbounds for pair in ws.dummy_line_pairs
        if pair.zone_i == zone_pos
            push!(sup, pair.from_vm_col_i); push!(sup, pair.to_vm_col_i)
            push!(sup, pair.from_theta_col_i); push!(sup, pair.to_theta_col_i)
        elseif pair.zone_j == zone_pos
            push!(sup, pair.from_vm_col_j); push!(sup, pair.to_vm_col_j)
            push!(sup, pair.from_theta_col_j); push!(sup, pair.to_theta_col_j)
        end
    end
    return nothing
end
