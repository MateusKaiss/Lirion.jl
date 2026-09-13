# ============================================================
# ALADIN-DC — Phase 2 line search (Armijo backtracking)
# ============================================================

const _DCPhase2Workspace = Union{
    DCStateWorkspaceALADIN,
    DCPowerWorkspaceALADIN,
    DCHybridWorkspaceALADIN,
}

# ============================================================
# Per-zone cost trial evaluations
# ============================================================

@inline function _dc_cost_at_trial(zone::AladinZoneData, alpha::Float64)
    acc = 0.0
    y = zone.y
    delta = zone.delta
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        xk = y[col] + alpha * delta[col]
        acc += 0.5 * zone.pg_c2_sbase2[k] * xk * xk + zone.pg_c1_sbase[k] * xk
    end
    return acc
end

@inline function _dc_cost_dir_deriv_0(zone::AladinZoneData)
    acc = 0.0
    y = zone.y
    delta = zone.delta
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        acc += (zone.pg_c2_sbase2[k] * y[col] + zone.pg_c1_sbase[k]) * delta[col]
    end
    return acc
end

# ============================================================
# Interzonal line flow trial evaluations
# ============================================================

@inline function _dc_pflow_at(zone::AladinZoneData, le::Int, alpha::Float64)
    from_col = zone.edge_from_cols[le]
    to_col = zone.edge_to_cols[le]
    b_val = zone.edge_b[le]
    yf = zone.y[from_col] + alpha * zone.delta[from_col]
    yt = zone.y[to_col] + alpha * zone.delta[to_col]
    return b_val * (yf - yt)
end

@inline function _dc_pflow_dir_deriv(zone::AladinZoneData, le::Int)
    from_col = zone.edge_from_cols[le]
    to_col = zone.edge_to_cols[le]
    b_val = zone.edge_b[le]
    return b_val * (zone.delta[from_col] - zone.delta[to_col])
end

# ============================================================
# Dummy line flow trial evaluations
# ============================================================

@inline function _dc_dummy_pflow_at(
    zone::AladinZoneData,
    from_col::Int, to_col::Int, b_val::Float64,
    alpha::Float64,
)
    yf = zone.y[from_col] + alpha * zone.delta[from_col]
    yt = zone.y[to_col] + alpha * zone.delta[to_col]
    return b_val * (yf - yt)
end

@inline function _dc_dummy_pflow_dir_deriv(
    zone::AladinZoneData,
    from_col::Int, to_col::Int, b_val::Float64,
)
    return b_val * (zone.delta[from_col] - zone.delta[to_col])
end

# ============================================================
# Merit function phi(alpha)
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function _dc_merit(ws::DCStateWorkspaceALADIN, alpha::Float64)
    phi = 0.0
    for zone in ws.zones
        phi += _dc_cost_at_trial(zone, alpha)
    end
    r_sq = 0.0
    lambda_r = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.col_i] + alpha * zi.delta[pair.col_i]
        xj = zj.y[pair.col_j] + alpha * zj.delta[pair.col_j]
        d = xi - xj
        r_sq += d * d
        lambda_r += ws.vars.theta.lambda[k] * d
    end
    return phi + lambda_r + 0.5 * ws.mu[] * r_sq
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function _dc_merit(ws::DCPowerWorkspaceALADIN, alpha::Float64)
    phi = 0.0
    for zone in ws.zones
        phi += _dc_cost_at_trial(zone, alpha)
    end
    r_sq = 0.0
    lambda_r = 0.0
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _dc_pflow_at(zi, pair.line_local_i, alpha)
        pj_val = _dc_pflow_at(zj, pair.line_local_j, alpha)
        d = pi_val - pj_val
        r_sq += d * d
        lambda_r += ws.vars.p_flow.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _dc_dummy_pflow_at(
            zi, pair.from_col_i, pair.to_col_i, pair.b_i, alpha,
        )
        pj_val = _dc_dummy_pflow_at(
            zj, pair.from_col_j, pair.to_col_j, pair.b_j, alpha,
        )
        d = pi_val - pj_val
        r_sq += d * d
        lambda_r += ws.vars.p_dummy.lambda[k] * d
    end
    return phi + lambda_r + 0.5 * ws.mu[] * r_sq
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function _dc_merit(ws::DCHybridWorkspaceALADIN, alpha::Float64)
    phi = 0.0
    for zone in ws.zones
        phi += _dc_cost_at_trial(zone, alpha)
    end
    r_sq = 0.0
    lambda_r = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.col_i] + alpha * zi.delta[pair.col_i]
        xj = zj.y[pair.col_j] + alpha * zj.delta[pair.col_j]
        d = xi - xj
        r_sq += d * d
        lambda_r += ws.vars.theta.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _dc_pflow_at(zi, pair.line_local_i, alpha)
        pj_val = _dc_pflow_at(zj, pair.line_local_j, alpha)
        d = pi_val - pj_val
        r_sq += d * d
        lambda_r += ws.vars.p_flow.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _dc_dummy_pflow_at(
            zi, pair.from_col_i, pair.to_col_i, pair.b_i, alpha,
        )
        pj_val = _dc_dummy_pflow_at(
            zj, pair.from_col_j, pair.to_col_j, pair.b_j, alpha,
        )
        d = pi_val - pj_val
        r_sq += d * d
        lambda_r += ws.vars.p_dummy.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        d = 0.0
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            col = pair.p_inj_cols[c]
            d += zone.y[col] + alpha * zone.delta[col]
        end
        r_sq += d * d
        lambda_r += ws.vars.p_inj.lambda[k] * d
    end
    return phi + lambda_r + 0.5 * ws.mu[] * r_sq
end

# ============================================================
# Directional derivative dphi/dalpha at alpha = 0
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function _dc_dphi_0(ws::DCStateWorkspaceALADIN)
    dphi = 0.0
    for zone in ws.zones
        dphi += _dc_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 = zi.y[pair.col_i] - zj.y[pair.col_j]
        dr = zi.delta[pair.col_i] - zj.delta[pair.col_j]
        dphi += (ws.vars.theta.lambda[k] + mu * r0) * dr
    end
    return dphi
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function _dc_dphi_0(ws::DCPowerWorkspaceALADIN)
    dphi = 0.0
    for zone in ws.zones
        dphi += _dc_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 =
            _dc_pflow_at(zi, pair.line_local_i, 0.0) -
            _dc_pflow_at(zj, pair.line_local_j, 0.0)
        dr =
            _dc_pflow_dir_deriv(zi, pair.line_local_i) -
            _dc_pflow_dir_deriv(zj, pair.line_local_j)
        dphi += (ws.vars.p_flow.lambda[k] + mu * r0) * dr
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 =
            _dc_dummy_pflow_at(zi, pair.from_col_i, pair.to_col_i, pair.b_i, 0.0) -
            _dc_dummy_pflow_at(zj, pair.from_col_j, pair.to_col_j, pair.b_j, 0.0)
        dr =
            _dc_dummy_pflow_dir_deriv(zi, pair.from_col_i, pair.to_col_i, pair.b_i) -
            _dc_dummy_pflow_dir_deriv(zj, pair.from_col_j, pair.to_col_j, pair.b_j)
        dphi += (ws.vars.p_dummy.lambda[k] + mu * r0) * dr
    end
    return dphi
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function _dc_dphi_0(ws::DCHybridWorkspaceALADIN)
    dphi = 0.0
    for zone in ws.zones
        dphi += _dc_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 = zi.y[pair.col_i] - zj.y[pair.col_j]
        dr = zi.delta[pair.col_i] - zj.delta[pair.col_j]
        dphi += (ws.vars.theta.lambda[k] + mu * r0) * dr
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 =
            _dc_pflow_at(zi, pair.line_local_i, 0.0) -
            _dc_pflow_at(zj, pair.line_local_j, 0.0)
        dr =
            _dc_pflow_dir_deriv(zi, pair.line_local_i) -
            _dc_pflow_dir_deriv(zj, pair.line_local_j)
        dphi += (ws.vars.p_flow.lambda[k] + mu * r0) * dr
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 =
            _dc_dummy_pflow_at(zi, pair.from_col_i, pair.to_col_i, pair.b_i, 0.0) -
            _dc_dummy_pflow_at(zj, pair.from_col_j, pair.to_col_j, pair.b_j, 0.0)
        dr =
            _dc_dummy_pflow_dir_deriv(zi, pair.from_col_i, pair.to_col_i, pair.b_i) -
            _dc_dummy_pflow_dir_deriv(zj, pair.from_col_j, pair.to_col_j, pair.b_j)
        dphi += (ws.vars.p_dummy.lambda[k] + mu * r0) * dr
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        r0 = 0.0
        dr = 0.0
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            col = pair.p_inj_cols[c]
            r0 += zone.y[col]
            dr += zone.delta[col]
        end
        dphi += (ws.vars.p_inj.lambda[k] + mu * r0) * dr
    end
    return dphi
end

# ============================================================
# HFD Algorithm 3 globalization — L1 exact penalty Φ, steps (a)/(b)/(c)
# ============================================================

@inline function _dc_cost(zone::AladinZoneData, acc::F) where {F}
    total = 0.0
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        xk = acc(zone, col)
        total += 0.5 * zone.pg_c2_sbase2[k] * xk * xk + zone.pg_c1_sbase[k] * xk
    end
    return total
end

@inline function _dc_pflow(zone::AladinZoneData, le::Int, acc::F) where {F}
    return zone.edge_b[le] *
           (acc(zone, zone.edge_from_cols[le]) - acc(zone, zone.edge_to_cols[le]))
end

@inline function _dc_dummy_pf(
    zone::AladinZoneData, from_col::Int, to_col::Int, b_val::Float64, acc::F,
) where {F}
    return b_val * (acc(zone, from_col) - acc(zone, to_col))
end

_dc_cost_sum(ws::_DCPhase2Workspace, acc::F) where {F} =
    sum(_dc_cost(zone, acc) for zone in ws.zones)

# ℓ1 consensus residual ‖Ax − b‖₁ at the point defined by `acc` (b = 0).
function _dc_resid_l1(ws::DCStateWorkspaceALADIN, acc::F) where {F}
    r = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(acc(zi, pair.col_i) - acc(zj, pair.col_j))
    end
    return r
end

function _dc_resid_l1(ws::DCPowerWorkspaceALADIN, acc::F) where {F}
    r = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(_dc_pflow(zi, pair.line_local_i, acc) -
                 _dc_pflow(zj, pair.line_local_j, acc))
    end
    @inbounds for pair in ws.dummy_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(_dc_dummy_pf(zi, pair.from_col_i, pair.to_col_i, pair.b_i, acc) -
                 _dc_dummy_pf(zj, pair.from_col_j, pair.to_col_j, pair.b_j, acc))
    end
    return r
end

function _dc_resid_l1(ws::DCHybridWorkspaceALADIN, acc::F) where {F}
    r = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(acc(zi, pair.col_i) - acc(zj, pair.col_j))
    end
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(_dc_pflow(zi, pair.line_local_i, acc) -
                 _dc_pflow(zj, pair.line_local_j, acc))
    end
    @inbounds for pair in ws.dummy_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(_dc_dummy_pf(zi, pair.from_col_i, pair.to_col_i, pair.b_i, acc) -
                 _dc_dummy_pf(zj, pair.from_col_j, pair.to_col_j, pair.b_j, acc))
    end
    @inbounds for pair in ws.split_bus_pairs
        d = 0.0
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            d += acc(zone, pair.p_inj_cols[c])
        end
        r += abs(d)
    end
    return r
end

function _dc_dirderiv(ws::DCStateWorkspaceALADIN, dir, acc::F) where {F}
    s = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += dir.theta[k] * (acc(zi, pair.col_i) - acc(zj, pair.col_j))
    end
    return s
end

function _dc_dirderiv(ws::DCPowerWorkspaceALADIN, dir, acc::F) where {F}
    s = 0.0
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += dir.p_flow[k] * (_dc_pflow(zi, pair.line_local_i, acc) -
                              _dc_pflow(zj, pair.line_local_j, acc))
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += dir.p_dummy[k] *
             (_dc_dummy_pf(zi, pair.from_col_i, pair.to_col_i, pair.b_i, acc) -
              _dc_dummy_pf(zj, pair.from_col_j, pair.to_col_j, pair.b_j, acc))
    end
    return s
end

function _dc_dirderiv(ws::DCHybridWorkspaceALADIN, dir, acc::F) where {F}
    s = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += dir.theta[k] * (acc(zi, pair.col_i) - acc(zj, pair.col_j))
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += dir.p_flow[k] * (_dc_pflow(zi, pair.line_local_i, acc) -
                              _dc_pflow(zj, pair.line_local_j, acc))
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += dir.p_dummy[k] *
             (_dc_dummy_pf(zi, pair.from_col_i, pair.to_col_i, pair.b_i, acc) -
              _dc_dummy_pf(zj, pair.from_col_j, pair.to_col_j, pair.b_j, acc))
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        d = 0.0
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            d += acc(zone, pair.p_inj_cols[c])
        end
        s += dir.p_inj[k] * d
    end
    return s
end

function _dual_slope(ws::_DCPhase2Workspace, dirs)
    dir = dirs isa NamedTuple ? dirs : NamedTuple{keys(ws.vars)}(Tuple(dirs))
    return _dc_dirderiv(ws, dir, (z, c) -> JuMP.value(z.all_vars[c]))
end

function _dc_residual_rows(ws::DCStateWorkspaceALADIN, acc::F) where {F}
    r_th = Vector{Float64}(undef, length(ws.bus_pairs))
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_th[k] = acc(zi, pair.col_i) - acc(zj, pair.col_j)
    end
    return NamedTuple{keys(ws.vars)}((theta = r_th,))
end

function _dc_residual_rows(ws::DCPowerWorkspaceALADIN, acc::F) where {F}
    r_p = Vector{Float64}(undef, length(ws.interzonal_line_pairs))
    r_pd = Vector{Float64}(undef, length(ws.dummy_line_pairs))
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_p[k] = _dc_pflow(zi, pair.line_local_i, acc) -
                 _dc_pflow(zj, pair.line_local_j, acc)
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_pd[k] = _dc_dummy_pf(zi, pair.from_col_i, pair.to_col_i, pair.b_i, acc) -
                  _dc_dummy_pf(zj, pair.from_col_j, pair.to_col_j, pair.b_j, acc)
    end
    return NamedTuple{keys(ws.vars)}((p_flow = r_p, p_dummy = r_pd))
end

function _dc_residual_rows(ws::DCHybridWorkspaceALADIN, acc::F) where {F}
    r_th = Vector{Float64}(undef, length(ws.bus_pairs))
    r_p = Vector{Float64}(undef, length(ws.interzonal_line_pairs))
    r_pd = Vector{Float64}(undef, length(ws.dummy_line_pairs))
    r_pi = Vector{Float64}(undef, length(ws.split_bus_pairs))
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_th[k] = acc(zi, pair.col_i) - acc(zj, pair.col_j)
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_p[k] = _dc_pflow(zi, pair.line_local_i, acc) -
                 _dc_pflow(zj, pair.line_local_j, acc)
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_pd[k] = _dc_dummy_pf(zi, pair.from_col_i, pair.to_col_i, pair.b_i, acc) -
                  _dc_dummy_pf(zj, pair.from_col_j, pair.to_col_j, pair.b_j, acc)
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        d = 0.0
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            d += acc(zone, pair.p_inj_cols[c])
        end
        r_pi[k] = d
    end
    return NamedTuple{keys(ws.vars)}(
        (theta = r_th, p_flow = r_p, p_inj = r_pi, p_dummy = r_pd),
    )
end

_dual_residual(ws::_DCPhase2Workspace) =
    _dc_residual_rows(ws, (z, c) -> JuMP.value(z.all_vars[c]))

# Σ ½ dᵀ(ρΣ)d = value of the proximal quadratic at d = y − x (mirrors _ac_prox).
function _dc_prox(ws::_DCPhase2Workspace)
    s = 0.0
    @inbounds for zone in ws.zones
        v2c = zone.var_to_col
        y = zone.y; x = zone.x_prev
        for (pair, a) in zone.quad_static.terms
            i = v2c[pair.a]; j = v2c[pair.b]
            s += a * (y[i] - x[i]) * (y[j] - x[j])
        end
    end
    return s
end

function _dc_local_infeas_l1(ws::_DCPhase2Workspace, acc::F) where {F}
    total = 0.0
    for zone in ws.zones
        model = zone.local_model
        v2c = zone.var_to_col
        valfn = v -> acc(zone, v2c[v])
        for (Fc, S) in JuMP.list_of_constraint_types(model)
            Fc <: JuMP.VariableRef && continue
            for ci in JuMP.all_constraints(model, Fc, S)
                cobj = JuMP.constraint_object(ci)
                g = try
                    JuMP.value(valfn, cobj.func)
                catch
                    return Inf
                end
                s = cobj.set
                if s isa MOI.EqualTo{Float64}
                    total += abs(g - s.value)
                elseif s isa MOI.LessThan{Float64}
                    total += max(0.0, g - s.upper)
                elseif s isa MOI.GreaterThan{Float64}
                    total += max(0.0, s.lower - g)
                elseif s isa MOI.Interval{Float64}
                    total += max(0.0, g - s.upper) + max(0.0, s.lower - g)
                end
            end
        end
        lb = zone.lb; ub = zone.ub
        @inbounds for c in 1:zone.n_vars
            v = acc(zone, c)
            isnan(lb[c]) || (total += max(0.0, lb[c] - v))
            isnan(ub[c]) || (total += max(0.0, v - ub[c]))
        end
    end
    return total
end

function _dc_max_local_dual(ws::_DCPhase2Workspace; include_bounds::Bool = false)
    maxkap = 0.0
    for zone in ws.zones
        model = zone.local_model
        (try; JuMP.has_duals(model); catch; false; end) || continue
        for (Fc, S) in JuMP.list_of_constraint_types(model)
            Fc <: JuMP.VariableRef && continue
            for ci in JuMP.all_constraints(model, Fc, S)
                d = try
                    JuMP.dual(ci)
                catch
                    continue
                end
                d isa Real || continue
                ad = abs(d)
                ad > maxkap && (maxkap = ad)
            end
        end
        include_bounds || continue
        @inbounds for v in zone.all_vars
            if JuMP.has_lower_bound(v)
                d = try; abs(JuMP.dual(JuMP.LowerBoundRef(v))); catch; 0.0; end
                d > maxkap && (maxkap = d)
            end
            if JuMP.has_upper_bound(v)
                d = try; abs(JuMP.dual(JuMP.UpperBoundRef(v))); catch; 0.0; end
                d > maxkap && (maxkap = d)
            end
        end
    end
    return maxkap
end

const _ALG3_ORIGINAL_LAM_BAR_FLOOR = 1.0

const _ALG3_ORIGINAL_META_KEYS = (
    :a, :b, :c, :c_full, :c_wolfe, :c_budget, :c_no_ascent,
    :c_full_grad, :c_wolfe_grad, :c_budget_grad, :lam_bar_violations,
    :kap_bar_violations,
)

function _alg3_original_meta(ws::_DCPhase2Workspace)
    d = get!(ws.hist.meta, :alg3_original) do
        Dict{Symbol, Int}(k => 0 for k in _ALG3_ORIGINAL_META_KEYS)
    end
    return d::Dict{Symbol, Int}
end

@inline function _alg3_original_count!(ws::_DCPhase2Workspace, key::Symbol)
    d = _alg3_original_meta(ws)
    d[key] = get(d, key, 0) + 1
    return nothing
end

@inline function _alg3_original_c_key(status::Symbol)
    status === :full_step && return :c_full
    status === :wolfe && return :c_wolfe
    status === :no_ascent && return :c_no_ascent
    status === :full_step_grad && return :c_full_grad
    status === :wolfe_grad && return :c_wolfe_grad
    (status === :armijo_grad || status === :ascent_grad) && return :c_budget_grad
    return :c_budget          # :armijo_only / :ascent_only
end

function _dc_alg3_step!(
    ws::_DCPhase2Workspace;
    gamma::Float64 = 1e-4,
    original::Bool = false,
)
    accx = (z, c) -> z.x_prev[c]
    accy = (z, c) -> z.y[c]
    accyd = (z, c) -> z.y[c] + z.delta[c]
    maxlam = 0.0
    for state in values(ws.vars)
        @inbounds for v in state.lambda
            av = abs(v)
            av > maxlam && (maxlam = av)
        end
    end
    use_viol = _ALG3_PHI_LOCAL_VIOL[]
    if original
        if ws.lam_bar <= 0.0
            scale = maxlam
            for state in values(ws.vars)
                @inbounds for v in state.nu
                    av = abs(v)
                    av > scale && (scale = av)
                end
            end
            use_viol && (scale = max(scale, _dc_max_local_dual(ws; include_bounds = true)))
            ws.lam_bar = max(10.0 * scale, _ALG3_ORIGINAL_LAM_BAR_FLOOR)
        elseif maxlam > ws.lam_bar
            _alg3_original_count!(ws, :lam_bar_violations)
            @warn "ALADIN :alg3_original — |λ|∞ = $(maxlam) exceeds the fixed " *
                  "L1 penalty weight λ̄ = $(ws.lam_bar) (iter $(ws.iter)); " *
                  "Φ is no longer guaranteed exact." maxlog = 1
        end
        if use_viol && ws.lam_bar > 0.0
            maxkap = _dc_max_local_dual(ws; include_bounds = true)
            if maxkap > ws.lam_bar
                _alg3_original_count!(ws, :kap_bar_violations)
                @warn "ALADIN :alg3_original — max|κ_local| = $(maxkap) exceeds the " *
                      "fixed local-penalty weight κ̄ = λ̄ = $(ws.lam_bar) (iter $(ws.iter))." maxlog = 1
            end
        end
        if !haskey(ws.hist.meta, :alg3_lam_bar)
            # first globalization call: record the λ̄ in force for the whole
            # run (automatic or user-fixed) and materialise the counters
            ws.hist.meta[:alg3_lam_bar] = ws.lam_bar
            ws.hist.meta[:alg3_lam_bar_fixed] = ws.lam_bar_fixed
            _alg3_original_meta(ws)
        end
    else
        use_viol && (maxlam = max(maxlam, _dc_max_local_dual(ws)))
        if ws.lam_bar <= 0.0 || ws.lam_bar < 1.5 * maxlam
            ws.lam_bar = max(10.0 * maxlam, 1.0)
        end
    end
    lam_bar = ws.lam_bar

    infeas = use_viol ?
        (acc -> lam_bar * _dc_local_infeas_l1(ws, acc)) : (acc -> 0.0)

    rhs = gamma * (_dc_prox(ws) + lam_bar * _dc_resid_l1(ws, accy))
    phi_x = _dc_cost_sum(ws, accx) + lam_bar * _dc_resid_l1(ws, accx) + infeas(accx)

    # branch (a): full step x⁺ = y + Δy
    phi_a = _dc_cost_sum(ws, accyd) + lam_bar * _dc_resid_l1(ws, accyd) + infeas(accyd)
    if isfinite(phi_a) && phi_x - phi_a >= rhs
        original && _alg3_original_count!(ws, :a)
        ws.alpha1 = 1.0; ws.alpha2 = 1.0; ws.alpha3 = 1.0
        return nothing
    end

    # branch (b): x⁺ = y (y is local-feasible by construction).
    phi_b = _dc_cost_sum(ws, accy) + lam_bar * _dc_resid_l1(ws, accy) + infeas(accy)
    if isfinite(phi_b) && phi_x - phi_b >= rhs
        original && _alg3_original_count!(ws, :b)
        ws.alpha1 = 1.0; ws.alpha2 = 0.0; ws.alpha3 = 0.0
        return nothing
    end

    # branch (c): dual-only step.
    #   :alg3          — shared core golden-section on ν−λ (6 probes).
    #   :alg3_original — footnote-5 inexact Wolfe search, α₃ = 1 reachable.
    if original
        _alg3_original_count!(ws, :c)
        status = _branch_c_dual_linesearch_wolfe!(ws)
        _alg3_original_count!(ws, _alg3_original_c_key(status))
        return status
    end
    return _branch_c_dual_linesearch!(ws)
end

# ============================================================
# Armijo backtracking step
# ============================================================

function _line_search_step!(
    ws::_DCPhase2Workspace;
    c1::Float64 = 1e-4,
    max_backtrack::Int = 25,
    alpha_min::Float64 = 1e-3,
)
    ws.line_search || return nothing
    ws.globalization === :alg3 && return _dc_alg3_step!(ws)
    ws.globalization === :alg3_original &&
        return _dc_alg3_step!(ws; original = true)
    # :full_step — HFD Alg.2 step 5 undamped step (α₁=α₂=1, α₃=alpha3_cfg).
    ws.globalization === :full_step && return nothing

    phi_0 = _dc_merit(ws, 0.0)
    dphi_0 = _dc_dphi_0(ws)

    if !isfinite(dphi_0) || dphi_0 >= 0.0
        ws.alpha2 = alpha_min
        return nothing
    end

    alpha = ws.alpha2_init
    @inbounds for _ in 1:max_backtrack
        phi_alpha = _dc_merit(ws, alpha)
        if isfinite(phi_alpha) && phi_alpha <= phi_0 + c1 * alpha * dphi_0
            ws.alpha2 = alpha
            return nothing
        end
        alpha *= 0.5
        if alpha < alpha_min
            break
        end
    end
    ws.alpha2 = alpha_min
    return nothing
end
