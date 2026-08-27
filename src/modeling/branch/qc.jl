# ============================================================
# QC — Quadratic-Convex relaxation of AC OPF (Hijazi-Coffrin-Van Hentenryck
# 2017, "Convex quadratic relaxations for mixed-integer nonlinear programs
# in power systems", Math. Prog. Comp.). Recursive-McCormick variant (QCRM).
#
# Bus variables (per time step):
#   v_magnitude[b]     — polar voltage magnitude (from AC bus vars)
#   theta[b]           — voltage angle (from AC bus vars)
#   v_magnitude_sq[b]  — w[b] = vm[b]² (SOCWR bus var)
#
# Buspair variables:
#   R[l]    — wr[bp] = vm_i · vm_j · cos(θ_i − θ_j)  (shared per (from,to) pair)
#   T_var[l]— wi[bp] = vm_i · vm_j · sin(θ_i − θ_j)  (shared per (from,to) pair)
#   td[l]   — θ_i − θ_j (angle difference)
#   vv[l]   — vm_i · vm_j
#   cs[l]   — cosine envelope, bounded by cos(angmin..angmax)
#   si[l]   — sine envelope, bounded by sin(angmin..angmax)
#   ccm[l]  — branch current magnitude squared (≥0, bounded by (rate_a·tap/vmin)²)
#
# Relaxations (all linear/quadratic):
#   w[i]  ≥ vm[i]²                 (relaxation_sqr lower)
#   w[i]  ≤ (vmin+vmax)·vm - vmin·vmax  (relaxation_sqr upper — secant)
#   θ_i − θ_j = td[l]             (exact coupling)
#   relaxation_sin(td, si)         (piecewise-linear envelope)
#   relaxation_cos(td, cs)         (quadratic upper + linear lower)
#   relaxation_product(vm_i, vm_j, vv)   (McCormick)
#   relaxation_product(vv, cs, wr)       (McCormick)
#   relaxation_product(vv, si, wi)       (McCormick)
#
# Branch constraints:
#   Ohm's law in W-space (reuse add_socwr_branch_expressions!)
#   p² + q² ≤ (w[fr]/tm²) · ccm    (constraint_power_magnitude_sqr)
#   ccm  == (g²+b²)·(w_fr/tm² + w_to − 2(tr·wr+ti·wi)/tm²)
#           − (g_fr²+b_fr²)·(w_fr/tm²) + 2(g_fr·p_fr − b_fr·q_fr)
#   Thermal SOC on (rate_a, p, q)
#   tan(angmin)·wr ≤ wi ≤ tan(angmax)·wr
#   Lifted Nonlinear Cuts (LNC) — Coffrin-Hijazi-Van Hentenryck 2015
# ============================================================

function add_qc_bus_variables!(model::JuMP.Model, bus::BusData)
    add_voltage_angle_variables!(model, bus)
    add_voltage_magnitude_variables!(model, bus)
    add_socwr_bus_variables!(model, bus)

    # Tight secant relaxation: w ≥ vm² and w ≤ (vmin+vmax)·vm − vmin·vmax
    T = get(model.ext, :time_periods, 1)
    id_to_idx = bus.id_to_idx
    bus_ids = bus.id
    v_min = bus.v_min
    v_max = bus.v_max

    for b in bus_ids
        idx = id_to_idx[b]
        lo = v_min[idx]
        hi = v_max[idx]
        for t in 1:T
            vm = model[:v_magnitude][b, t]
            w = model[:v_magnitude_sq][b, t]
            @constraint(model, w >= vm^2)
            @constraint(model, w <= (lo + hi) * vm - lo * hi)
        end
    end

    return nothing
end

function add_qc_branch_variables!(model::JuMP.Model, bus::BusData, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active
    id_to_idx = bus.id_to_idx
    v_min = bus.v_min
    v_max = bus.v_max

    # wr = R, wi = T_var with PM-equivalent bounds (_rt_product_bounds is the
    # branch-wise version of ref_calc_voltage_product_bounds).
    add_rt_variables_with_bounds!(model, bus, br; bounded = true)

    @variable(model, td[l in active, t in 1:T])
    @variable(model, vv[l in active, t in 1:T])
    @variable(model, cs[l in active, t in 1:T])
    @variable(model, si[l in active, t in 1:T])
    @variable(model, ccm[l in active, t in 1:T], lower_bound = 0.0)

    for l in active
        amin = br.ang_min[l]
        amax = br.ang_max[l]
        idx_i = id_to_idx[br.from[l]]
        idx_j = id_to_idx[br.to[l]]
        vi_min, vi_max = v_min[idx_i], v_max[idx_i]
        vj_min, vj_max = v_min[idx_j], v_max[idx_j]

        # td bounds
        set_lower_bound.(td[l, :], amin)
        set_upper_bound.(td[l, :], amax)

        # vv bounds
        set_lower_bound.(vv[l, :], vi_min * vj_min)
        set_upper_bound.(vv[l, :], vi_max * vj_max)

        # cs bounds
        if amin >= 0
            cs_up = cos(amin)
            cs_lo = cos(amax)
        elseif amax <= 0
            cs_up = cos(amax)
            cs_lo = cos(amin)
        else
            cs_up = 1.0
            cs_lo = min(cos(amin), cos(amax))
        end
        set_lower_bound.(cs[l, :], cs_lo)
        set_upper_bound.(cs[l, :], cs_up)
        set_start_value.(cs[l, :], 1.0)

        # si bounds
        set_lower_bound.(si[l, :], sin(amin))
        set_upper_bound.(si[l, :], sin(amax))

        # ccm upper bound using PM formula: ((rate_a · tap)/vmin_fr)²
        # PM uses vm_fr_min = bus[f_bus]["vmin"], here vi_min.
        tap = br.tap_position[l]
        ra = br.rate_a[l]
        if isfinite(ra) && ra > 0
            ub = ((ra * tap) / vi_min)^2
            set_upper_bound.(ccm[l, :], ub)
        end

        set_start_value.(vv[l, :], 1.0)
    end

    return nothing
end

function add_qc_envelope_constraints!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    td = model[:td]
    si = model[:si]
    cs = model[:cs]
    vv = model[:vv]
    R = model[:R]
    T_var = model[:T_var]
    vm = model[:v_magnitude]

    # td coupling: theta_i - theta_j == td
    theta = model[:theta]
    @constraint(model, [l in active, t in 1:T],
        theta[br.from[l], t] - theta[br.to[l], t] == td[l, t]
    )

    for l in active, t in 1:T
        # relaxation_sin(td, si)
        _qc_relaxation_sin!(model, td[l, t], si[l, t])
        # relaxation_cos(td, cs)
        _qc_relaxation_cos!(model, td[l, t], cs[l, t])
        # vv = vm_i · vm_j (McCormick)
        _qc_relaxation_product!(model, vm[br.from[l], t], vm[br.to[l], t], vv[l, t])
        # wr = vv · cs (McCormick)
        _qc_relaxation_product!(model, vv[l, t], cs[l, t], R[l, t])
        # wi = vv · si (McCormick)
        _qc_relaxation_product!(model, vv[l, t], si[l, t], T_var[l, t])
    end

    return nothing
end

# --- envelope helpers (mirror PM/IM forms) ---

@inline function _qc_variable_domain(v::JuMP.VariableRef)
    lb = JuMP.has_lower_bound(v) ? JuMP.lower_bound(v) : -Inf
    ub = JuMP.has_upper_bound(v) ? JuMP.upper_bound(v) : Inf
    return lb, ub
end

# Piecewise-linear envelope for y = sin(x), x ∈ [-π/2, π/2].
function _qc_relaxation_sin!(model::JuMP.Model, x::JuMP.VariableRef, y::JuMP.VariableRef)
    lb, ub = _qc_variable_domain(x)
    @assert lb >= -pi / 2 && ub <= pi / 2
    max_ad = max(abs(lb), abs(ub))

    if lb < 0 && ub > 0
        @constraint(model, y <= cos(max_ad / 2) * (x - max_ad / 2) + sin(max_ad / 2))
        @constraint(model, y >= cos(max_ad / 2) * (x + max_ad / 2) - sin(max_ad / 2))
    elseif ub <= 0
        @constraint(model, y <= (sin(lb) - sin(ub)) / (lb - ub) * (x - lb) + sin(lb))
        @constraint(model, y >= cos(max_ad / 2) * (x + max_ad / 2) - sin(max_ad / 2))
    else  # lb >= 0
        @constraint(model, y <= cos(max_ad / 2) * (x - max_ad / 2) + sin(max_ad / 2))
        @constraint(model, y >= (sin(lb) - sin(ub)) / (lb - ub) * (x - lb) + sin(lb))
    end
    return nothing
end

# Quadratic-upper + linear-lower envelope for y = cos(x), x ∈ [-π/2, π/2].
function _qc_relaxation_cos!(model::JuMP.Model, x::JuMP.VariableRef, y::JuMP.VariableRef)
    lb, ub = _qc_variable_domain(x)
    @assert lb >= -pi / 2 && ub <= pi / 2
    max_ad = max(abs(lb), abs(ub))

    @constraint(model, y <= 1 - (1 - cos(max_ad)) / (max_ad * max_ad) * (x^2))
    @constraint(model, y >= (cos(lb) - cos(ub)) / (lb - ub) * (x - lb) + cos(lb))
    return nothing
end

# McCormick envelope for z = x·y.
function _qc_relaxation_product!(
    model::JuMP.Model, x::JuMP.VariableRef, y::JuMP.VariableRef, z::JuMP.VariableRef,
)
    x_lb, x_ub = _qc_variable_domain(x)
    y_lb, y_ub = _qc_variable_domain(y)

    @constraint(model, z >= x_lb * y + y_lb * x - x_lb * y_lb)
    @constraint(model, z >= x_ub * y + y_ub * x - x_ub * y_ub)
    @constraint(model, z <= x_lb * y + y_ub * x - x_lb * y_ub)
    @constraint(model, z <= x_ub * y + y_lb * x - x_ub * y_lb)
    return nothing
end

# Same LNC block used by SDPWRM — strengthens the W-space relaxation.
function add_qc_lnc_cuts!(model::JuMP.Model, bus::BusData, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active
    id_to_idx = bus.id_to_idx
    v_min = bus.v_min
    v_max = bus.v_max

    for l in active
        amin = br.ang_min[l]
        amax = br.ang_max[l]
        (isnan(amin) || isnan(amax)) && continue
        amin >= -pi / 2 && amin <= pi / 2 || continue
        amax >= -pi / 2 && amax <= pi / 2 || continue
        amin < amax || continue

        idx_i = id_to_idx[br.from[l]]
        idx_j = id_to_idx[br.to[l]]
        vf_lb, vf_ub = v_min[idx_i], v_max[idx_i]
        vt_lb, vt_ub = v_min[idx_j], v_max[idx_j]

        phi = (amax + amin) / 2
        d = (amax - amin) / 2
        sf = vf_lb + vf_ub
        st = vt_lb + vt_ub

        for t in 1:T
            wf = model[:v_magnitude_sq][br.from[l], t]
            wt = model[:v_magnitude_sq][br.to[l], t]
            wr = model[:R][l, t]
            wi = model[:T_var][l, t]

            @constraint(model,
                sf * st * (cos(phi) * wr + sin(phi) * wi) -
                vt_ub * cos(d) * st * wf -
                vf_ub * cos(d) * sf * wt >=
                vf_ub * vt_ub * cos(d) * (vf_lb * vt_lb - vf_ub * vt_ub)
            )
            @constraint(model,
                sf * st * (cos(phi) * wr + sin(phi) * wi) -
                vt_lb * cos(d) * st * wf -
                vf_lb * cos(d) * sf * wt >=
                -vf_lb * vt_lb * cos(d) * (vf_lb * vt_lb - vf_ub * vt_ub)
            )
        end
    end
    return nothing
end

function add_qc_branch_constraints!(model::JuMP.Model, bus::BusData, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active
    id_to_idx = bus.id_to_idx

    # Envelopes (sqr + sin + cos + McCormick chain)
    add_qc_envelope_constraints!(model, br)

    # Thermal (quadratic — QC is an NLP, use Ipopt-friendly form)
    @constraint(model, [l in active, t in 1:T],
        model[:p_flow_branches][l, t]^2 + model[:q_flow_branches][l, t]^2 <=
        br.rate_a[l]^2
    )
    @constraint(model, [l in active, t in 1:T],
        model[:p_reverse_flow_branches][l, t]^2 +
        model[:q_reverse_flow_branches][l, t]^2 <= br.rate_a[l]^2
    )

    # Angle-difference tangent cuts on (wr, wi)
    add_rt_tangent_angle_cuts!(model, br)

    # Lifted Nonlinear Cuts
    add_qc_lnc_cuts!(model, bus, br)

    # Power-magnitude link: p² + q² ≤ w_fr/tm² · ccm
    ccm = model[:ccm]
    for l in active, t in 1:T
        w_fr = model[:v_magnitude_sq][br.from[l], t]
        tm = br.tap_position[l]
        p_fr = model[:p_flow_branches][l, t]
        q_fr = model[:q_flow_branches][l, t]
        @constraint(model, p_fr^2 + q_fr^2 <= w_fr / (tm^2) * ccm[l, t])
    end

    # ccm algebraic link:
    # ccm == (g²+b²)·(w_fr/tm² + w_to − 2(tr·wr + ti·wi)/tm²)
    #        − (g_fr²+b_fr²)·(w_fr/tm²) + 2(g_fr·p_fr − b_fr·q_fr)
    for l in active, t in 1:T
        w_fr = model[:v_magnitude_sq][br.from[l], t]
        w_to = model[:v_magnitude_sq][br.to[l], t]
        wr = model[:R][l, t]
        wi = model[:T_var][l, t]
        p_fr = model[:p_flow_branches][l, t]
        q_fr = model[:q_flow_branches][l, t]
        g = br.g[l]
        b = br.b[l]
        g_fr = br.g_fr[l]
        b_fr = br.b_fr[l]
        tr = br.tap_real[l]
        ti = br.tap_imag[l]
        tm = br.tap_position[l]
        ym_sh_sqr = g_fr^2 + b_fr^2

        @constraint(model,
            ccm[l, t] ==
            (g^2 + b^2) * (w_fr / tm^2 + w_to - 2 * (tr * wr + ti * wi) / tm^2) -
            ym_sh_sqr * (w_fr / tm^2) +
            2 * (g_fr * p_fr - b_fr * q_fr)
        )
    end

    return nothing
end
