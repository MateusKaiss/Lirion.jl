# =======================================================================
# SDPWRM — Semidefinite W-matrix relaxation (Bai-Wei-Fujisawa 2008).
#
# Variables per time step t:
#   WR_sdp[t]  — n×n Symmetric real matrix, WR[i,i] ≥ 0 (v² on diagonal)
#   WI_sdp[t]  — n×n real matrix, antisymmetric by the PSD embedding
#
# PSD constraint (per t):
#   [ WR   WI ]
#   [ -WI  WR ]  ⪰ 0
#
# The branch-flow equations reuse the SOCWR Ohm's law; we expose
# `v_magnitude_sq[b,t]`, `R[l,t]`, `T_var[l,t]` as JuMP expressions that
# alias the matrix entries, so `add_socwr_branch_expressions!` can be
# reused verbatim.
# =======================================================================

function _sdpwrm_buspair_bounds(bus::BusData, br::BranchData)
    n = bus.n
    active = br.active
    id_to_idx = bus.id_to_idx

    WR_min = fill(-Inf, n, n)
    WR_max = fill(Inf, n, n)
    WI_min = fill(-Inf, n, n)
    WI_max = fill(Inf, n, n)

    v_min = bus.v_min
    v_max = bus.v_max

    # diagonal: v_min² ≤ WR[i,i] ≤ v_max²
    @inbounds for i in 1:n
        WR_min[i, i] = v_min[i]^2
        WR_max[i, i] = v_max[i]^2
    end

    # off-diagonal bounds via voltage-product bounds for (f_bus, t_bus) pairs
    @inbounds for l in active
        bf = br.from[l]
        bt = br.to[l]
        idx_i = id_to_idx[bf]
        idx_j = id_to_idx[bt]
        amin = br.ang_min[l]
        amax = br.ang_max[l]
        vi_min, vi_max = v_min[idx_i], v_max[idx_i]
        vj_min, vj_max = v_min[idx_j], v_max[idx_j]

        if amin >= 0
            wr_up = vi_max * vj_max * cos(amin)
            wr_lo = vi_min * vj_min * cos(amax)
            wi_up = vi_max * vj_max * sin(amax)
            wi_lo = vi_min * vj_min * sin(amin)
        elseif amax <= 0
            wr_up = vi_max * vj_max * cos(amax)
            wr_lo = vi_min * vj_min * cos(amin)
            wi_up = vi_min * vj_min * sin(amax)
            wi_lo = vi_max * vj_max * sin(amin)
        else
            wr_up = vi_max * vj_max * 1.0
            wr_lo = vi_min * vj_min * min(cos(amin), cos(amax))
            wi_up = vi_max * vj_max * sin(amax)
            wi_lo = vi_max * vj_max * sin(amin)
        end

        # Intersect (tighter bounds when multiple branches connect the same pair)
        WR_max[idx_i, idx_j] = min(WR_max[idx_i, idx_j], wr_up)
        WR_min[idx_i, idx_j] = max(WR_min[idx_i, idx_j], wr_lo)
        WI_max[idx_i, idx_j] = min(WI_max[idx_i, idx_j], wi_up)
        WI_min[idx_i, idx_j] = max(WI_min[idx_i, idx_j], wi_lo)

        # symmetric / antisymmetric mirror
        WR_max[idx_j, idx_i] = WR_max[idx_i, idx_j]
        WR_min[idx_j, idx_i] = WR_min[idx_i, idx_j]
        WI_max[idx_j, idx_i] = -WI_min[idx_i, idx_j]
        WI_min[idx_j, idx_i] = -WI_max[idx_i, idx_j]
    end

    return WR_min, WR_max, WI_min, WI_max
end

function add_sdpwrm_bus_variables!(model::JuMP.Model, bus::BusData, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    n = bus.n
    bus_ids = bus.id
    id_to_idx = bus.id_to_idx

    WR_mats = Vector{Matrix{JuMP.VariableRef}}(undef, T)
    WI_mats = Vector{Matrix{JuMP.VariableRef}}(undef, T)

    WR_min, WR_max, WI_min, WI_max = _sdpwrm_buspair_bounds(bus, br)

    for t in 1:T
        WR_t = JuMP.@variable(model, [1:n, 1:n], Symmetric, base_name = "WR_sdp_$t")
        WI_t = JuMP.@variable(model, [1:n, 1:n], base_name = "WI_sdp_$t")

        # PSD constraint — JuMP enforces matrix symmetry, which in turn forces
        # WI_t to be antisymmetric (M[i,j] == M[j,i] across the -WI/WI blocks).
        JuMP.@constraint(model, [WR_t WI_t; -WI_t WR_t] in JuMP.PSDCone())

        # Identity warm-start on WR_t
        for i in 1:n, j in 1:n
            JuMP.set_start_value(WR_t[i, j], i == j ? 1.0 : 0.0)
            JuMP.set_start_value(WI_t[i, j], 0.0)
        end

        # Diagonal bounds
        for i in 1:n
            JuMP.set_lower_bound(WR_t[i, i], WR_min[i, i])
            JuMP.set_upper_bound(WR_t[i, i], WR_max[i, i])
        end

        # Off-diagonal bounds only for buspairs with a branch (finite bounds)
        for i in 1:n, j in 1:n
            i == j && continue
            if isfinite(WR_min[i, j])
                JuMP.set_lower_bound(WR_t[i, j], WR_min[i, j])
            end
            if isfinite(WR_max[i, j])
                JuMP.set_upper_bound(WR_t[i, j], WR_max[i, j])
            end
            if isfinite(WI_min[i, j])
                JuMP.set_lower_bound(WI_t[i, j], WI_min[i, j])
            end
            if isfinite(WI_max[i, j])
                JuMP.set_upper_bound(WI_t[i, j], WI_max[i, j])
            end
        end

        WR_mats[t] = Matrix(WR_t)
        WI_mats[t] = WI_t
    end

    model.ext[:WR_sdp] = WR_mats
    model.ext[:WI_sdp] = WI_mats

    # Register v_magnitude_sq[b,t] as an expression aliasing WR_sdp[t][idx,idx]
    JuMP.@expression(model, v_magnitude_sq[b in bus_ids, t in 1:T],
        WR_mats[t][id_to_idx[b], id_to_idx[b]]
    )

    return nothing
end

function add_sdpwrm_branch_variables!(model::JuMP.Model, bus::BusData, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active
    id_to_idx = bus.id_to_idx
    WR_mats = model.ext[:WR_sdp]::Vector{Matrix{JuMP.VariableRef}}
    WI_mats = model.ext[:WI_sdp]::Vector{Matrix{JuMP.VariableRef}}

    # Expose R[l,t] = WR[idx_from, idx_to, t] and T_var[l,t] = WI[idx_from, idx_to, t]
    # so add_socwr_branch_expressions! can be reused as-is.
    JuMP.@expression(model, R[l in active, t in 1:T],
        WR_mats[t][id_to_idx[br.from[l]], id_to_idx[br.to[l]]]
    )
    JuMP.@expression(model, T_var[l in active, t in 1:T],
        WI_mats[t][id_to_idx[br.from[l]], id_to_idx[br.to[l]]]
    )

    return nothing
end

function add_sdpwrm_branch_constraints!(model::JuMP.Model, bus::BusData, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    # Thermal SOC limits on (p_from, q_from) and (p_to, q_to)
    @constraint(model, [l in active, t in 1:T],
        [
            br.rate_a[l],
            model[:p_flow_branches][l, t],
            model[:q_flow_branches][l, t],
        ] in JuMP.SecondOrderCone()
    )

    @constraint(model, [l in active, t in 1:T],
        [
            br.rate_a[l],
            model[:p_reverse_flow_branches][l, t],
            model[:q_reverse_flow_branches][l, t],
        ] in JuMP.SecondOrderCone()
    )

    # tan(angmin)*wr ≤ wi ≤ tan(angmax)*wr
    add_rt_tangent_angle_cuts!(model, br)

    # Lifted Nonlinear Cuts (LNCs) — Coffrin, Hijazi, Van Hentenryck 2015
    # (arXiv:1512.04644). Strengthen the SDP relaxation; PM adds these for all
    # W-models via constraint_voltage_angle_difference(pm::AbstractWModels,...).
    add_sdpwrm_lnc_cuts!(model, bus, br)

    return nothing
end

function add_sdpwrm_lnc_cuts!(model::JuMP.Model, bus::BusData, br::BranchData)
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

            JuMP.@constraint(model,
                sf * st * (cos(phi) * wr + sin(phi) * wi) -
                vt_ub * cos(d) * st * wf -
                vf_ub * cos(d) * sf * wt >=
                vf_ub * vt_ub * cos(d) * (vf_lb * vt_lb - vf_ub * vt_ub)
            )
            JuMP.@constraint(model,
                sf * st * (cos(phi) * wr + sin(phi) * wi) -
                vt_lb * cos(d) * st * wf -
                vf_lb * cos(d) * sf * wt >=
                -vf_lb * vt_lb * cos(d) * (vf_lb * vt_lb - vf_ub * vt_ub)
            )
        end
    end

    return nothing
end
