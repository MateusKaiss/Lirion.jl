# ============================================================
# ADMM Rho Adaptation & Refresh — unified for DC · AC · SOCWR
# ============================================================

# ############################################################
#  CORE: adaptive penalty parameter
# ############################################################

# Boyd et al. adaptive scheme:
#   if r > mu_incr * s  →  rho *= tau_incr
#   if s > mu_decr * r  →  rho *= tau_decr
function _adaptive_rho!(
    rho::Base.RefValue{Float64}, r::Float64, s::Float64;
    mu_incr::Float64 = 3.0, mu_decr::Float64 = 3.0,
    tau_incr::Float64 = 1.05, tau_decr::Float64 = 0.98,
)
    if r > mu_incr * s
        rho[] *= tau_incr
    elseif s > mu_decr * r
        rho[] *= tau_decr
    end
    return nothing
end

@inline function _adapt!(rho, r, s, params::ADMMParams)
    return _adaptive_rho!(rho, r, s;
        mu_incr = params.mu_incr, mu_decr = params.mu_decr,
        tau_incr = params.tau_incr, tau_decr = params.tau_decr,
    )
end

# ############################################################
#  REFRESH: rebuild quadratic terms when rho changes
# ############################################################

# ---- DC bus ----
# Returns `true` iff ρ changed (so the caller can rebuild obj_static).
function _refresh_rho_bus_dc!(parts::Partitions{PartitionDC}, rho, prev_rho)::Bool
    rho_val = rho[]
    rho_val == prev_rho[] && return false
    @inbounds for part in parts.parts
        build_quad_bus_expr!(part, rho_val)
    end
    prev_rho[] = rho_val
    return true
end

# ---- DC branch ----
function _refresh_rho_branch_dc!(parts::Partitions{PartitionDC}, rho, prev_rho)::Bool
    rho_val = rho[]
    rho_val == prev_rho[] && return false
    @inbounds for part in parts.parts
        build_quad_branch_expr!(part, rho_val, part.zone_edges)
    end
    prev_rho[] = rho_val
    return true
end

# ---- AC bus ----
function _refresh_rho_bus_ac!(
    parts::Partitions{PartitionAC},
    rho_theta,
    prev_rho_theta,
    rho_v_mag,
    prev_rho_v_mag,
)
    rho_theta_val = rho_theta[]
    rho_v_mag_val = rho_v_mag[]
    if rho_theta_val != prev_rho_theta[] || rho_v_mag_val != prev_rho_v_mag[]
        @inbounds for part in parts.parts
            build_quad_bus_expr!(part, rho_theta_val, rho_v_mag_val)
        end
        prev_rho_theta[] = rho_theta_val
        prev_rho_v_mag[] = rho_v_mag_val
    end
    return nothing
end

# ---- AC branch ----
function _refresh_rho_branch_ac!(
    parts::Partitions{PartitionAC},
    rho_p,
    prev_rho_p,
    rho_q,
    prev_rho_q,
)
    rho_p_val = rho_p[]
    rho_q_val = rho_q[]
    if rho_p_val != prev_rho_p[] || rho_q_val != prev_rho_q[]
        @inbounds for part in parts.parts
            build_quad_branch_expr!(part, rho_p_val, rho_q_val, part.zone_edges)
        end
        prev_rho_p[] = rho_p_val
        prev_rho_q[] = rho_q_val
    end
    return nothing
end

# ---- SOCWR bus (V²) ----
function _refresh_rho_bus_socwr!(parts::Partitions{PartitionSOCWR}, rho, prev_rho)
    rho_val = rho[]
    if rho_val != prev_rho[]
        @inbounds for part in parts.parts
            build_quad_bus_expr!(part, rho_val)
            part.obj_static =
                part.obj_base + part.quad_bus_expr + part.quad_rt_expr
        end
        prev_rho[] = rho_val
    end
    return nothing
end

# ---- SOCWR branch RT ----
function _refresh_rho_branch_rt_socwr!(
    parts::Partitions{PartitionSOCWR},
    rho_r,
    prev_rho_r,
    rho_t,
    prev_rho_t,
)
    rho_r_val = rho_r[]
    rho_t_val = rho_t[]
    if rho_r_val != prev_rho_r[] || rho_t_val != prev_rho_t[]
        @inbounds for part in parts.parts
            build_quad_branch_rt_expr!(part, rho_r_val, rho_t_val, part.zone_edges)
            part.obj_static =
                part.obj_base + part.quad_bus_expr + part.quad_rt_expr
        end
        prev_rho_r[] = rho_r_val
        prev_rho_t[] = rho_t_val
    end
    return nothing
end

# ---- SOCWR branch PQ ----
function _refresh_rho_branch_pq_socwr!(
    parts::Partitions{PartitionSOCWR},
    rho_p,
    prev_rho_p,
    rho_q,
    prev_rho_q,
)
    rho_p_val = rho_p[]
    rho_q_val = rho_q[]
    if rho_p_val != prev_rho_p[] || rho_q_val != prev_rho_q[]
        @inbounds for part in parts.parts
            build_quad_branch_pq_expr!(part, rho_p_val, rho_q_val, part.zone_edges)
            part.obj_static =
                part.obj_base + part.quad_pq_expr + part.soft_anchor_expr
        end
        prev_rho_p[] = rho_p_val
        prev_rho_q[] = rho_q_val
    end
    return nothing
end

# ---- SOCWR Hybrid (all groups) ----
# NOTE: unified refresh (instead of calling 3 separate helpers like DC/AC Hybrid)
# because obj_static depends on v2, pq, AND rt — doing it in one pass avoids
# rebuilding obj_static up to 3 times.
function _refresh_rho_hybrid_socwr!(
    parts,
    rho_v2, prev_rho_v2,
    rho_p, prev_rho_p,
    rho_q, prev_rho_q,
    rho_r, prev_rho_r,
    rho_t, prev_rho_t,
)
    rho_v2_val = rho_v2[]
    rho_p_val = rho_p[]
    rho_q_val = rho_q[]
    rho_r_val = rho_r[]
    rho_t_val = rho_t[]

    changed =
        rho_v2_val != prev_rho_v2[] ||
        rho_p_val != prev_rho_p[] ||
        rho_q_val != prev_rho_q[] ||
        rho_r_val != prev_rho_r[] ||
        rho_t_val != prev_rho_t[]
    changed || return nothing

    @inbounds for part in parts.parts
        rho_v2_val != prev_rho_v2[] && build_quad_bus_expr!(part, rho_v2_val)
        (rho_p_val != prev_rho_p[] || rho_q_val != prev_rho_q[]) &&
            build_quad_branch_pq_expr!(part, rho_p_val, rho_q_val, part.zone_edges)
        (rho_r_val != prev_rho_r[] || rho_t_val != prev_rho_t[]) &&
            build_quad_branch_rt_expr!(part, rho_r_val, rho_t_val, part.zone_edges)
        part.obj_static =
            part.obj_base + part.quad_bus_expr + part.quad_pq_expr +
            part.quad_rt_expr
    end

    prev_rho_v2[] = rho_v2_val
    prev_rho_p[] = rho_p_val
    prev_rho_q[] = rho_q_val
    prev_rho_r[] = rho_r_val
    prev_rho_t[] = rho_t_val
    return nothing
end

# ############################################################
#  DISPATCH: refresh_rho!(ws)
# ############################################################

# ---- DC ----
function refresh_rho!(ws::DCStateWorkspace)
    if _refresh_rho_bus_dc!(ws.partitions, ws.rho_theta, ws.prev_rho_theta)
        @inbounds for part in ws.partitions.parts
            part.obj_static = part.obj_base + part.quad_bus_expr
        end
    end
    return nothing
end

function refresh_rho!(ws::DCPowerWorkspace)
    if _refresh_rho_branch_dc!(ws.partitions, ws.rho_p, ws.prev_rho_p)
        @inbounds for part in ws.partitions.parts
            part.obj_static = part.obj_base + part.quad_branch_expr
        end
    end
    return nothing
end

function refresh_rho!(ws::DCHybridWorkspace)
    bus_changed =
        _refresh_rho_bus_dc!(ws.partitions, ws.rho_theta, ws.prev_rho_theta)
    branch_changed =
        _refresh_rho_branch_dc!(ws.partitions, ws.rho_p, ws.prev_rho_p)
    if bus_changed || branch_changed
        @inbounds for part in ws.partitions.parts
            part.obj_static =
                part.obj_base + part.quad_bus_expr + part.quad_branch_expr
        end
    end
    return nothing
end

# ---- AC ----
refresh_rho!(ws::ACStateWorkspace) = _refresh_rho_bus_ac!(
    ws.partitions,
    ws.rho_theta, ws.prev_rho_theta,
    ws.rho_v, ws.prev_rho_v,
)

refresh_rho!(ws::ACPowerWorkspace) = _refresh_rho_branch_ac!(
    ws.partitions,
    ws.rho_p, ws.prev_rho_p,
    ws.rho_q, ws.prev_rho_q,
)

function refresh_rho!(ws::ACHybridWorkspace)
    _refresh_rho_bus_ac!(
        ws.partitions,
        ws.rho_theta, ws.prev_rho_theta,
        ws.rho_v, ws.prev_rho_v,
    )
    _refresh_rho_branch_ac!(
        ws.partitions,
        ws.rho_p, ws.prev_rho_p,
        ws.rho_q, ws.prev_rho_q,
    )
    return nothing
end

# ---- SOCWR ----
function refresh_rho!(ws::SOCWRStateWorkspace)
    _refresh_rho_bus_socwr!(ws.partitions, ws.rho_v2, ws.prev_rho_v2)
    _refresh_rho_branch_rt_socwr!(
        ws.partitions,
        ws.rho_r, ws.prev_rho_r,
        ws.rho_t, ws.prev_rho_t,
    )
    return nothing
end

refresh_rho!(ws::SOCWRPowerWorkspace) = _refresh_rho_branch_pq_socwr!(
    ws.partitions,
    ws.rho_p, ws.prev_rho_p,
    ws.rho_q, ws.prev_rho_q,
)

function refresh_rho!(ws::SOCWRHybridWorkspace)
    return _refresh_rho_hybrid_socwr!(
        ws.partitions,
        ws.rho_v2, ws.prev_rho_v2,
        ws.rho_p, ws.prev_rho_p,
        ws.rho_q, ws.prev_rho_q,
        ws.rho_r, ws.prev_rho_r,
        ws.rho_t, ws.prev_rho_t,
    )
end

# Fallback for any workspace type not explicitly handled
refresh_rho!(::AbstractADMMWorkspace) = nothing

# ############################################################
#  DISPATCH: rho_update!(ws, params)
# ############################################################

# ---- DC State ----
function rho_update!(ws::DCStateWorkspace, params::ADMMParams)
    params.adaptive || return nothing
    _adapt!(ws.rho_theta, ws.r_primal, ws.r_dual, params)
    return nothing
end

# ---- DC Power ----
function rho_update!(ws::DCPowerWorkspace, params::ADMMParams)
    params.adaptive || return nothing
    _adapt!(ws.rho_p, ws.r_primal, ws.r_dual, params)
    return nothing
end

# ---- DC Hybrid (per-component adaptation) ----
function rho_update!(ws::DCHybridWorkspace, params::ADMMParams)
    params.adaptive || return nothing

    _adapt!(ws.rho_theta, ws.r_primal_bus, ws.r_dual_bus, params)
    _adapt!(ws.rho_p, ws.r_primal_branch, ws.r_dual_branch, params)
    return nothing
end

# ---- AC State (per-component theta/V adaptation) ----
function rho_update!(ws::ACStateWorkspace, params::ADMMParams)
    params.adaptive || return nothing

    _adapt!(ws.rho_theta, ws.r_primal_theta, ws.r_dual_theta, params)
    _adapt!(ws.rho_v, ws.r_primal_v_mag, ws.r_dual_v_mag, params)
    return nothing
end

# ---- AC Power (P/Q adaptation) ----
function rho_update!(ws::ACPowerWorkspace, params::ADMMParams)
    params.adaptive || return nothing

    _adapt!(ws.rho_p, ws.r_primal_p, ws.r_dual_p, params)
    _adapt!(ws.rho_q, ws.r_primal_q, ws.r_dual_q, params)
    return nothing
end

# ---- AC Hybrid ----
function rho_update!(ws::ACHybridWorkspace, params::ADMMParams)
    params.adaptive || return nothing

    _adapt!(ws.rho_theta, ws.r_primal_theta, ws.r_dual_theta, params)
    _adapt!(ws.rho_v, ws.r_primal_v_mag, ws.r_dual_v_mag, params)
    _adapt!(ws.rho_p, ws.r_primal_p, ws.r_dual_p, params)
    _adapt!(ws.rho_q, ws.r_primal_q, ws.r_dual_q, params)
    return nothing
end

# ---- SOCWR State (V² + RT) ----
function rho_update!(ws::SOCWRStateWorkspace, params::ADMMParams)
    params.adaptive || return nothing

    _adapt!(ws.rho_v2, ws.r_primal_bus, ws.r_dual_bus, params)
    _adapt!(ws.rho_r, ws.r_primal_r, ws.r_dual_r, params)
    _adapt!(ws.rho_t, ws.r_primal_t, ws.r_dual_t, params)
    return nothing
end

# ---- SOCWR Power (PQ) ----
function rho_update!(ws::SOCWRPowerWorkspace, params::ADMMParams)
    params.adaptive || return nothing

    _adapt!(ws.rho_p, ws.r_primal_p, ws.r_dual_p, params)
    _adapt!(ws.rho_q, ws.r_primal_q, ws.r_dual_q, params)
    return nothing
end

# ---- SOCWR Hybrid (V² + PQ + RT) ----
function rho_update!(ws::SOCWRHybridWorkspace, params::ADMMParams)
    params.adaptive || return nothing

    _adapt!(ws.rho_v2, ws.r_primal_bus, ws.r_dual_bus, params)
    _adapt!(ws.rho_p, ws.r_primal_p, ws.r_dual_p, params)
    _adapt!(ws.rho_q, ws.r_primal_q, ws.r_dual_q, params)
    _adapt!(ws.rho_r, ws.r_primal_r, ws.r_dual_r, params)
    _adapt!(ws.rho_t, ws.r_primal_t, ws.r_dual_t, params)
    return nothing
end
