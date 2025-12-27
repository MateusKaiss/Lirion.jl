# ============================================================
# Adaptive rho + refreshers (ADMM)
# ============================================================
# ------------------------------------------------------------
# AC refreshers: rebuild cached quadratic terms when rho changes
# ------------------------------------------------------------

@inline function refresh_rho_bus_ac!(
    partitions,
    rho_theta::Base.RefValue{Float64},
    prev_rho_theta::Base.RefValue{Float64},
    rho_v::Base.RefValue{Float64},
    prev_rho_v::Base.RefValue{Float64},
)
    rho_theta_val = rho_theta[]
    rho_v_val = rho_v[]

    if (rho_theta_val != prev_rho_theta[]) || (rho_v_val != prev_rho_v[])
        for part in values(partitions)
            _build_quad_bus_expr_ac!(part, rho_theta_val, rho_v_val)
        end
        prev_rho_theta[] = rho_theta_val
        prev_rho_v[] = rho_v_val
    end
    return nothing
end

@inline function refresh_rho_branch_ac!(
    partitions,
    rho_p::Base.RefValue{Float64},
    prev_rho_p::Base.RefValue{Float64},
    rho_q::Base.RefValue{Float64},
    prev_rho_q::Base.RefValue{Float64},
)
    rho_p_val = rho_p[]
    rho_q_val = rho_q[]

    if (rho_p_val != prev_rho_p[]) || (rho_q_val != prev_rho_q[])
        for part in values(partitions)
            _build_quad_branch_expr_ac!(part, rho_p_val, rho_q_val, part.zone_edges)
        end
        prev_rho_p[] = rho_p_val
        prev_rho_q[] = rho_q_val
    end
    return nothing
end

# ------------------------------------------------------------
# DC refreshers: rebuild cached quadratic terms when rho changes
# ------------------------------------------------------------

@inline function refresh_rho_bus!(
    partitions,
    rho::Base.RefValue{Float64},
    prev_rho::Base.RefValue{Float64},
)
    rho_val = rho[]
    if rho_val != prev_rho[]
        for part in values(partitions)
            _build_quad_bus_expr!(part, rho_val)
        end
        prev_rho[] = rho_val
    end
    return nothing
end

@inline function refresh_rho_branch!(
    partitions,
    rho::Base.RefValue{Float64},
    prev_rho::Base.RefValue{Float64},
)
    rho_val = rho[]
    if rho_val != prev_rho[]
        for part in values(partitions)
            _build_quad_branch_expr!(part, rho_val, part.zone_edges)
        end
        prev_rho[] = rho_val
    end
    return nothing
end

# ------------------------------------------------------------
# SOCP refreshers: rebuild cached objective pieces when rho changes
# ------------------------------------------------------------

@inline function refresh_rho_v2_socp!(
    partitions,
    rhoV2::Base.RefValue{Float64},
    prev_rhoV2::Base.RefValue{Float64},
)
    rho_v2_val = rhoV2[]
    if rho_v2_val != prev_rhoV2[]
        for sub in values(partitions)
            _build_quad_bus_expr_socp!(sub, rho_v2_val)
            sub.obj_static = sub.obj_base + sub.quad_bus_expr_socp + sub.quad_rt_expr_socp
        end
        prev_rhoV2[] = rho_v2_val
    end
    return nothing
end

@inline function refresh_rho_rt_socp!(
    partitions,
    rhoR::Base.RefValue{Float64},
    prev_rhoR::Base.RefValue{Float64},
    rhoT::Base.RefValue{Float64},
    prev_rhoT::Base.RefValue{Float64},
)
    rho_r_val = rhoR[]
    rho_t_val = rhoT[]

    if (rho_r_val != prev_rhoR[]) || (rho_t_val != prev_rhoT[])
        for sub in values(partitions)
            _build_quad_branch_rt_expr!(sub, rho_r_val, rho_t_val, sub.zone_edges)
            sub.obj_static = sub.obj_base + sub.quad_bus_expr_socp + sub.quad_rt_expr_socp
        end
        prev_rhoR[] = rho_r_val
        prev_rhoT[] = rho_t_val
    end
    return nothing
end

@inline function refresh_rho_pq_socp!(
    partitions,
    rhoP::Base.RefValue{Float64},
    prev_rhoP::Base.RefValue{Float64},
    rhoQ::Base.RefValue{Float64},
    prev_rhoQ::Base.RefValue{Float64},
)
    rho_p_val = rhoP[]
    rho_q_val = rhoQ[]

    if (rho_p_val != prev_rhoP[]) || (rho_q_val != prev_rhoQ[])
        for sub in values(partitions)
            _build_quad_branch_pq_expr!(sub, rho_p_val, rho_q_val, sub.zone_edges)
            sub.obj_static = sub.obj_base + sub.quad_pq_expr_socp + sub.soft_anchor_expr
        end
        prev_rhoP[] = rho_p_val
        prev_rhoQ[] = rho_q_val
    end
    return nothing
end

@inline function refresh_rho_hybrid_socp!(
    partitions,
    rhoV2::Base.RefValue{Float64},
    prev_rhoV2::Base.RefValue{Float64},
    rhoP::Base.RefValue{Float64},
    prev_rhoP::Base.RefValue{Float64},
    rhoQ::Base.RefValue{Float64},
    prev_rhoQ::Base.RefValue{Float64},
    rhoR::Base.RefValue{Float64},
    prev_rhoR::Base.RefValue{Float64},
    rhoT::Base.RefValue{Float64},
    prev_rhoT::Base.RefValue{Float64},
)
    rho_v2_val = rhoV2[]
    rho_p_val = rhoP[]
    rho_q_val = rhoQ[]
    rho_r_val = rhoR[]
    rho_t_val = rhoT[]

    if (rho_v2_val != prev_rhoV2[]) || (rho_p_val != prev_rhoP[]) ||
       (rho_q_val != prev_rhoQ[]) ||
       (rho_r_val != prev_rhoR[]) || (rho_t_val != prev_rhoT[])
        for sub in values(partitions)
            (rho_v2_val != prev_rhoV2[]) && _build_quad_bus_expr_socp!(sub, rho_v2_val)

            ((rho_p_val != prev_rhoP[]) || (rho_q_val != prev_rhoQ[])) &&
                _build_quad_branch_pq_expr!(sub, rho_p_val, rho_q_val, sub.zone_edges)

            ((rho_r_val != prev_rhoR[]) || (rho_t_val != prev_rhoT[])) &&
                _build_quad_branch_rt_expr!(sub, rho_r_val, rho_t_val, sub.zone_edges)

            sub.obj_static =
                sub.obj_base + sub.quad_bus_expr_socp + sub.quad_pq_expr_socp +
                sub.quad_rt_expr_socp
        end

        prev_rhoV2[] = rho_v2_val
        prev_rhoP[] = rho_p_val
        prev_rhoQ[] = rho_q_val
        prev_rhoR[] = rho_r_val
        prev_rhoT[] = rho_t_val
    end
    return nothing
end

# ------------------------------------------------------------
# Core: adaptive rho update rule
# ------------------------------------------------------------

# ---------- rho adaptativo (paridade com DC) ----------
function adaptive_rho!(
    rho::Base.RefValue{Float64},
    r::Float64,
    s::Float64;
    mu_incr = 2.0,
    mu_decr = 2.0,
    tau_incr = 1.05,
    tau_decr = 0.98,
)
    if r > mu_incr * s
        rho[] *= tau_incr
    elseif s > mu_decr * r
        rho[] *= tau_decr
    end
end

# ------------------------------------------------------------
# Adaptive rho kernels: AC / DC
# ------------------------------------------------------------

function adaptive_rho_ac!(
    rho::Base.RefValue{Float64},
    r::Float64,
    s::Float64;
    mu_incr = 2.0,
    mu_decr = 2.0,
    tau_incr = 1.05,
    tau_decr = 0.98,
)
    if r > mu_incr * s
        rho[] *= tau_incr
    elseif s > mu_decr * r
        rho[] *= tau_decr
    end
end

function adaptive_rho_dc!(
    rho::Base.RefValue{Float64},
    r::Float64,
    s::Float64;
    mu_incr = 3.0,
    mu_decr = 3.0,
    tau_incr = 1.05,
    tau_decr = 0.98,
)
    if r > mu_incr * s
        rho[] *= tau_incr
    elseif s > mu_decr * r
        rho[] *= tau_decr
    end
end

# ------------------------------------------------------------
# AC adaptive rho: bus (θ,V via components) + branch (P,Q)
# ------------------------------------------------------------

# --- ADAPTIVE RHO (bus AC: θ e V separados via componentes) ---
@inline function maybe_adapt_rho!(
    adaptive::Bool,
    rho_theta::Base.RefValue{Float64},
    rho_v::Base.RefValue{Float64},
    partitions,
    theta_consensus,
    vmag_consensus,
    theta_consensus_prev,
    vmag_consensus_prev,
    shared_buses_sorted,
    shared_buses_idx,
    config,
)
    adaptive || return nothing

    rθ, rV = primal_residual_bus_ac_components(
        partitions,
        theta_consensus,
        vmag_consensus,
        shared_buses_sorted,
        shared_buses_idx,
    )
    sθ, sV = dual_residual_bus_ac_components(
        theta_consensus,
        vmag_consensus,
        theta_consensus_prev,
        vmag_consensus_prev,
    )

    adaptive_rho_ac!(
        rho_theta, rθ, sθ;
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    adaptive_rho_ac!(
        rho_v, rV, sV;
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    return nothing
end

# --- ADAPTIVE RHO (branch AC/SOCP-style: usa r,s do branch pros dois) ---
@inline function maybe_adapt_rho_socp!(
    adaptive::Bool,
    rho_p::Base.RefValue{Float64},
    rho_q::Base.RefValue{Float64},
    r::Float64,
    s::Float64,
    config,
)
    adaptive || return nothing
    adaptive_rho_ac!(
        rho_p, r, s;
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    adaptive_rho_ac!(
        rho_q, r, s;
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    return nothing
end

# --- ADAPTIVE RHO (hybrid AC: θ,V por componentes + P,Q via r_br,s_br) ---
@inline function maybe_adapt_rho!(
    adaptive::Bool,
    rho_theta::Base.RefValue{Float64},
    rho_v::Base.RefValue{Float64},
    partitions,
    theta_consensus,
    vmag_consensus,
    theta_consensus_prev,
    vmag_consensus_prev,
    shared_buses_sorted,
    shared_buses_idx,
    rho_p::Base.RefValue{Float64},
    rho_q::Base.RefValue{Float64},
    r_br::Float64,
    s_br::Float64,
    config,
)
    adaptive || return nothing

    maybe_adapt_rho!(
        true,
        rho_theta,
        rho_v,
        partitions,
        theta_consensus,
        vmag_consensus,
        theta_consensus_prev,
        vmag_consensus_prev,
        shared_buses_sorted,
        shared_buses_idx,
        config,
    )
    maybe_adapt_rho_socp!(true, rho_p, rho_q, r_br, s_br, config)

    return nothing
end

# ------------------------------------------------------------
# DC helpers: stop + maybe_adapt_rho overloads
# ------------------------------------------------------------

@inline function stop_if_converged!(r::Float64, s::Float64, eps_r::Float64, eps_s::Float64)
    return (r < eps_r) && (s < eps_s)
end

@inline function maybe_adapt_rho!(
    adaptive::Bool,
    rho::Base.RefValue{Float64},
    r::Float64,
    s::Float64,
    config::DistributedOPFModel,
)
    adaptive || return nothing
    adaptive_rho_dc!(
        rho, r, s;
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    return nothing
end

@inline function maybe_adapt_rho!(
    adaptive::Bool,
    rho_theta::Base.RefValue{Float64},
    r_bus::Float64,
    s_bus::Float64,
    rho_flow::Base.RefValue{Float64},
    r_br::Float64,
    s_br::Float64,
    config::DistributedOPFModel,
)
    adaptive || return nothing
    adaptive_rho_dc!(
        rho_theta, r_bus, s_bus;
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    adaptive_rho_dc!(
        rho_flow, r_br, s_br;
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    return nothing
end

# ------------------------------------------------------------
# Helpers: ADMMHistory rho tracking
# ------------------------------------------------------------

@inline function _push_rho!(hist::ADMMHistory, key::Symbol, val::Real)
    v = get!(hist.rhos, key, Float64[])
    push!(v, Float64(val))
    return nothing
end

# ------------------------------------------------------------
# SOCP adaptive rho: overloads
# ------------------------------------------------------------

@inline function maybe_adapt_rho!(
    adaptive::Bool,
    rhoV2::Base.RefValue{Float64},
    r_bus::Float64,
    s_bus::Float64,
    rhoR::Base.RefValue{Float64},
    rhoT::Base.RefValue{Float64},
    r_rt::Float64,
    s_rt::Float64,
    config,
)
    adaptive || return nothing
    adaptive_rho_dc!(
        rhoV2, r_bus, s_bus;
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    adaptive_rho_dc!(
        rhoR, r_rt / sqrt(2), s_rt / sqrt(2);
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    adaptive_rho_dc!(
        rhoT, r_rt / sqrt(2), s_rt / sqrt(2);
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    return nothing
end

@inline function maybe_adapt_rho!(
    adaptive::Bool,
    rhoP::Base.RefValue{Float64},
    rhoQ::Base.RefValue{Float64},
    r::Float64,
    s::Float64,
    config,
)
    adaptive || return nothing
    adaptive_rho_dc!(
        rhoP, r / sqrt(2), s / sqrt(2);
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    adaptive_rho_dc!(
        rhoQ, r / sqrt(2), s / sqrt(2);
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    return nothing
end

@inline function maybe_adapt_rho!(
    adaptive::Bool,
    rhoV2::Base.RefValue{Float64},
    r_bus::Float64,
    s_bus::Float64,
    rhoP::Base.RefValue{Float64},
    rhoQ::Base.RefValue{Float64},
    r_pq::Float64,
    s_pq::Float64,
    rhoR::Base.RefValue{Float64},
    rhoT::Base.RefValue{Float64},
    r_rt::Float64,
    s_rt::Float64,
    config,
)
    adaptive || return nothing

    adaptive_rho_dc!(
        rhoV2, r_bus, s_bus;
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )

    adaptive_rho_dc!(
        rhoP, r_pq / sqrt(2), s_pq / sqrt(2);
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    adaptive_rho_dc!(
        rhoQ, r_pq / sqrt(2), s_pq / sqrt(2);
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )

    adaptive_rho_dc!(
        rhoR, r_rt / sqrt(2), s_rt / sqrt(2);
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )
    adaptive_rho_dc!(
        rhoT, r_rt / sqrt(2), s_rt / sqrt(2);
        mu_incr = config.mu_incr,
        mu_decr = config.mu_decr,
        tau_incr = config.tau_incr,
        tau_decr = config.tau_decr,
    )

    return nothing
end
