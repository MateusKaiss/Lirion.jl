function dual_update!(
    ws::ACStateWorkspaceALADIN,
    nu::@NamedTuple{
        nu_theta::Vector{Float64},
        nu_vm::Vector{Float64},
    },
)
    copyto!(ws.lambda_theta_prev, ws.lambda_theta)
    copyto!(ws.lambda_vm_prev, ws.lambda_vm)
    a3 = ws.alpha3
    @inbounds for k in eachindex(ws.lambda_theta)
        ws.lambda_theta[k] += a3 * (nu.nu_theta[k] - ws.lambda_theta[k])
        ws.lambda_vm[k] += a3 * (nu.nu_vm[k] - ws.lambda_vm[k])
    end
    return nothing
end

function dual_update!(
    ws::ACPowerWorkspaceALADIN,
    nu::@NamedTuple{
        nu_p::Vector{Float64},
        nu_q::Vector{Float64},
    },
)
    copyto!(ws.lambda_p_prev, ws.lambda_p)
    copyto!(ws.lambda_q_prev, ws.lambda_q)
    a3 = ws.alpha3
    @inbounds for k in eachindex(ws.lambda_p)
        ws.lambda_p[k] += a3 * (nu.nu_p[k] - ws.lambda_p[k])
        ws.lambda_q[k] += a3 * (nu.nu_q[k] - ws.lambda_q[k])
    end
    return nothing
end

function dual_update!(
    ws::ACHybridWorkspaceALADIN,
    nu::@NamedTuple{
        nu_theta::Vector{Float64},
        nu_vm::Vector{Float64},
        nu_p::Vector{Float64},
        nu_q::Vector{Float64},
    },
)
    copyto!(ws.lambda_theta_prev, ws.lambda_theta)
    copyto!(ws.lambda_vm_prev, ws.lambda_vm)
    copyto!(ws.lambda_p_prev, ws.lambda_p)
    copyto!(ws.lambda_q_prev, ws.lambda_q)
    a3 = ws.alpha3
    @inbounds for k in eachindex(ws.lambda_theta)
        ws.lambda_theta[k] += a3 * (nu.nu_theta[k] - ws.lambda_theta[k])
        ws.lambda_vm[k] += a3 * (nu.nu_vm[k] - ws.lambda_vm[k])
    end
    @inbounds for k in eachindex(ws.lambda_p)
        ws.lambda_p[k] += a3 * (nu.nu_p[k] - ws.lambda_p[k])
        ws.lambda_q[k] += a3 * (nu.nu_q[k] - ws.lambda_q[k])
    end
    return nothing
end

# ------------------------------------------------------------
# Helpers for residuals and tolerances
# ------------------------------------------------------------

@inline function _ac_pf_at(zone::AladinZoneDataAC, le::Int)
    θf = zone.y[zone.edge_from_theta_col[le]]
    θt = zone.y[zone.edge_to_theta_col[le]]
    vf = zone.y[zone.edge_from_vm_col[le]]
    vt = zone.y[zone.edge_to_vm_col[le]]
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    tau = zone.edge_inv_tap[le]
    tap_angle = zone.edge_tap_angle[le]
    Δθ = θf - θt - tap_angle
    return g * (tau * vf)^2 - tau * vf * vt * (g * cos(Δθ) + b * sin(Δθ))
end

@inline function _ac_qf_at(zone::AladinZoneDataAC, le::Int)
    θf = zone.y[zone.edge_from_theta_col[le]]
    θt = zone.y[zone.edge_to_theta_col[le]]
    vf = zone.y[zone.edge_from_vm_col[le]]
    vt = zone.y[zone.edge_to_vm_col[le]]
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    b_sh = zone.edge_b_sh[le]
    tau = zone.edge_inv_tap[le]
    tap_angle = zone.edge_tap_angle[le]
    Δθ = θf - θt - tap_angle
    return -(b + b_sh / 2) * (tau * vf)^2 -
           tau * vf * vt * (-b * cos(Δθ) + g * sin(Δθ))
end

function _progress_sq_ac(zone::AladinZoneDataAC)::Float64
    s = 0.0
    y = zone.y
    x = zone.x_prev
    @inbounds for k in 1:zone.n_vars
        d = y[k] - x[k]
        s += d * d
    end
    return s
end

# ------------------------------------------------------------
# Residuals (HFD 2016 Theorem 4.1)
# ------------------------------------------------------------

function compute_residuals!(ws::ACStateWorkspaceALADIN)
    r_pri_sq = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        θi = zi.y[pair.theta_col_i]
        θj = zj.y[pair.theta_col_j]
        vi = zi.y[pair.vm_col_i]
        vj = zj.y[pair.vm_col_j]
        r_pri_sq += (θi - θj)^2 + (vi - vj)^2
    end
    ws.r_primal = sqrt(r_pri_sq)

    prog_sq = 0.0
    @inbounds for zone in ws.zones
        prog_sq += _progress_sq_ac(zone)
    end
    ws.r_progress = sqrt(prog_sq)

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.lambda_theta)
        dt = ws.lambda_theta[k] - ws.lambda_theta_prev[k]
        dv = ws.lambda_vm[k] - ws.lambda_vm_prev[k]
        lam_sq += dt * dt + dv * dv
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

function compute_residuals!(ws::ACPowerWorkspaceALADIN)
    r_pri_sq = 0.0
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _ac_pf_at(zi, pair.edge_local_i)
        pfj = _ac_pf_at(zj, pair.edge_local_j)
        qfi = _ac_qf_at(zi, pair.edge_local_i)
        qfj = _ac_qf_at(zj, pair.edge_local_j)
        r_pri_sq += (pfi - pfj)^2 + (qfi - qfj)^2
    end
    ws.r_primal = sqrt(r_pri_sq)

    prog_sq = 0.0
    @inbounds for zone in ws.zones
        prog_sq += _progress_sq_ac(zone)
    end
    ws.r_progress = sqrt(prog_sq)

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.lambda_p)
        dp = ws.lambda_p[k] - ws.lambda_p_prev[k]
        dq = ws.lambda_q[k] - ws.lambda_q_prev[k]
        lam_sq += dp * dp + dq * dq
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

function compute_residuals!(ws::ACHybridWorkspaceALADIN)
    r_pri_sq = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        θi = zi.y[pair.theta_col_i]
        θj = zj.y[pair.theta_col_j]
        vi = zi.y[pair.vm_col_i]
        vj = zj.y[pair.vm_col_j]
        r_pri_sq += (θi - θj)^2 + (vi - vj)^2
    end
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _ac_pf_at(zi, pair.edge_local_i)
        pfj = _ac_pf_at(zj, pair.edge_local_j)
        qfi = _ac_qf_at(zi, pair.edge_local_i)
        qfj = _ac_qf_at(zj, pair.edge_local_j)
        r_pri_sq += (pfi - pfj)^2 + (qfi - qfj)^2
    end
    ws.r_primal = sqrt(r_pri_sq)

    prog_sq = 0.0
    @inbounds for zone in ws.zones
        prog_sq += _progress_sq_ac(zone)
    end
    ws.r_progress = sqrt(prog_sq)

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.lambda_theta)
        dt = ws.lambda_theta[k] - ws.lambda_theta_prev[k]
        dv = ws.lambda_vm[k] - ws.lambda_vm_prev[k]
        lam_sq += dt * dt + dv * dv
    end
    @inbounds for k in eachindex(ws.lambda_p)
        dp = ws.lambda_p[k] - ws.lambda_p_prev[k]
        dq = ws.lambda_q[k] - ws.lambda_q_prev[k]
        lam_sq += dp * dp + dq * dq
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

# ------------------------------------------------------------
# Tolerances (Boyd-style)
# ------------------------------------------------------------

function _progress_tol_ac(
    ws::AbstractALADINWorkspace,
    eps_abs::Float64,
    eps_rel::Float64,
)::Float64
    n = 0
    ny_sq = 0.0
    nx_sq = 0.0
    for zone in ws.zones
        n += zone.n_vars
        @inbounds for k in 1:zone.n_vars
            ny_sq += zone.y[k]^2
            nx_sq += zone.x_prev[k]^2
        end
    end
    return sqrt(max(n, 1)) * eps_abs + eps_rel * max(sqrt(ny_sq), sqrt(nx_sq))
end

function compute_tolerances!(
    ws::ACStateWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 2 * length(ws.bus_pairs)
    norm_i = 0.0
    norm_j = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        θi = zi.y[pair.theta_col_i]
        θj = zj.y[pair.theta_col_j]
        vi = zi.y[pair.vm_col_i]
        vj = zj.y[pair.vm_col_j]
        norm_i += θi * θi + vi * vi
        norm_j += θj * θj + vj * vj
    end
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i), sqrt(norm_j))
    ws.eps_progress = _progress_tol_ac(ws, eps_abs, eps_rel)
    return nothing
end

function compute_tolerances!(
    ws::ACPowerWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 2 * length(ws.edge_pairs)
    norm_i = 0.0
    norm_j = 0.0
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _ac_pf_at(zi, pair.edge_local_i)
        pfj = _ac_pf_at(zj, pair.edge_local_j)
        qfi = _ac_qf_at(zi, pair.edge_local_i)
        qfj = _ac_qf_at(zj, pair.edge_local_j)
        norm_i += pfi * pfi + qfi * qfi
        norm_j += pfj * pfj + qfj * qfj
    end
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i), sqrt(norm_j))
    ws.eps_progress = _progress_tol_ac(ws, eps_abs, eps_rel)
    return nothing
end

function compute_tolerances!(
    ws::ACHybridWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 2 * length(ws.bus_pairs) + 2 * length(ws.edge_pairs)
    norm_i = 0.0
    norm_j = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        θi = zi.y[pair.theta_col_i]
        θj = zj.y[pair.theta_col_j]
        vi = zi.y[pair.vm_col_i]
        vj = zj.y[pair.vm_col_j]
        norm_i += θi * θi + vi * vi
        norm_j += θj * θj + vj * vj
    end
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _ac_pf_at(zi, pair.edge_local_i)
        pfj = _ac_pf_at(zj, pair.edge_local_j)
        qfi = _ac_qf_at(zi, pair.edge_local_i)
        qfj = _ac_qf_at(zj, pair.edge_local_j)
        norm_i += pfi * pfi + qfi * qfi
        norm_j += pfj * pfj + qfj * qfj
    end
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i), sqrt(norm_j))
    ws.eps_progress = _progress_tol_ac(ws, eps_abs, eps_rel)
    return nothing
end

# ------------------------------------------------------------
# Adaptive ρ/μ
# ------------------------------------------------------------

function adapt_rho_mu!(
    ws::Union{
        ACStateWorkspaceALADIN,
        ACPowerWorkspaceALADIN,
        ACHybridWorkspaceALADIN,
    },
    params::ALADINParams,
)
    params.adaptive || return nothing
    _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    return nothing
end

# ------------------------------------------------------------
# Sensitivity refresh hook (AC dispatch)
# ------------------------------------------------------------

@inline function refresh_lagrangian_hessian!(
    zone::AladinZoneDataAC,
    ws::Union{
        ACStateWorkspaceALADIN,
        ACPowerWorkspaceALADIN,
        ACHybridWorkspaceALADIN,
    },
)
    refresh_sensitivities_ac!(
        zone;
        rho = ws.rho[],
        sigma_theta = ws.sigma_theta,
        sigma_vm = ws.sigma_vm,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )
    return nothing
end

# ------------------------------------------------------------
# Logging tags
# ------------------------------------------------------------

_aladin_print_tag(::ACStateWorkspaceALADIN) = "[ALADIN-state/AC]"
_aladin_print_tag(::ACPowerWorkspaceALADIN) = "[ALADIN-power/AC]"
_aladin_print_tag(::ACHybridWorkspaceALADIN) = "[ALADIN-hybrid/AC]"
