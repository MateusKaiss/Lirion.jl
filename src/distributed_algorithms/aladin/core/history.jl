# ============================================================
# ALADIN — Diagnostics history (canonical HFD 2016 metrics)
# ------------------------------------------------------------
# Fields tracked per iteration:
#
#   r_primal        ‖∑_i A_i y_i‖₂          (consensus residual)
#   r_progress      ‖y_i - x_iᵏ‖₂           (Phase-1 progress / KKT proxy)
#   eps_primal      Boyd-style tolerance for r_primal
#   eps_progress    Boyd-style tolerance for r_progress
#   lambda_step     ‖λᵏ⁺¹ - λᵏ‖₂            (dual movement, diagnostic only)
#   rhos            adaptive ρ trajectory (one entry per rho key, e.g. :rho_bus)
#   mus             adaptive μ trajectory
# ============================================================

Base.@kwdef struct ALADINHistory <: AbstractRunDiagnostics
    r_primal::Vector{Float64} = Float64[]
    r_progress::Vector{Float64} = Float64[]
    eps_primal::Vector{Float64} = Float64[]
    eps_progress::Vector{Float64} = Float64[]
    lambda_step::Vector{Float64} = Float64[]
    rhos::Dict{Symbol, Vector{Float64}} = Dict{Symbol, Vector{Float64}}()
    mus::Vector{Float64} = Float64[]
    meta::Dict{Symbol, Any} = Dict{Symbol, Any}()
end

@inline function _push_rho_aladin!(hist::ALADINHistory, key::Symbol, val::Real)
    v = get!(hist.rhos, key, Float64[])
    push!(v, Float64(val))
    return nothing
end

@inline function _record_iter!(ws::AbstractALADINWorkspace)
    push!(ws.hist.r_primal, ws.r_primal)
    push!(ws.hist.r_progress, ws.r_progress)
    push!(ws.hist.eps_primal, ws.eps_primal)
    push!(ws.hist.eps_progress, ws.eps_progress)
    push!(ws.hist.lambda_step, ws.lambda_step)
    push!(ws.hist.mus, ws.mu[])
    for (k, v) in pairs(_get_rhos(ws))
        _push_rho_aladin!(ws.hist, k, v)
    end
    return nothing
end
