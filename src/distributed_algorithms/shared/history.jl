# ============================================================
# ADMM
# ============================================================

@inline function _push_rho!(hist::ADMMHistory, key::Symbol, val::Real)
    v = get!(hist.rhos, key, Float64[])
    push!(v, Float64(val))
    return nothing
end

@inline function _record_iter!(ws::AbstractADMMWorkspace)
    push!(ws.hist.r_primal, ws.r_primal)
    push!(ws.hist.r_dual, ws.r_dual)
    push!(ws.hist.eps_primal, ws.eps_primal)
    push!(ws.hist.eps_dual, ws.eps_dual)
    for (k, v) in pairs(_get_rhos(ws))
        _push_rho!(ws.hist, k, v)
    end
    return nothing
end

# ============================================================
# ALADIN
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
