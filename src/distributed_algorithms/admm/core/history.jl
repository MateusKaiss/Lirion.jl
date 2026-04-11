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
    for (k, v) in _get_rhos(ws)
        _push_rho!(ws.hist, k, v)
    end
    return nothing
end
