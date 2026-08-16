@inline function _apply_dual_step!(state::LambdaState, a3::Float64)
    copyto!(state.lambda_prev, state.lambda)
    @inbounds for k in eachindex(state.lambda)
        state.lambda[k] += a3 * (state.nu[k] - state.lambda[k])
    end
    return nothing
end

function dual_update!(ws::AbstractALADINWorkspace)
    a3 = ws.alpha3
    for state in values(ws.vars)
        _apply_dual_step!(state, a3)
    end
    return nothing
end
