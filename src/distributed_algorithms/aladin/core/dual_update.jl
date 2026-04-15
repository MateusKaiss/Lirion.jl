@inline function _apply_dual_step!(state::LambdaState, a3::Float64)
    copyto!(state.lambda_prev, state.lambda)
    @inbounds for k in eachindex(state.lambda)
        state.lambda[k] += a3 * (state.nu[k] - state.lambda[k])
    end
    return nothing
end

function dual_update!(ws::DCStateWorkspaceALADIN)
    _apply_dual_step!(ws.vars.theta, ws.alpha3)
    return nothing
end

function dual_update!(ws::DCPowerWorkspaceALADIN)
    _apply_dual_step!(ws.vars.p_flow, ws.alpha3)
    return nothing
end

function dual_update!(ws::DCHybridWorkspaceALADIN)
    a3 = ws.alpha3
    _apply_dual_step!(ws.vars.theta, a3)
    _apply_dual_step!(ws.vars.p_flow, a3)
    return nothing
end

function dual_update!(ws::ACStateWorkspaceALADIN)
    a3 = ws.alpha3
    _apply_dual_step!(ws.vars.theta, a3)
    _apply_dual_step!(ws.vars.v_mag, a3)
    return nothing
end

function dual_update!(ws::ACPowerWorkspaceALADIN)
    a3 = ws.alpha3
    _apply_dual_step!(ws.vars.p_flow, a3)
    _apply_dual_step!(ws.vars.q_flow, a3)
    return nothing
end

function dual_update!(ws::ACHybridWorkspaceALADIN)
    a3 = ws.alpha3
    _apply_dual_step!(ws.vars.theta, a3)
    _apply_dual_step!(ws.vars.v_mag, a3)
    _apply_dual_step!(ws.vars.p_flow, a3)
    _apply_dual_step!(ws.vars.q_flow, a3)
    return nothing
end

function dual_update!(ws::SOCWRStateWorkspaceALADIN)
    a3 = ws.alpha3
    _apply_dual_step!(ws.vars.v2, a3)
    _apply_dual_step!(ws.vars.r, a3)
    _apply_dual_step!(ws.vars.t, a3)
    return nothing
end

function dual_update!(ws::SOCWRPowerWorkspaceALADIN)
    a3 = ws.alpha3
    _apply_dual_step!(ws.vars.p_flow, a3)
    _apply_dual_step!(ws.vars.q_flow, a3)
    return nothing
end

function dual_update!(ws::SOCWRHybridWorkspaceALADIN)
    a3 = ws.alpha3
    _apply_dual_step!(ws.vars.v2, a3)
    _apply_dual_step!(ws.vars.r, a3)
    _apply_dual_step!(ws.vars.t, a3)
    _apply_dual_step!(ws.vars.p_flow, a3)
    _apply_dual_step!(ws.vars.q_flow, a3)
    return nothing
end
