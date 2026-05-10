# ============================================================
# DC
# ============================================================

function save_prev!(ws::DCStateWorkspace)
    copyto!(ws.vars.theta.consensus_prev, ws.vars.theta.consensus)
    return nothing
end

function save_prev!(ws::DCPowerWorkspace)
    copyto!(ws.vars.p_flow.consensus_prev, ws.vars.p_flow.consensus)
    copyto!(ws.vars.p_dummy.consensus_prev, ws.vars.p_dummy.consensus)
    return nothing
end

function save_prev!(ws::DCHybridWorkspace{<:Union{TieLine, DummyBus}})
    copyto!(ws.vars.theta.consensus_prev, ws.vars.theta.consensus)
    copyto!(ws.vars.p_flow.consensus_prev, ws.vars.p_flow.consensus)
    copyto!(ws.vars.p_dummy.consensus_prev, ws.vars.p_dummy.consensus)
    return nothing
end

function save_prev!(ws::DCHybridWorkspace{SplitBus})
    copyto!(ws.vars.theta.consensus_prev, ws.vars.theta.consensus)
    copyto!(ws.vars.p_flow.consensus_prev, ws.vars.p_flow.consensus)
    copyto!(ws.vars.p_inj.consensus_prev, ws.vars.p_inj.consensus)
    return nothing
end

# ============================================================
# AC
# ============================================================

function save_prev!(ws::ACStateWorkspace)
    copyto!(ws.vars.theta.consensus_prev, ws.vars.theta.consensus)
    copyto!(ws.vars.v_mag.consensus_prev, ws.vars.v_mag.consensus)
    return nothing
end

function save_prev!(ws::ACPowerWorkspace)
    copyto!(ws.vars.p_flow.consensus_prev, ws.vars.p_flow.consensus)
    copyto!(ws.vars.q_flow.consensus_prev, ws.vars.q_flow.consensus)
    copyto!(ws.vars.p_dummy.consensus_prev, ws.vars.p_dummy.consensus)
    copyto!(ws.vars.q_dummy.consensus_prev, ws.vars.q_dummy.consensus)
    return nothing
end

function save_prev!(ws::ACHybridWorkspace{<:Union{TieLine, DummyBus}})
    copyto!(ws.vars.theta.consensus_prev, ws.vars.theta.consensus)
    copyto!(ws.vars.v_mag.consensus_prev, ws.vars.v_mag.consensus)
    copyto!(ws.vars.p_flow.consensus_prev, ws.vars.p_flow.consensus)
    copyto!(ws.vars.q_flow.consensus_prev, ws.vars.q_flow.consensus)
    copyto!(ws.vars.p_dummy.consensus_prev, ws.vars.p_dummy.consensus)
    copyto!(ws.vars.q_dummy.consensus_prev, ws.vars.q_dummy.consensus)
    return nothing
end

function save_prev!(ws::ACHybridWorkspace{SplitBus})
    copyto!(ws.vars.theta.consensus_prev, ws.vars.theta.consensus)
    copyto!(ws.vars.v_mag.consensus_prev, ws.vars.v_mag.consensus)
    copyto!(ws.vars.p_flow.consensus_prev, ws.vars.p_flow.consensus)
    copyto!(ws.vars.q_flow.consensus_prev, ws.vars.q_flow.consensus)
    copyto!(ws.vars.p_inj.consensus_prev, ws.vars.p_inj.consensus)
    copyto!(ws.vars.q_inj.consensus_prev, ws.vars.q_inj.consensus)
    return nothing
end

# ============================================================
# SOCWR
# ============================================================

function save_prev!(ws::SOCWRStateWorkspace)
    copyto!(ws.vars.v2.consensus_prev, ws.vars.v2.consensus)
    copyto!(ws.vars.r.consensus_prev, ws.vars.r.consensus)
    copyto!(ws.vars.t.consensus_prev, ws.vars.t.consensus)
    return nothing
end

function save_prev!(ws::SOCWRPowerWorkspace)
    copyto!(ws.vars.p_flow.consensus_prev, ws.vars.p_flow.consensus)
    copyto!(ws.vars.q_flow.consensus_prev, ws.vars.q_flow.consensus)
    copyto!(ws.vars.p_dummy.consensus_prev, ws.vars.p_dummy.consensus)
    copyto!(ws.vars.q_dummy.consensus_prev, ws.vars.q_dummy.consensus)
    return nothing
end

function save_prev!(ws::SOCWRHybridWorkspace{<:Union{TieLine, DummyBus}})
    copyto!(ws.vars.v2.consensus_prev, ws.vars.v2.consensus)
    copyto!(ws.vars.p_flow.consensus_prev, ws.vars.p_flow.consensus)
    copyto!(ws.vars.q_flow.consensus_prev, ws.vars.q_flow.consensus)
    copyto!(ws.vars.p_dummy.consensus_prev, ws.vars.p_dummy.consensus)
    copyto!(ws.vars.q_dummy.consensus_prev, ws.vars.q_dummy.consensus)
    copyto!(ws.vars.r.consensus_prev, ws.vars.r.consensus)
    copyto!(ws.vars.t.consensus_prev, ws.vars.t.consensus)
    return nothing
end

function save_prev!(ws::SOCWRHybridWorkspace{SplitBus})
    copyto!(ws.vars.v2.consensus_prev, ws.vars.v2.consensus)
    copyto!(ws.vars.p_flow.consensus_prev, ws.vars.p_flow.consensus)
    copyto!(ws.vars.q_flow.consensus_prev, ws.vars.q_flow.consensus)
    copyto!(ws.vars.r.consensus_prev, ws.vars.r.consensus)
    copyto!(ws.vars.t.consensus_prev, ws.vars.t.consensus)
    copyto!(ws.vars.p_inj.consensus_prev, ws.vars.p_inj.consensus)
    copyto!(ws.vars.q_inj.consensus_prev, ws.vars.q_inj.consensus)
    return nothing
end
