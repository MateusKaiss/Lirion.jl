# ============================================================
# SAVE — PREV
# ============================================================

# ============================================================
# DC
# ============================================================

function save_prev!(ws::DCStateWorkspace)
    copyto!(ws.vars.theta.consensus_prev, ws.vars.theta.consensus)
    return nothing
end

function save_prev!(ws::DCPowerWorkspace)
    copyto!(ws.vars.p_flow.consensus_prev, ws.vars.p_flow.consensus)
    return nothing
end

function save_prev!(ws::DCHybridWorkspace)
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
    return nothing
end

function save_prev!(ws::ACHybridWorkspace)
    copyto!(ws.vars.theta.consensus_prev, ws.vars.theta.consensus)
    copyto!(ws.vars.v_mag.consensus_prev, ws.vars.v_mag.consensus)
    copyto!(ws.vars.p_flow.consensus_prev, ws.vars.p_flow.consensus)
    copyto!(ws.vars.q_flow.consensus_prev, ws.vars.q_flow.consensus)
    if !isempty(ws.bus_pairs)
        copyto!(ws.vars.p_inj.consensus_prev, ws.vars.p_inj.consensus)
        copyto!(ws.vars.q_inj.consensus_prev, ws.vars.q_inj.consensus)
    end
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
    return nothing
end

function save_prev!(ws::SOCWRHybridWorkspace)
    copyto!(ws.vars.v2.consensus_prev, ws.vars.v2.consensus)
    copyto!(ws.vars.p_flow.consensus_prev, ws.vars.p_flow.consensus)
    copyto!(ws.vars.q_flow.consensus_prev, ws.vars.q_flow.consensus)
    copyto!(ws.vars.r.consensus_prev, ws.vars.r.consensus)
    copyto!(ws.vars.t.consensus_prev, ws.vars.t.consensus)
    return nothing
end
