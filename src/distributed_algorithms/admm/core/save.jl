# ============================================================
# SAVE — PREV (consensus snapshot for dual residual)
# ============================================================

# ============================================================
# DC
# ============================================================

function save_prev!(ws::DCStateWorkspace)
    copyto!(ws.theta_consensus_prev, ws.theta_consensus)
    return nothing
end

function save_prev!(ws::DCPowerWorkspace)
    copyto!(ws.p_flow_consensus_prev, ws.p_flow_consensus)
    return nothing
end

function save_prev!(ws::DCHybridWorkspace)
    copyto!(ws.theta_consensus_prev, ws.theta_consensus)
    copyto!(ws.p_flow_consensus_prev, ws.p_flow_consensus)
    return nothing
end

# ============================================================
# AC
# ============================================================

function save_prev!(ws::ACStateWorkspace)
    copyto!(ws.theta_consensus_prev, ws.theta_consensus)
    copyto!(ws.v_mag_consensus_prev, ws.v_mag_consensus)
    return nothing
end

function save_prev!(ws::ACPowerWorkspace)
    copyto!(ws.p_flow_consensus_prev, ws.p_flow_consensus)
    copyto!(ws.q_flow_consensus_prev, ws.q_flow_consensus)
    return nothing
end

function save_prev!(ws::ACHybridWorkspace)
    copyto!(ws.theta_consensus_prev, ws.theta_consensus)
    copyto!(ws.v_mag_consensus_prev, ws.v_mag_consensus)
    copyto!(ws.p_flow_consensus_prev, ws.p_flow_consensus)
    copyto!(ws.q_flow_consensus_prev, ws.q_flow_consensus)
    return nothing
end

# ============================================================
# SOCWR
# ============================================================

function save_prev!(ws::SOCWRStateWorkspace)
    copyto!(ws.v2_consensus_prev, ws.v2_consensus)
    copyto!(ws.r_consensus_prev, ws.r_consensus)
    copyto!(ws.t_consensus_prev, ws.t_consensus)
    return nothing
end

function save_prev!(ws::SOCWRPowerWorkspace)
    copyto!(ws.p_flow_consensus_prev, ws.p_flow_consensus)
    copyto!(ws.q_flow_consensus_prev, ws.q_flow_consensus)
    return nothing
end

function save_prev!(ws::SOCWRHybridWorkspace)
    copyto!(ws.v2_consensus_prev, ws.v2_consensus)
    copyto!(ws.p_flow_consensus_prev, ws.p_flow_consensus)
    copyto!(ws.q_flow_consensus_prev, ws.q_flow_consensus)
    copyto!(ws.r_consensus_prev, ws.r_consensus)
    copyto!(ws.t_consensus_prev, ws.t_consensus)
    return nothing
end
