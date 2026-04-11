# ============================================================
# SAVE — PREV (consensus snapshot for dual residual)
# ============================================================

# ============================================================
# DC
# ============================================================

function save_prev!(ws::DCStateWorkspace)
    ws.theta_consensus_prev .= ws.theta_consensus
    return nothing
end

function save_prev!(ws::DCPowerWorkspace)
    for k in ws.interzonal
        ws.p_flow_consensus_prev[k] = ws.p_flow_consensus[k]
    end
    return nothing
end

function save_prev!(ws::DCHybridWorkspace)
    ws.theta_consensus_prev .= ws.theta_consensus
    for k in ws.interzonal
        ws.p_flow_consensus_prev[k] = ws.p_flow_consensus[k]
    end
    return nothing
end

# ============================================================
# AC
# ============================================================

function save_prev!(ws::ACStateWorkspace)
    ws.theta_consensus_prev .= ws.theta_consensus
    ws.v_mag_consensus_prev .= ws.v_mag_consensus
    return nothing
end

function save_prev!(ws::ACPowerWorkspace)
    for k in ws.interzonal
        ws.p_flow_consensus_prev[k] = ws.p_flow_consensus[k]
        ws.q_flow_consensus_prev[k] = ws.q_flow_consensus[k]
    end
    return nothing
end

function save_prev!(ws::ACHybridWorkspace)
    ws.theta_consensus_prev .= ws.theta_consensus
    ws.v_mag_consensus_prev .= ws.v_mag_consensus

    for k in ws.interzonal
        ws.p_flow_consensus_prev[k] = ws.p_flow_consensus[k]
        ws.q_flow_consensus_prev[k] = ws.q_flow_consensus[k]
    end
    return nothing
end

# ============================================================
# SOCP
# ============================================================

function save_prev!(ws::SOCPStateWorkspace)
    ws.v2_consensus_prev .= ws.v2_consensus

    for k in ws.interzonal
        ws.r_consensus_prev[k] = ws.r_consensus[k]
        ws.t_consensus_prev[k] = ws.t_consensus[k]
    end
    return nothing
end

function save_prev!(ws::SOCPPowerWorkspace)
    for k in ws.interzonal
        ws.p_flow_consensus_prev[k] = ws.p_flow_consensus[k]
        ws.q_flow_consensus_prev[k] = ws.q_flow_consensus[k]
    end
    return nothing
end

function save_prev!(ws::SOCPHybridWorkspace)
    ws.v2_consensus_prev .= ws.v2_consensus

    for k in ws.interzonal
        ws.p_flow_consensus_prev[k] = ws.p_flow_consensus[k]
        ws.q_flow_consensus_prev[k] = ws.q_flow_consensus[k]
        ws.r_consensus_prev[k] = ws.r_consensus[k]
        ws.t_consensus_prev[k] = ws.t_consensus[k]
    end
    return nothing
end
