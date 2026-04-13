# ============================================================
# ALADIN — Phase 2 dual update (line search, HFD Eq. 13b)
# ------------------------------------------------------------
#
#   λᵏ⁺¹ = λᵏ + α₃ (ν - λᵏ)
#
# where ν is the Lagrange multiplier of the consensus row in the
# coordinator QP. With α₃ = 1 this collapses to λ ← ν.
# ============================================================

function dual_update!(ws::DCStateWorkspaceALADIN, nu::Vector{Float64})
    copyto!(ws.lambda_bus_prev, ws.lambda_bus)
    a3 = ws.alpha3
    @inbounds for k in eachindex(ws.lambda_bus)
        ws.lambda_bus[k] = ws.lambda_bus[k] + a3 * (nu[k] - ws.lambda_bus[k])
    end
    return nothing
end

function dual_update!(ws::DCPowerWorkspaceALADIN, nu::Vector{Float64})
    copyto!(ws.lambda_edge_prev, ws.lambda_edge)
    a3 = ws.alpha3
    @inbounds for k in eachindex(ws.lambda_edge)
        ws.lambda_edge[k] = ws.lambda_edge[k] + a3 * (nu[k] - ws.lambda_edge[k])
    end
    return nothing
end

function dual_update!(
    ws::DCHybridWorkspaceALADIN,
    nu::@NamedTuple{nu_bus::Vector{Float64}, nu_edge::Vector{Float64}}
)
    copyto!(ws.lambda_bus_prev, ws.lambda_bus)
    copyto!(ws.lambda_edge_prev, ws.lambda_edge)
    a3 = ws.alpha3
    @inbounds for k in eachindex(ws.lambda_bus)
        ws.lambda_bus[k] = ws.lambda_bus[k] + a3 * (nu.nu_bus[k] - ws.lambda_bus[k])
    end
    @inbounds for k in eachindex(ws.lambda_edge)
        ws.lambda_edge[k] = ws.lambda_edge[k] + a3 * (nu.nu_edge[k] - ws.lambda_edge[k])
    end
    return nothing
end
