@inline function _ac_pf_at(zone::AladinZoneDataAC, le::Int)
    θf = zone.y[zone.edge_from_theta_col[le]]
    θt = zone.y[zone.edge_to_theta_col[le]]
    vf = zone.y[zone.edge_from_vm_col[le]]
    vt = zone.y[zone.edge_to_vm_col[le]]
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    tau = zone.edge_inv_tap[le]
    tap_angle = zone.edge_tap_angle[le]
    Δθ = θf - θt - tap_angle
    return g * (tau * vf)^2 - tau * vf * vt * (g * cos(Δθ) + b * sin(Δθ))
end

@inline function _ac_qf_at(zone::AladinZoneDataAC, le::Int)
    θf = zone.y[zone.edge_from_theta_col[le]]
    θt = zone.y[zone.edge_to_theta_col[le]]
    vf = zone.y[zone.edge_from_vm_col[le]]
    vt = zone.y[zone.edge_to_vm_col[le]]
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    b_sh = zone.edge_b_sh[le]
    tau = zone.edge_inv_tap[le]
    tap_angle = zone.edge_tap_angle[le]
    Δθ = θf - θt - tap_angle
    return -(b + b_sh / 2) * (tau * vf)^2 -
           tau * vf * vt * (-b * cos(Δθ) + g * sin(Δθ))
end

@inline function _socwr_flow_value_at_y(
    zone::AladinZoneDataSOCWR,
    le::Int,
    which::Symbol,
)::Float64
    cols = zone.edge_flow_cols[le]
    coefs = which === :p ? zone.edge_pflow_coefs[le] : zone.edge_qflow_coefs[le]
    const_val = which === :p ? zone.edge_pflow_const[le] : zone.edge_qflow_const[le]
    acc = const_val
    y = zone.y
    @inbounds for i in eachindex(cols)
        acc += coefs[i] * y[cols[i]]
    end
    return acc
end
