# ============================================================
# SOCWR storage constraints
# ============================================================

function add_socwr_storage_constraints!(model::JuMP.Model, stor::StorageData)
    T = get(model.ext, :time_periods, 1)
    B = stor.n

    @variable(model, z_b[b in 1:B, t in 1:T], Bin)

    @variable(
        model,
        ccms[b in 1:B, t in 1:T],
        lower_bound = 0.0,
        start = 0.0,
    )

    vmin_global = get(model.ext, :v_min, 0.9)

    @inbounds for b in 1:B, t in 1:T
        ub = (stor.thermal_limit[b] / vmin_global)^2
        set_upper_bound(ccms[b, t], ub)
    end

    @constraint(model, [b in 1:B, t in 1:T],
        model[:charge][b, t] <= z_b[b, t] * stor.charge_ramp[b]
    )

    @constraint(model, [b in 1:B, t in 1:T],
        model[:discharge][b, t] <= (1 - z_b[b, t]) * stor.discharge_ramp[b]
    )

    _storage_energy_balance!(model, stor)

    @constraint(model, [b in 1:B, t in 1:T],
        model[:p_storage_injected][b, t]^2 +
        model[:q_storage_injected][b, t]^2 <=
        model[:v_magnitude_sq][stor.bus[b], t] * ccms[b, t]
    )

    @constraint(model, [b in 1:B, t in 1:T],
        model[:p_storage_injected][b, t] +
        (model[:discharge][b, t] - model[:charge][b, t]) ==
        stor.p_loss[b] + stor.r[b] * ccms[b, t]
    )

    @constraint(model, [b in 1:B, t in 1:T],
        model[:q_storage_injected][b, t] ==
        model[:qsc][b, t] + stor.q_loss[b] + stor.x[b] * ccms[b, t]
    )

    return nothing
end
