# ---------------------------------------------------------------------
# Getting Reactive Power Balance Parameters
# ---------------------------------------------------------------------
function get_reactive_power_parameters(matpower_data)
    return (
        bus_ids = matpower_data["bus"]["bus_i"],
        generator_bus_index = matpower_data["gen"]["gen_bus"],
        battery_bus_index = get(matpower_data, "storage", Dict("storage_bus" => Int[]))["storage_bus"],
        prosumer_bus_index = get(matpower_data, "prosumer", Dict("bus" => Int[]))["bus"],
        demand_bus_index = matpower_data["bus"]["bus_i"],
        branch_bus_from = matpower_data["branch"]["f_bus"],
        branch_bus_to = matpower_data["branch"]["t_bus"],
        shunt_susceptance = matpower_data["bus"]["bs"],
        generator_status = matpower_data["gen"]["gen_status"],
        has_timeseries = haskey(matpower_data, "bus_timeseries"),
    )
end

function get_reactive_demand_at_time(data::Dict{String, Any}, bus_index::Int, t::Int)
    if haskey(data, "bus_timeseries")
        df = data["bus_timeseries"]
        return sum(
            row.qd for row in eachrow(df) if row.bus == bus_index && row.time == t;
            init = 0.0,
        )
    else
        return sum(
            data["bus"]["qd"][d] for
            (d, b) in enumerate(data["bus"]["bus_i"]) if b == bus_index;
            init = 0.0,
        )
    end
end

function compute_reactive_power_balance(model, matpower_data, bus_index, t)
    params = get_reactive_power_parameters(matpower_data)

    return sum(
        model[:q_generation][g, t] for (g, b) in enumerate(params.generator_bus_index)
        if b == bus_index && params.generator_status[g] == 1; init = 0.0,
    ) -
           sum(
        model[:q_battery_injected][bat, t] for
        (bat, b) in enumerate(params.battery_bus_index)
        if b == bus_index; init = 0.0,
    ) +
           sum(
        model[:q_prosumer_injected][p, t] for
        (p, b) in enumerate(params.prosumer_bus_index)
        if b == bus_index; init = 0.0,
    ) -
           get_reactive_demand_at_time(matpower_data, bus_index, t) -
           sum(
        model[:q_flow_branches][l, t] for (l, f) in enumerate(params.branch_bus_from)
        if f == bus_index; init = 0.0,
    ) -
           sum(
        model[:q_reverse_flow_branches][l, t] for
        (l, to) in enumerate(params.branch_bus_to)
        if to == bus_index; init = 0.0,
    )
end

function reactive_power_balance!(model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    return nothing
end

function reactive_power_balance!(model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = get_reactive_power_parameters(matpower_data)
    T = get(model.ext, :time_periods, 1)

    for (idx, b) in enumerate(params.bus_ids), t in 1:T
        balance = compute_reactive_power_balance(model, matpower_data, b, t)
        @constraint(
            model,
            balance + params.shunt_susceptance[idx] * model[:v_magnitude][b, t]^2 == 0
        )
    end
end

function reactive_power_balance!(
    model,
    matpower_data,
    ::Type{SOCP},
    ::Type{SubproblemBuild},
)
    params = get_reactive_power_parameters(matpower_data)
    T = get(model.ext, :time_periods, 1)

    for (idx, b) in enumerate(params.bus_ids), t in 1:T
        balance = compute_reactive_power_balance(model, matpower_data, b, t)
        @constraint(
            model,
            balance + params.shunt_susceptance[idx] * model[:v_magnitude_sq][b, t] == 0
        )
    end
end
