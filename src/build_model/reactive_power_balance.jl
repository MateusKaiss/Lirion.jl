# ---------------------------------------------------------------------
# Getting Reactive Power Balance Parameters (with preprocessed maps)
# ---------------------------------------------------------------------
function get_reactive_power_parameters(matpower_data)
    return (
        bus_ids = matpower_data["bus"]["bus_i"],
        shunt_susceptance = matpower_data["bus"]["bs"],
        generator_status = matpower_data["gen"]["gen_status"],
        has_timeseries = haskey(matpower_data, "bus_timeseries"),
        bus_idx_map = matpower_data["bus_idx_map"],
        bus_to_generators = get(matpower_data, "bus_to_gens", Dict{Int, Vector{Int}}()),
        bus_to_batteries = get(matpower_data, "bus_to_batteries", Dict{Int, Vector{Int}}()),
        bus_to_prosumers = get(matpower_data, "bus_to_prosumers", Dict{Int, Vector{Int}}()),
        bus_to_branches_from = get(
            matpower_data,
            "bus_to_branches_from",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_branches_to = get(
            matpower_data,
            "bus_to_branches_to",
            Dict{Int, Vector{Int}}(),
        ),
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
            data["bus"]["qd"][i] for i in eachindex(data["bus"]["bus_i"])
            if data["bus"]["bus_i"][i] == bus_index;
            init = 0.0,
        )
    end
end

function compute_reactive_power_balance(model, matpower_data, params, bus_index, t)
    q_gen = sum(
        model[:q_generation][g, t] for
        g in get(params.bus_to_generators, bus_index, Int[])
        if params.generator_status[g] == 1;
        init = 0.0,
    )
    q_bat = sum(
        model[:q_battery_injected][b, t] for
        b in get(params.bus_to_batteries, bus_index, Int[]);
        init = 0.0,
    )
    q_pros = sum(
        model[:q_prosumer_injected][p, t] for
        p in get(params.bus_to_prosumers, bus_index, Int[]);
        init = 0.0,
    )
    q_out = sum(
        model[:q_flow_branches][l, t] for
        l in get(params.bus_to_branches_from, bus_index, Int[]);
        init = 0.0,
    )
    q_in = sum(
        model[:q_reverse_flow_branches][l, t] for
        l in get(params.bus_to_branches_to, bus_index, Int[]);
        init = 0.0,
    )
    qd = get_reactive_demand_at_time(matpower_data, bus_index, t)

    return q_gen - q_bat + q_pros - qd - q_out - q_in
end

function reactive_power_balance!(model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    return nothing
end

function reactive_power_balance!(model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = get_reactive_power_parameters(matpower_data)
    T = get(model.ext, :time_periods, 1)

    for (idx, b) in enumerate(params.bus_ids), t in 1:T
        balance = compute_reactive_power_balance(model, matpower_data, params, b, t)
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
        balance = compute_reactive_power_balance(model, matpower_data, params, b, t)
        @constraint(
            model,
            balance + params.shunt_susceptance[idx] * model[:v_magnitude_sq][b, t] == 0
        )
    end
end
