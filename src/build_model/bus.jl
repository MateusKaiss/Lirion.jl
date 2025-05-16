# ---------------------------------------------------------------------
# Creating Bus Variables (Vectorized and Optimized)
# ---------------------------------------------------------------------
function create_voltage_angle_variables!(model, bus_data)
    bus_ids = bus_data["bus_i"]
    T = get(model.ext, :time_periods, 1)
    reference_bus_idx = findfirst(==(3), bus_data["bus_type"])
    reference_bus = bus_ids[reference_bus_idx]

    @variable(model, theta[b in bus_ids, t = 1:T], start = 0.0)
    @constraint(model, theta_ref[t = 1:T], theta[reference_bus, t] == 0)

    return nothing
end

function create_voltage_magnitude_variables!(model, bus_data)
    bus_ids = bus_data["bus_i"]
    T = get(model.ext, :time_periods, 1)
    vmin = Dict(bus_ids[i] => bus_data["vmin"][i] for i in eachindex(bus_ids))
    vmax = Dict(bus_ids[i] => bus_data["vmax"][i] for i in eachindex(bus_ids))
    vstart = Dict(
        bus_ids[i] => isnan(bus_data["vm"][i]) ? nothing : bus_data["vm"][i] for
        i in eachindex(bus_ids)
    )

    @variable(
        model,
        v_magnitude[b in bus_ids, t = 1:T],
        lower_bound = vmin[b],
        upper_bound = vmax[b],
        start = vstart[b]
    )

    return nothing
end

function create_socp_variables!(model, bus_data)
    bus_ids = bus_data["bus_i"]
    T = get(model.ext, :time_periods, 1)

    create_voltage_magnitude_variables!(model, bus_data)

    @variable(model, v_magnitude_sq[b in bus_ids, t = 1:T])
    @constraint(model, voltage_sq_def[b in bus_ids, t in 1:T],
        v_magnitude_sq[b, t] == model[:v_magnitude][b, t]^2
    )

    return nothing
end

# ---------------------------------------------------------------------
# Dispatching Bus Variables by Model Type
# ---------------------------------------------------------------------
function buses!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    create_voltage_angle_variables!(model, matpower_data["bus"])
    return nothing
end

function buses!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    create_voltage_angle_variables!(model, matpower_data["bus"])
    create_voltage_magnitude_variables!(model, matpower_data["bus"])
    return nothing
end

function buses!(model::JuMP.Model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    create_voltage_angle_variables!(model, matpower_data["bus"])
    create_socp_variables!(model, matpower_data["bus"])
    return nothing
end
