# ---------------------------------------------------------------------
# Creating Bus Variables
# ---------------------------------------------------------------------
function create_voltage_angle_variables!(model, bus_data)
    number_of_buses = length(bus_data["bus_i"])
    reference_bus_index = findfirst(b -> b == 3, bus_data["bus_type"])
    T = get(model.ext, :time_periods, 1)

    @variable(model, theta[b = 1:number_of_buses, t = 1:T], start = 0.0)
    @constraint(model, [t = 1:T], theta[reference_bus_index, t] == 0)

    return nothing
end

function create_voltage_magnitude_variables!(model, bus_data)
    number_of_buses = length(bus_data["bus_i"])
    T = get(model.ext, :time_periods, 1)

    @variable(
        model,
        v_magnitude[b = 1:number_of_buses, t = 1:T],
        lower_bound = bus_data["vmin"][b],
        upper_bound = bus_data["vmax"][b],
        start = isnan(bus_data["vm"][b]) ? nothing : bus_data["vm"][b]
    )

    return nothing
end

function create_socp_variables!(model, bus_data)
    number_of_buses = length(bus_data["bus_i"])
    T = get(model.ext, :time_periods, 1)

    create_voltage_magnitude_variables!(model, bus_data)

    @variable(model, v_magnitude_sq[b = 1:number_of_buses, t = 1:T],
    lower_bound = bus_data["vmin"][b]^2,
    upper_bound = bus_data["vmax"][b]^2,
    )
    @constraint(
        model,
        [b = 1:number_of_buses, t = 1:T],
        v_magnitude_sq[b, t] >= model[:v_magnitude][b, t]^2
    )

    return nothing
end

# ---------------------------------------------------------------------
# Dispatching Bus Variables by Model Type
# ---------------------------------------------------------------------
function buses!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    bus_data = matpower_data["bus"]
    create_voltage_angle_variables!(model, bus_data)
    return nothing
end

function buses!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    bus_data = matpower_data["bus"]
    create_voltage_angle_variables!(model, bus_data)
    create_voltage_magnitude_variables!(model, bus_data)
    return nothing
end

function buses!(model::JuMP.Model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    bus_data = matpower_data["bus"]
    create_voltage_angle_variables!(model, bus_data)
    create_socp_variables!(model, bus_data)
    return nothing
end