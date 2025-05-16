# ---------------------------------------------------------------------
# Getting Bus Parameters (preprocessed)
# ---------------------------------------------------------------------
function get_bus_parameters(bus_data)
    bus_ids = bus_data["bus_i"]
    reference_bus_idx = findfirst(==(3), bus_data["bus_type"])
    reference_bus = bus_ids[reference_bus_idx]

    return (
        bus_ids = bus_ids,
        reference_bus = reference_bus,
        vmin = Dict(bus_ids[i] => bus_data["vmin"][i] for i in eachindex(bus_ids)),
        vmax = Dict(bus_ids[i] => bus_data["vmax"][i] for i in eachindex(bus_ids)),
        vstart = Dict(bus_ids[i] => bus_data["vm"][i] for i in eachindex(bus_ids)),
    )
end

# ---------------------------------------------------------------------
# Creating Bus Variables (Vectorized and Optimized)
# ---------------------------------------------------------------------
function create_voltage_angle_variables!(model, params)
    T = get(model.ext, :time_periods, 1)
    @variable(model, theta[b in params.bus_ids, t = 1:T], start = 0.0)
    @constraint(model, theta_ref[t in 1:T], theta[params.reference_bus, t] == 0)
    return nothing
end

function create_voltage_magnitude_variables!(model, params)
    T = get(model.ext, :time_periods, 1)
    @variable(
        model,
        v_magnitude[b in params.bus_ids, t = 1:T],
        lower_bound = params.vmin[b],
        upper_bound = params.vmax[b],
        start = isnan(params.vstart[b]) ? nothing : params.vstart[b]
    )
    return nothing
end

function create_socp_variables!(model, params)
    T = get(model.ext, :time_periods, 1)

    create_voltage_magnitude_variables!(model, params)

    @variable(model, v_magnitude_sq[b in params.bus_ids, t = 1:T])
    @constraint(model, voltage_sq_def[b in params.bus_ids, t in 1:T],
        v_magnitude_sq[b, t] == model[:v_magnitude][b, t]^2
    )
    return nothing
end

# ---------------------------------------------------------------------
# Dispatching Bus Variables by Model Type
# ---------------------------------------------------------------------
function buses!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    params = get_bus_parameters(matpower_data["bus"])
    create_voltage_angle_variables!(model, params)
    return nothing
end

function buses!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = get_bus_parameters(matpower_data["bus"])
    create_voltage_angle_variables!(model, params)
    create_voltage_magnitude_variables!(model, params)
    return nothing
end

function buses!(model::JuMP.Model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    params = get_bus_parameters(matpower_data["bus"])
    create_voltage_angle_variables!(model, params)
    create_socp_variables!(model, params)
    return nothing
end
