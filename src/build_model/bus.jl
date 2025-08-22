# ---------------------------------------------------------------------
# Creating Bus Variables (Vectorized and Optimized)
# ---------------------------------------------------------------------
function create_voltage_angle_variables!(model, params)
    T = get(model.ext, :time_periods, 1)
    @variable(model, theta[b in params.bus_ids, t = 1:T], start = 0.0)
    if !isnothing(params.reference_bus) && params.reference_bus in params.bus_ids
        @constraint(model, theta_ref[t in 1:T], theta[params.reference_bus, t] == 0.0)
        # else
        #     local_ref_bus = first(params.bus_ids)
        #     @constraint(model, theta_local_ref[t in 1:T], theta[local_ref_bus, t] == 0.0)
        #     @info "Using local bar $(local_ref_bus) as auxiliary reference (not the global reference)."
    end

    return nothing
end

function create_voltage_magnitude_variables!(model, params)
    T = get(model.ext, :time_periods, 1)
    @variable(
        model,
        v_magnitude[b in params.bus_ids, t = 1:T],
        lower_bound = params.vmin[b],
        upper_bound = params.vmax[b],
        start = params.vstart[b]
    )
    return nothing
end

function create_socp_bus_variables!(model, bus_params)
    T = get(model.ext, :time_periods, 1)
    @variable(model,
        v_magnitude_sq[b in bus_params.bus_ids, t in 1:T],
        lower_bound = bus_params.vmin[b]^2,
        upper_bound = bus_params.vmax[b]^2,
        start = isnothing(bus_params.vstart[b]) ? 1.0 : (bus_params.vstart[b]^2)
    )
    return nothing
end

# ---------------------------------------------------------------------
# Dispatching Bus Variables by Model Type
# ---------------------------------------------------------------------
function buses!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    params = matpower_data["bus_params"]
    create_voltage_angle_variables!(model, params)
    return nothing
end

function buses!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = matpower_data["bus_params"]
    create_voltage_angle_variables!(model, params)
    create_voltage_magnitude_variables!(model, params)
    return nothing
end

function buses!(model::JuMP.Model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    bus_params = matpower_data["bus_params"]
    create_voltage_angle_variables!(model, bus_params)
    create_socp_bus_variables!(model, bus_params)
    return nothing
end
