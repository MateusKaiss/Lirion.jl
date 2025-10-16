# ---------------------------------------------------------------------
# Creating Bus Variables (Vectorized and Optimized)
# ---------------------------------------------------------------------
function create_voltage_angle_variables!(model, params)
    T = get(model.ext, :time_periods, 1)

    @variable(model, theta[b in params.bus_ids, t = 1:T])

    if !isnothing(params.reference_bus) && params.reference_bus in params.bus_ids
        @constraint(model, [t in 1:T], theta[params.reference_bus, t] == 0.0)
    end

    ref = params.reference_bus
    ref_va =
        isnothing(ref) ? 0.0 :
        (
            params.initial_voltage_angle[ref] === nothing ? 0.0 :
            params.initial_voltage_angle[ref]
        )

    for b in params.bus_ids, t in 1:T
        va_b = params.initial_voltage_angle[b]
        if va_b !== nothing && !isnan(va_b)
            JuMP.set_start_value(theta[b, t], va_b - ref_va)
        else
            JuMP.set_start_value(theta[b, t], 0.0)
        end
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
