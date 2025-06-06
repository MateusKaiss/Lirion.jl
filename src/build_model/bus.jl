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
        start = params.vstart[b]
    )
    return nothing
end

function create_socp_variables!(model, params)
    T = get(model.ext, :time_periods, 1)

    create_voltage_magnitude_variables!(model, params)

    @variable(model,
        v_magnitude_sq[b in params.bus_ids, t in 1:T],
        lower_bound = params.vmin[b]^2,
        upper_bound = params.vmax[b]^2,
        start = params.vstart[b]^2
    )

    scale = get(model.ext, :selected_scale, nothing)
    scale = isnothing(scale) ? select_cone_scale(params) : scale

    @constraint(model, voltage_sq_def[b in params.bus_ids, t in 1:T],
        [
            sqrt(scale),
            model[:v_magnitude][b, t] * sqrt(scale),
            model[:v_magnitude_sq][b, t] * sqrt(scale),
        ]
        in
        JuMP.RotatedSecondOrderCone()
    )

    # Alternatives (not used, but kept for reference):
    # Strict equality (usually slower and less numerically stable)
    # @constraint(model, voltage_sq_def[b in params.bus_ids, t in 1:T],
    #     v_magnitude_sq[b, t] == model[:v_magnitude][b, t]^2
    # )

    # Relaxed inequality (faster but physically inaccurate)
    # @constraint(model, voltage_sq_def[b in params.bus_ids, t in 1:T],
    #     v_magnitude_sq[b, t] >= model[:v_magnitude][b, t]^2
    # )

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
    params = matpower_data["bus_params"]
    create_voltage_angle_variables!(model, params)
    create_socp_variables!(model, params)
    return nothing
end