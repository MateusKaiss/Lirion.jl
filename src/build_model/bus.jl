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

# SOCWR (SOCP): usar w = |V|^2 (sem theta)
function create_socp_variables!(model, params)
    T = get(model.ext, :time_periods, 1)

    @variable(
        model,
        v_magnitude_sq[b in params.bus_ids, t in 1:T],
        lower_bound = params.vmin[b]^2,
        upper_bound = params.vmax[b]^2,
        start = max(params.vstart[b]^2, params.vmin[b]^2 + 1.0e-3),
        base_name = "0_w"  
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
    params = matpower_data["bus_params"]
    # SOCWR: não criar theta; somente w = |V|^2
    create_socp_variables!(model, params)
    return nothing
end