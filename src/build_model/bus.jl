# ---------------------------------------------------------------------
# Getting Bus Parameters (preprocessed)
# ---------------------------------------------------------------------
function get_bus_parameters(bus_data)
    bus_ids = bus_data["bus_i"]
    reference_bus_idx = findfirst(==(3), bus_data["bus_type"])
    reference_bus = bus_ids[reference_bus_idx]
    vmin = Dict{Int, Float64}()
    vmax = Dict{Int, Float64}()
    vstart = Dict{Int, Union{Nothing, Float64}}()
    @inbounds for i in eachindex(bus_ids)
        b = bus_ids[i]
        vmin[b] = bus_data["vmin"][i]
        vmax[b] = bus_data["vmax"][i]
        vstart[b] = isnan(bus_data["vm"][i]) ? nothing : bus_data["vm"][i]
    end

    return (
        bus_ids = bus_ids,
        reference_bus = reference_bus,
        vmin = vmin,
        vmax = vmax,
        vstart = vstart,
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
        start = 1.001
    )

    scale = select_cone_scale(params)

    @constraint(model, voltage_sq_def[b in params.bus_ids, t in 1:T],
        [sqrt(scale), model[:v_magnitude][b, t] * sqrt(scale), model[:v_magnitude_sq][b, t] * sqrt(scale)]
        in JuMP.RotatedSecondOrderCone()
    )
    
    # Alternativas (não usadas, mas deixadas para referência):
    # Igualdade exata (costuma deixar o modelo mais lento e menos estável)
    # @constraint(model, voltage_sq_def[b in params.bus_ids, t in 1:T],
    #     v_magnitude_sq[b, t] == model[:v_magnitude][b, t]^2
    # )

    # Desigualdade relaxada (mais rápida, mas fisicamente imprecisa)
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

function select_cone_scale(params)
    vmin_vals = collect(values(params.vmin))
    vmax_vals = collect(values(params.vmax))

    vmin = clamp(minimum(vmin_vals), 0.80, 1.05)
    vmax = clamp(maximum(vmax_vals), 0.90, 1.15)
    delta = vmax / vmin

    if delta < 1.11
        base_scale = 1.0 / max(vmin^2, 1e-4)
        correction = 1.0 + 0.5 * log(delta)^2
        raw_scale = base_scale * correction
        return clamp(raw_scale, 25.0, 150.0)
    else
        return 100.0
    end
end

function select_cone_scale(params)
    vmin_vals = collect(values(params.vmin))
    vmax_vals = collect(values(params.vmax))

    vmin = clamp(minimum(vmin_vals), 0.80, 1.05)
    vmax = clamp(maximum(vmax_vals), 0.90, 1.15)
    delta = vmax / vmin

    if delta < 1.11
        return 30
    elseif delta < 1.13
        return 20
    else
        return 100.0
    end
end