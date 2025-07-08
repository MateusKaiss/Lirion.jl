function build_model(data, model_type; optimizer)
    if model_type == SOCP
        return try_with_fallback_scale(data, model_type; optimizer)
    else
        return build_and_optimize_model(data, model_type; optimizer)
    end
end

function build_and_optimize_model(
    matpower_data,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
)
    optimizer === nothing && (optimizer = Ipopt.Optimizer)
    times = get_model_times(matpower_data)
    if length(unique(times)) > 1
        return build_and_optimize_multiperiod_model(matpower_data, model_type; optimizer, times=times)
    else
        return build_and_optimize_singleperiod_model(matpower_data, model_type; optimizer)
    end
end

"""
Attempts to solve the model using a sequence of scaling factors
to improve numerical stability for conic models (e.g., SOCP).

This function tries up to 5 different scale values:
  1. A default scale selected by `select_cone_scale`
  2. A fallback value of 150
  3. A minimal scale of 1
  4. A larger fallback of 300
  5. A final extreme attempt with scale = 500

Each attempt uses a deepcopy of the original input data to avoid side effects
from previous runs. Returns the model from the last successful attempt
(or the final attempt if all fail).
"""
function try_with_fallback_scale(
    matpower_data::Dict,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
)
    # First attempt using the scale selected from bus data
    scale = select_cone_scale(matpower_data["bus_params"])
    matpower_data_attempt = deepcopy(matpower_data)
    matpower_data_attempt["selected_scale"] = scale
    model = build_and_optimize_model(matpower_data_attempt, model_type; optimizer)
    status = JuMP.termination_status(model)

    # Second attempt with scale = 150
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Model with scale=$scale failed (status: $status). Trying scale=150..."
        matpower_data_fallback = deepcopy(matpower_data)
        matpower_data_fallback["selected_scale"] = 150
        model = build_and_optimize_model(matpower_data_fallback, model_type; optimizer)
        status = JuMP.termination_status(model)
    end

    # Third attempt with scale = 1
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Model still failed with scale=150 (status: $status). Trying scale=1..."
        matpower_data_final_attempt = deepcopy(matpower_data)
        matpower_data_final_attempt["selected_scale"] = 1
        model = build_and_optimize_model(matpower_data_final_attempt, model_type; optimizer)
        status = JuMP.termination_status(model)
    end

    # Fourth attempt with scale = 300
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Model still failed with scale=1 (status: $status). Trying scale=300..."
        matpower_data_last_attempt = deepcopy(matpower_data)
        matpower_data_last_attempt["selected_scale"] = 300
        model = build_and_optimize_model(matpower_data_last_attempt, model_type; optimizer)
        status = JuMP.termination_status(model)
    end

    # Fifth and final attempt with scale = 500
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Model still failed with scale=300 (status: $status). Trying scale=500 (final attempt)..."
        matpower_data_extreme_attempt = deepcopy(matpower_data)
        matpower_data_extreme_attempt["selected_scale"] = 500
        model = build_and_optimize_model(matpower_data_extreme_attempt, model_type; optimizer)
    end

    return model
end

function build_and_optimize_singleperiod_model(
    matpower_data,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
)
    model = JuMP.Model(optimizer)
    model.ext[:time_periods] = 1
    model.ext[:times] = [1]
    model.ext[:selected_scale] = get(matpower_data, "selected_scale", nothing)

    process_subproblem(
        model,
        matpower_data,
        model_type,
        SubproblemBuild,
    )

    JuMP.optimize!(model)
    termination_status = JuMP.termination_status(model)
    if termination_status == MOI.OPTIMAL || termination_status == MOI.LOCALLY_SOLVED
        @show JuMP.objective_value(model)
    else
        @warn "Modelo não resolveu com sucesso. Status: $termination_status"
    end
    return model
end

function build_and_optimize_multiperiod_model(
    matpower_data,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
    times::Vector{Any},
)
    model = JuMP.Model(optimizer)
    model.ext[:time_periods] = length(times)
    model.ext[:times] = times

    process_subproblem(model, matpower_data, model_type, SubproblemBuild)

    JuMP.optimize!(model)
    @show JuMP.objective_value(model)
    return model
end

function process_subproblem(args...)
    buses!(args...)
    branches!(args...)
    generators!(args...)
    batteries!(args...)
    prosumer(args...)
    active_power_balance!(args...)
    reactive_power_balance!(args...)
    objective!(args...)

    return nothing
end