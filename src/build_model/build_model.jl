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

function build_and_optimize_model(
    matpower_data,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
)
    optimizer === nothing && (optimizer = Ipopt.Optimizer)
    has_time_series = haskey(matpower_data, "bus_timeseries")
    if has_time_series
        return build_and_optimize_multiperiod_model(matpower_data, model_type; optimizer)
    else
        return build_and_optimize_singleperiod_model(matpower_data, model_type; optimizer)
    end
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
)
    times = unique(
        vcat(
            haskey(matpower_data, "bus_timeseries") ?
            matpower_data["bus_timeseries"][!, :time] : [],
            haskey(matpower_data, "bus_timeseries.csv") ?
            matpower_data["bus_timeseries.csv"][!, :time] : [],
        ),
    )
    sort!(times)

    model = JuMP.Model(optimizer)
    model.ext[:time_periods] = length(times)
    model.ext[:times] = times

    process_subproblem(model, matpower_data, model_type, SubproblemBuild)

    JuMP.optimize!(model)
    @show JuMP.objective_value(model)
    return model
end

function try_with_fallback_scale(
    matpower_data::Dict,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
)
    scale = select_cone_scale(matpower_data["bus_params"])

    # Primeira tentativa
    matpower_data_attempt = deepcopy(matpower_data)
    matpower_data_attempt["selected_scale"] = scale
    model = build_and_optimize_model(matpower_data_attempt, model_type; optimizer)
    status = JuMP.termination_status(model)

    # Segunda tentativa com scale = 150
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Modelo com scale=$scale falhou (status: $status). Tentando scale=150..."
        matpower_data_fallback = deepcopy(matpower_data)
        matpower_data_fallback["selected_scale"] = 150
        model = build_and_optimize_model(matpower_data_fallback, model_type; optimizer)
        status = JuMP.termination_status(model)
    end

    # Terceira tentativa com scale = 1
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Modelo ainda falhou com scale=150 (status: $status). Tentando scale=1..."
        matpower_data_final_attempt = deepcopy(matpower_data)
        matpower_data_final_attempt["selected_scale"] = 1
        model = build_and_optimize_model(matpower_data_final_attempt, model_type; optimizer)
        status = JuMP.termination_status(model)
    end

    # Quarta tentativa com scale = 300
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Modelo ainda falhou com scale=1 (status: $status). Tentando scale=300..."
        matpower_data_last_attempt = deepcopy(matpower_data)
        matpower_data_last_attempt["selected_scale"] = 300
        model = build_and_optimize_model(matpower_data_last_attempt, model_type; optimizer)
        status = JuMP.termination_status(model)
    end

    # Quinta tentativa com scale = 500
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Modelo ainda falhou com scale=300 (status: $status). Tentando scale=500 (última tentativa)..."
        matpower_data_extreme_attempt = deepcopy(matpower_data)
        matpower_data_extreme_attempt["selected_scale"] = 500
        model =
            build_and_optimize_model(matpower_data_extreme_attempt, model_type; optimizer)
    end

    return model
end
