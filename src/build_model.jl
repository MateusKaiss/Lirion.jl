function build_and_optimize_model(
    matpower_data,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
)
    convert_to_pu!(matpower_data)
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

function build_result_dict(model::JuMP.Model, data::Dict)
    T = get(model.ext, :time_periods, 1)
    result = Dict{String, Any}()

    result["objective"] = objective_value(model)
    bus_ids = data["bus"]["bus_i"]
    result["bus"] = Dict(
        "theta" => [
            value(model[:theta][b, t]) for b in bus_ids, t in 1:T
        ],
        "vm" =>
            haskey(model, :v_magnitude) ?
            [value(model[:v_magnitude][b, t]) for b in bus_ids, t in 1:T] : nothing,
        "vm_sq" =>
            haskey(model, :v_magnitude_sq) ?
            [value(model[:v_magnitude_sq][b, t]) for b in bus_ids, t in 1:T] : nothing,
    )

    if haskey(model, :p_generation)
        gen_status = data["gen"]["gen_status"]
        gen_ids = [g for g in 1:length(gen_status) if gen_status[g] == 1]

        result["gen"] = Dict(
            "p" => [value(model[:p_generation][g, t]) for g in gen_ids, t in 1:T],
            "q" =>
                haskey(model, :q_generation) ?
                [value(model[:q_generation][g, t]) for g in gen_ids, t in 1:T] :
                nothing,
            "curtailment" =>
                haskey(model, :curtailment) ?
                [
                    try
                        value(model[:curtailment][g, t])
                    catch
                        0.0
                    end
                    for g in gen_ids, t in 1:T
                ] : nothing)
    end

    if haskey(model, :p_flow_branches)
        branch_status = data["branch"]["br_status"]
        branch_ids = [g for g in 1:length(branch_status) if branch_status[g] == 1]
        result["branch"] = Dict(
            "p_from" =>
                [value(model[:p_flow_branches][l, t]) for l in branch_ids, t in 1:T],
            "p_to" => [
                value(model[:p_reverse_flow_branches][l, t]) for l in branch_ids, t in 1:T
            ],
            "q_from" =>
                haskey(model, :q_flow_branches) ?
                [value(model[:q_flow_branches][l, t]) for l in branch_ids, t in 1:T] :
                nothing,
            "q_to" =>
                haskey(model, :q_reverse_flow_branches) ?
                [
                    value(model[:q_reverse_flow_branches][l, t]) for l in branch_ids,
                    t in 1:T
                ] : nothing,
        )
    end

    if haskey(model, :storage)
        nb = length(data["storage"]["storage_bus"])
        result["storage"] = Dict(
            "storage" => [value(model[:storage][b, t]) for b in 1:nb, t in 1:T],
            "charge" => [value(model[:charge][b, t]) for b in 1:nb, t in 1:T],
            "discharge" => [value(model[:discharge][b, t]) for b in 1:nb, t in 1:T],
            "qsc" => [value(model[:qsc][b, t]) for b in 1:nb, t in 1:T],
            "ccms" =>
                haskey(model, :ccms) ?
                [value(model[:ccms][b, t]) for b in 1:nb, t in 1:T] : nothing,
            "p_injected" =>
                [value(model[:p_battery_injected][b, t]) for b in 1:nb, t in 1:T],
            "q_injected" =>
                [value(model[:q_battery_injected][b, t]) for b in 1:nb, t in 1:T],
        )
    end

    if haskey(model, :p_prosumer_injected)
        np = length(data["prosumer"]["bus"])
        result["prosumer"] = Dict(
            "p_injected" =>
                [value(model[:p_prosumer_injected][p, t]) for p in 1:np, t in 1:T],
            "q_injected" =>
                [value(model[:q_prosumer_injected][p, t]) for p in 1:np, t in 1:T],
            "pv_used" =>
                [value(model[:pv_generation_used][p, t]) for p in 1:np, t in 1:T],
            "curtailment" =>
                [value(model[:pv_curtailment][p, t]) for p in 1:np, t in 1:T],
        )
    end

    return result
end

function add_case_index_maps!(case::Dict{String, Any})
    if haskey(case, "bus")
        bus_ids = case["bus"]["bus_i"]
        case["bus_idx_map"] = Dict(b => i for (i, b) in enumerate(bus_ids))
    end

    if haskey(case, "gen")
        gen_bus = case["gen"]["gen_bus"]
        case["gen_idx_map"] = Dict(i => b for (i, b) in enumerate(gen_bus))

        bus_to_gens = Dict{Int, Vector{Int}}()
        for (i, b) in enumerate(gen_bus)
            push!(get!(bus_to_gens, b, Int[]), i)
        end
        case["bus_to_gens"] = bus_to_gens
    end

    if haskey(case, "branch")
        f_bus = case["branch"]["f_bus"]
        t_bus = case["branch"]["t_bus"]

        bus_to_branches_from = Dict{Int, Vector{Int}}()
        bus_to_branches_to = Dict{Int, Vector{Int}}()

        for (i, f) in enumerate(f_bus)
            push!(get!(bus_to_branches_from, f, Int[]), i)
        end
        for (i, t) in enumerate(t_bus)
            push!(get!(bus_to_branches_to, t, Int[]), i)
        end

        case["bus_to_branches_from"] = bus_to_branches_from
        case["bus_to_branches_to"] = bus_to_branches_to
    end

    if haskey(case, "prosumer")
        prosumer_bus = case["prosumer"]["bus"]
        bus_to_prosumers = Dict{Int, Vector{Int}}()
        for (i, b) in enumerate(prosumer_bus)
            push!(get!(bus_to_prosumers, b, Int[]), i)
        end
        case["bus_to_prosumers"] = bus_to_prosumers
    end

    if haskey(case, "storage")
        storage_bus = case["storage"]["storage_bus"]
        bus_to_batteries = Dict{Int, Vector{Int}}()
        for (i, b) in enumerate(storage_bus)
            push!(get!(bus_to_batteries, b, Int[]), i)
        end
        case["bus_to_batteries"] = bus_to_batteries
    end

    case["bus_params"] = get_bus_parameters(case["bus"])
    case["branch_params"] = get_branch_parameters(case["branch"], case)
    case["generator_params"] = get_generator_parameters(case["gen"])
    case["power_balance_params"] = get_power_balance_parameters(case)
    case["reactive_balance_params"] = get_reactive_power_parameters(case)
    case["objective_params"] = get_objective_parameters(
        case["gen"],
        get(case, "storage", Dict()),
        case["gencost"],
        get(case, "prosumer", Dict()),
    )

    case["objective_params"] = merge(
        case["objective_params"],
        (; sbase = case["baseMVA"]),
    )

    return case
end

function convert_to_pu!(case::Dict{String, Any})
    Sbase = case["baseMVA"]

    for field in ["pd", "qd", "gs", "bs"]
        case["bus"][field] ./= Sbase
    end

    for field in ["pg", "qg", "pmax", "pmin", "qmax", "qmin"]
        case["gen"][field] ./= Sbase
    end

    for field in ["rate_a", "rate_b", "rate_c"]
        case["branch"][field] ./= Sbase
    end

    # if haskey(case, "storage")
    #     for field in ["ps", "qs", "energy", "energy_rating", "charge_rating", "discharge_rating",
    #                 "thermal_rating", "qmin", "qmax", "q_loss", "p_loss"]
    #         case["storage"][field] ./= Sbase
    #     end
    # end
    if haskey(case, "storage")
        for field in
            ["ps", "qs", "energy", "energy_rating", "charge_rating", "discharge_rating",
            "thermal_rating", "qmin", "qmax", "q_loss", "p_loss"]
            case["storage"][field] ./= Sbase
        end
    end

    return case
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

function parse_case_from_csv_silent(dir::String)
    # Silencia logs dentro desta função
    original_logger = current_logger()
    logger = NullLogger()
    global_logger(logger)
    result = parse_case_from_csv(dir)
    global_logger(original_logger)
    return result
end

function parse_matpower_string_silent(data_string::String)
    original_logger = current_logger()
    global_logger(NullLogger())
    case = _parse_matpower_string(data_string)
    global_logger(original_logger)
    return case
end
