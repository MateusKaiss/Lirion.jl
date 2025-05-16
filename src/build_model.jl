abstract type OPFModel end

struct DC <: OPFModel end
struct AC <: OPFModel end
struct SOCP <: OPFModel end

abstract type SubproblemAction end

struct SubproblemBuild <: SubproblemAction end
struct SubproblemUpdate <: SubproblemAction end
struct SubproblemGetResult <: SubproblemAction end
struct SubproblemNameElementsDimensions <: SubproblemAction end
struct SubproblemWriteResult <: SubproblemAction end

function build_and_optimize_model(
    matpower_data,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
)
    has_time_series =
        haskey(matpower_data, "bus_timeseries")

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

    process_subproblem(
        model,
        matpower_data,
        model_type,
        SubproblemBuild,
    )

    JuMP.optimize!(model)
    @show JuMP.objective_value(model)
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
            [
                value(model[:v_magnitude][b, t]) for b in bus_ids, t in 1:T
            ] : nothing,
        "vm_sq" =>
            haskey(model, :v_magnitude_sq) ?
            [
                value(model[:v_magnitude_sq][b, t]) for b in bus_ids, t in 1:T
            ] : nothing,
    )

    if haskey(model, :p_generation)
        gen_ids = 1:length(data["gen"]["gen_bus"])
        result["gen"] = Dict(
            "p" => [
                value(model[:p_generation][g, t]) for g in gen_ids, t in 1:T
            ],
            "q" =>
                haskey(model, :q_generation) ?
                [
                    value(model[:q_generation][g, t]) for g in gen_ids, t in 1:T
                ] : nothing,
            "curtailment" =>
                haskey(model, :curtailment) ?
                [
                    value(model[:curtailment][g, t]) for g in gen_ids, t in 1:T
                ] : nothing,
        )
    end

    if haskey(model, :p_flow_branches)
        branch_ids = 1:length(data["branch"]["f_bus"])
        result["branch"] = Dict(
            "p_from" => [
                value(model[:p_flow_branches][l, t]) for l in branch_ids, t in 1:T
            ],
            "p_to" => [
                value(model[:p_reverse_flow_branches][l, t]) for l in branch_ids, t in 1:T
            ],
            "q_from" =>
                haskey(model, :q_flow_branches) ?
                [
                    value(model[:q_flow_branches][l, t]) for l in branch_ids, t in 1:T
                ] : nothing,
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
            "qcharge" => [value(model[:qcharge][b, t]) for b in 1:nb, t in 1:T],
            "qdischarge" => [value(model[:qdischarge][b, t]) for b in 1:nb, t in 1:T],
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
    # ------------------------------
    # Bus index map: bus ID → index
    # ------------------------------
    if haskey(case, "bus")
        bus_ids = case["bus"]["bus_i"]
        case["bus_idx_map"] = Dict(b => i for (i, b) in enumerate(bus_ids))
    end

    # --------------------------------------------------
    # Generator bus map: generator index → bus index/id
    # --------------------------------------------------
    if haskey(case, "gen")
        gen_bus = case["gen"]["gen_bus"]
        case["gen_idx_map"] = Dict(i => b for (i, b) in enumerate(gen_bus))

        # Inverso: bus → list of generator indices
        bus_to_gens = Dict{Int, Vector{Int}}()
        for (i, b) in enumerate(gen_bus)
            push!(get!(bus_to_gens, b, Int[]), i)
        end
        case["bus_to_gens"] = bus_to_gens
    end

    # --------------------------------------
    # Branch maps: bus → branches (from/to)
    # --------------------------------------
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

    # ---------------------------------------------------
    # Prosumers and storage maps (if available)
    # ---------------------------------------------------
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

    return case
end
