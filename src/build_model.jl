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

function build_and_optimize_model(matpower_data, model_type::Type{<:OPFModel})
    has_time_series =
        haskey(matpower_data, "bus_timeseries")

    if has_time_series
        return build_and_optimize_multiperiod_model(matpower_data, model_type)
    else
        return build_and_optimize_singleperiod_model(matpower_data, model_type)
    end
end

function build_and_optimize_singleperiod_model(matpower_data, model_type::Type{<:OPFModel})
    model = JuMP.Model(Gurobi.Optimizer)
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

function build_and_optimize_multiperiod_model(matpower_data, model_type::Type{<:OPFModel})
    times = unique(
        vcat(
            haskey(matpower_data, "bus_timeseries") ?
            matpower_data["bus_timeseries"][!, :time] : [],
            haskey(matpower_data, "bus_timeseries.csv") ?
            matpower_data["bus_timeseries.csv"][!, :time] : [],
        ),
    )
    sort!(times)

    model = JuMP.Model(Gurobi.Optimizer)
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
