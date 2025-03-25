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

function build_and_optimize_model(matpower_data, model_type)
    # problem_results = ProblemResults()
    model = JuMP.Model(Gurobi.Optimizer)

    process_subproblem(
        model,
        matpower_data,
        model_type,
        SubproblemBuild,
    )

    JuMP.optimize!(model)

    @show JuMP.objective_value(model)

    return model

    return nothing
end

function process_subproblem(args...)
    buses!(args...)
    branches!(args...)
    generators!(args...)
    batteries!(args...)
    active_power_balance!(args...)
    reactive_power_balance!(args...)
    objective!(args...)

    return nothing
end
