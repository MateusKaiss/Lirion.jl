abstract type OPFModel end

struct DC <: OPFModel end
struct AC <: OPFModel end
struct SOCP <: OPFModel end

function build_and_optimize_model(matpower_data, model_type)
    # problem_results = ProblemResults()
    model = JuMP.Model(Ipopt.Optimizer)

    process_subproblem(
        model,
        matpower_data,
        model_type,
        SubproblemBuild,
    )

    JuMP.optimize!(model)

    return nothing
end

function process_subproblem(args...)
    # buses!(args...)
    # branches!(args...)
    # generators!(args...)
    # batteries!(args...)
    # active_power_balance!(args...)
    # reactive_power_balance!(args...)
    # objective!(args...)

    return nothing
end