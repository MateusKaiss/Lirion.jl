# ---------------------------------------------------------------------------
# Generator Parameters (precomputed once and passed to all helper functions)
# ---------------------------------------------------------------------------
function get_generator_parameters(gen_data)
    number_of_generators = length(gen_data["gen_bus"])
    start_pg = Dict{Int, Union{Nothing, Float64}}()
    start_qg = Dict{Int, Union{Nothing, Float64}}()
    for g in 1:number_of_generators
        start_pg[g] = isnan(gen_data["pg"][g]) ? nothing : gen_data["pg"][g]
        start_qg[g] = isnan(gen_data["qg"][g]) ? nothing : gen_data["qg"][g]
    end

    return (
        number_of_generators = number_of_generators,
        existing = gen_data["gen_status"] .== 1,
        pmn = gen_data["pmin"],
        pmx = gen_data["pmax"],
        qmn = gen_data["qmin"],
        qmx = gen_data["qmax"],
        initial_pg = start_pg,
        initial_qg = start_qg,
        plant_type = gen_data["plant_type"],
        curtailment_cost = gen_data["curtailment_cost"],
    )
end

function has_curtailment(params, g)
    return params.existing[g] && params.plant_type[g] == "R"
end

function create_active_generation_variables!(model, params)
    T = get(model.ext, :time_periods, 1)

    @variable(
        model,
        p_generation[g = 1:params.number_of_generators, t = 1:T; params.existing[g]],
        lower_bound = params.pmn[g],
        upper_bound = params.pmx[g],
        start = params.initial_pg[g]
    )
end

function create_reactive_generation_variables!(model, params)
    T = get(model.ext, :time_periods, 1)

    @variable(
        model,
        q_generation[g = 1:params.number_of_generators, t = 1:T; params.existing[g]],
        lower_bound = params.qmn[g],
        upper_bound = params.qmx[g],
        start = params.initial_qg[g]
    )
end

function create_curtailment_variables!(model, params)
    T = get(model.ext, :time_periods, 1)
    @variable(model, curtailment[g = 1:params.number_of_generators, t = 1:T] >= 0)
end

function create_curtailment_constraints!(model, params)
    T = get(model.ext, :time_periods, 1)
    curtailable = filter(g -> has_curtailment(params, g), 1:params.number_of_generators)

    @constraint(model, [g in curtailable, t in 1:T],
        model[:p_generation][g, t] + model[:curtailment][g, t] == params.pmx[g]
    )
end

# ---------------------------------------------------------------------------
# Dispatching Generator Variables, Expressions and Constraints by Model Type
# ---------------------------------------------------------------------------
function generators!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    params = matpower_data["generator_params"]
    create_active_generation_variables!(model, params)
    create_curtailment_variables!(model, params)
    create_curtailment_constraints!(model, params)
    return nothing
end

function generators!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = matpower_data["generator_params"]
    create_active_generation_variables!(model, params)
    create_reactive_generation_variables!(model, params)
    create_curtailment_variables!(model, params)
    create_curtailment_constraints!(model, params)
    return nothing
end

function generators!(
    model::JuMP.Model,
    matpower_data,
    ::Type{SOCP},
    ::Type{SubproblemBuild},
)
    params = matpower_data["generator_params"]
    create_active_generation_variables!(model, params)
    create_reactive_generation_variables!(model, params)
    create_curtailment_variables!(model, params)
    create_curtailment_constraints!(model, params)
    return nothing
end
