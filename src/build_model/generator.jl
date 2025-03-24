function get_generator_parameters(gen_data)
    return (
        number_of_generators = length(gen_data["gen_bus"]),
        existing = gen_data["gen_status"],
        pmn = gen_data["pmin"],
        pmx = gen_data["pmax"],
        qmn = gen_data["qmin"],
        qmx = gen_data["qmax"],
        initial_p_generation = gen_data["pg"],
        initial_q_generation = gen_data["qg"],
        plant_type = gen_data["plant_type"],
        curtailment_cost = gen_data["curtailment_cost"],
    )
end

function has_curtailment(gen_data, gen_idx)
    return gen_data["gen_status"][gen_idx] == 1 && gen_data["plant_type"][gen_idx] == "R"
end

function create_active_generation_variables!(model, matpower_data)
    gen_data = matpower_data["gen"]
    params = get_generator_parameters(gen_data)

    @variable(
        model,
        p_generation[
            g = 1:params.number_of_generators;
            params.existing[g] == 1
        ],
        lower_bound = params.pmn[g],
        upper_bound = params.pmx[g],
        start =
            isnan(params.initial_p_generation[g]) ? nothing :
            params.initial_p_generation[g]
    )
end

function create_reactive_generation_variables!(model, matpower_data)
    gen_data = matpower_data["gen"]
    params = get_generator_parameters(gen_data)

    @variable(
        model,
        q_generation[
            g = 1:params.number_of_generators;
            params.existing[g] == 1
        ],
        lower_bound = params.qmn[g],
        upper_bound = params.qmx[g],
        start =
            isnan(params.initial_q_generation[g]) ? nothing :
            params.initial_q_generation[g]
    )
end

function create_curtailment_variables!(model, matpower_data)
    gen_data = matpower_data["gen"]
    params = get_generator_parameters(gen_data)

    @variable(
        model,
        curtailment[g = 1:params.number_of_generators],
        lower_bound = 0
    )
end

function create_curtailment_constraints!(model, matpower_data)
    gen_data = matpower_data["gen"]
    params = get_generator_parameters(gen_data)

    for g in 1:params.number_of_generators
        if has_curtailment(gen_data, g)
            @constraint(
                model,
                model[:p_generation][g] + model[:curtailment][g] == params.pmx[g]
            )
        end
    end
end

# ---------------------------------------------------------------------------
# Dispatching Generator Variables, Expressions and Constraints by Model Type
# ---------------------------------------------------------------------------
function generators!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    create_active_generation_variables!(model, matpower_data)
    create_curtailment_variables!(model, matpower_data)
    create_curtailment_constraints!(model, matpower_data)
    return nothing
end

function generators!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    create_active_generation_variables!(model, matpower_data)
    create_reactive_generation_variables!(model, matpower_data)
    create_curtailment_variables!(model, matpower_data)
    create_curtailment_constraints!(model, matpower_data)
    return nothing
end

function generators!(model::JuMP.Model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    create_active_generation_variables!(model, matpower_data)
    create_reactive_generation_variables!(model, matpower_data)
    create_curtailment_variables!(model, matpower_data)
    create_curtailment_constraints!(model, matpower_data)
    return nothing
end
