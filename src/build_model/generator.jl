function get_generator_parameters(inputs)
    return (
        number_of_generators = length(inputs["gen_bus"]),
        existing = inputs["gen_status"],
        pmn = inputs["pmin"],
        pmx = inputs["pmax"],
        qmn = inputs["qmin"],
        qmx = inputs["qmax"],
        initial_p_generation = inputs["pg"],
        initial_q_generation = inputs["qg"],
        plant_type = inputs["plant_type"],
        curtailment_cost = inputs["curtailment_cost"],
    )
end

function has_curtailment(inputs, gen_idx)
    return inputs["gen_status"][gen_idx] == 1 && inputs["plant_type"][gen_idx] == "R"
end

function create_active_generation_variables!(model, inputs)
    params = get_generator_parameters(inputs)

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

function create_reactive_generation_variables!(model, inputs)
    params = get_generator_parameters(inputs)

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

function create_curtailment_variables!(model, inputs)
    params = get_generator_parameters(inputs)

    @variable(
        model,
        curtailment[g = 1:params.number_of_generators],
        lower_bound = 0
    )
end

function create_curtailment_constraints!(model, inputs)
    params = get_generator_parameters(inputs)

    for g in 1:params.number_of_generators
        if has_curtailment(inputs, g)
            @constraint(
                model,
                model[:p_generation][g] + model[:curtailment][g] == params.pmx[g]
            )
        end
    end
end

# ---------------------------------------------------------------------
# Dispatching Generator Variables by Model Type
# ---------------------------------------------------------------------
function generators!(model::JuMP.Model, inputs, ::Type{DC}, ::Type{SubproblemBuild})
    create_active_generation_variables!(model, inputs)
    create_curtailment_variables!(model, inputs)
    create_curtailment_constraints!(model, inputs)
    return nothing
end

function generators!(model::JuMP.Model, inputs, ::Type{AC}, ::Type{SubproblemBuild})
    create_active_generation_variables!(model, inputs)
    create_reactive_generation_variables!(model, inputs)
    create_curtailment_variables!(model, inputs)
    create_curtailment_constraints!(model, inputs)
    return nothing
end

function generators!(model::JuMP.Model, inputs, ::Type{SOCP}, ::Type{SubproblemBuild})
    create_active_generation_variables!(model, inputs)
    create_reactive_generation_variables!(model, inputs)
    create_curtailment_variables!(model, inputs)
    create_curtailment_constraints!(model, inputs)
    return nothing
end
