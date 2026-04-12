function generators!(
    model::JuMP.Model,
    sys::SystemData,
    model_type::Type{<:OPFModel},
)
    gen = sys.gen

    if gen.n == 0
        @info "No generators in SystemData. Skipping generator modeling."
        return nothing
    end

    add_generator_variables!(model, gen, model_type)
    add_generator_constraints!(model, gen, model_type)

    return nothing
end

# =======================================================================
# VARIABLES
# =======================================================================

function add_generator_variables!(model::JuMP.Model, gen::GenData, ::Type{DC})
    add_active_generation_variables!(model, gen)
    return add_curtailment_variables!(model, gen)
end

function add_generator_variables!(model::JuMP.Model, gen::GenData, ::Type{AC})
    add_active_generation_variables!(model, gen)
    add_reactive_generation_variables!(model, gen)
    return add_curtailment_variables!(model, gen)
end

function add_generator_variables!(model::JuMP.Model, gen::GenData, ::Type{ACR})
    add_active_generation_variables!(model, gen)
    add_reactive_generation_variables!(model, gen)
    return add_curtailment_variables!(model, gen)
end

function add_generator_variables!(model::JuMP.Model, gen::GenData, ::Type{SOCWR})
    add_active_generation_variables!(model, gen)
    add_reactive_generation_variables!(model, gen)
    return add_curtailment_variables!(model, gen)
end

# =======================================================================
# CONSTRAINTS
# =======================================================================

function add_generator_constraints!(model::JuMP.Model, gen::GenData, ::Type{DC})
    return add_curtailment_constraints!(model, gen)
end

function add_generator_constraints!(model::JuMP.Model, gen::GenData, ::Type{AC})
    return add_curtailment_constraints!(model, gen)
end

function add_generator_constraints!(model::JuMP.Model, gen::GenData, ::Type{ACR})
    return add_curtailment_constraints!(model, gen)
end

function add_generator_constraints!(model::JuMP.Model, gen::GenData, ::Type{SOCWR})
    return add_curtailment_constraints!(model, gen)
end
