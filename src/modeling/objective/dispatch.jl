function objective!(
    model::JuMP.Model,
    sys::SystemData,
    model_type::Type{<:OPFModel},
)
    obj = sys.objective

    add_objective_expressions!(model, obj, model_type)
    add_objective_function!(model)

    return nothing
end

# =======================================================================
# EXPRESSIONS
# =======================================================================

function add_objective_expressions!(model::JuMP.Model, obj::ObjectiveData, ::Type{DC})
    return add_generation_cost!(model, obj)
end

function add_objective_expressions!(model::JuMP.Model, obj::ObjectiveData, ::Type{AC})
    return add_generation_cost!(model, obj)
end

function add_objective_expressions!(model::JuMP.Model, obj::ObjectiveData, ::Type{ACR})
    return add_generation_cost!(model, obj)
end

function add_objective_expressions!(model::JuMP.Model, obj::ObjectiveData, ::Type{IVR})
    return add_generation_cost!(model, obj)
end

function add_objective_expressions!(model::JuMP.Model, obj::ObjectiveData, ::Type{LPAC})
    return add_generation_cost!(model, obj)
end

function add_objective_expressions!(model::JuMP.Model, obj::ObjectiveData, ::Type{SOCWR})
    return add_generation_cost!(model, obj)
end

# =======================================================================
# OBJECTIVE FUNCTION
# =======================================================================

function add_objective_function!(model::JuMP.Model)
    @objective(model, Min,
        model[:generation_cost]
    )
    return nothing
end
