function branches!(
    model::JuMP.Model,
    sys::SystemData,
    model_type::Type{<:OPFModel},
)
    add_branch_variables!(model, sys, model_type)
    add_branch_expressions!(model, sys, model_type)
    add_branch_constraints!(model, sys, model_type)
    return nothing
end

# =======================================================================
# VARIABLES
# =======================================================================

function add_branch_variables!(model::JuMP.Model, sys::SystemData, ::Type{DC})
    return nothing
end

function add_branch_variables!(model::JuMP.Model, sys::SystemData, ::Type{AC})
    return nothing
end

function add_branch_variables!(model::JuMP.Model, sys::SystemData, ::Type{ACR})
    return nothing
end

function add_branch_variables!(model::JuMP.Model, sys::SystemData, ::Type{SOCWR})
    return add_rt_variables_with_bounds!(model, sys.bus, sys.branch)
end

# =======================================================================
# EXPRESSIONS
# =======================================================================

function add_branch_expressions!(model::JuMP.Model, sys::SystemData, ::Type{DC})
    return add_dc_branch_expressions!(model, sys.branch)
end

function add_branch_expressions!(model::JuMP.Model, sys::SystemData, ::Type{AC})
    return add_ac_branch_expressions!(model, sys.branch)
end

function add_branch_expressions!(model::JuMP.Model, sys::SystemData, ::Type{ACR})
    return add_acr_branch_expressions!(model, sys.branch)
end

function add_branch_expressions!(model::JuMP.Model, sys::SystemData, ::Type{SOCWR})
    return add_socwr_branch_expressions!(model, sys.branch)
end

# =======================================================================
# CONSTRAINTS
# =======================================================================

function add_branch_constraints!(model::JuMP.Model, sys::SystemData, ::Type{DC})
    return add_dc_branch_constraints!(model, sys.branch)
end

function add_branch_constraints!(model::JuMP.Model, sys::SystemData, ::Type{AC})
    return add_ac_branch_constraints!(model, sys.branch)
end

function add_branch_constraints!(model::JuMP.Model, sys::SystemData, ::Type{ACR})
    return add_acr_branch_constraints!(model, sys.branch)
end

function add_branch_constraints!(model::JuMP.Model, sys::SystemData, ::Type{SOCWR})
    return add_socwr_branch_constraints!(model, sys.branch)
end
