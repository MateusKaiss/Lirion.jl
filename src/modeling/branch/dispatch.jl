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

function add_branch_variables!(model::JuMP.Model, sys::SystemData, ::Type{IVR})
    return add_ivr_branch_variables!(model, sys.branch)
end

function add_branch_variables!(model::JuMP.Model, sys::SystemData, ::Type{LPAC})
    return add_lpac_branch_variables!(model, sys.branch)
end

function add_branch_variables!(model::JuMP.Model, sys::SystemData, ::Type{SOCWR})
    return add_rt_variables_with_bounds!(model, sys.bus, sys.branch)
end

function add_branch_variables!(model::JuMP.Model, sys::SystemData, ::Type{SDPWRM})
    # SDP needs the bus matrix variables first (they own v_magnitude_sq).
    add_sdpwrm_bus_variables!(model, sys.bus, sys.branch)
    return add_sdpwrm_branch_variables!(model, sys.bus, sys.branch)
end

function add_branch_variables!(model::JuMP.Model, sys::SystemData, ::Type{SOCBF})
    return add_socbf_branch_variables!(model, sys.bus, sys.branch)
end

function add_branch_variables!(model::JuMP.Model, sys::SystemData, ::Type{BFA})
    return add_bfa_branch_variables!(model, sys.bus, sys.branch)
end

function add_branch_variables!(model::JuMP.Model, sys::SystemData, ::Type{NFA})
    return add_nfa_branch_variables!(model, sys.branch)
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

function add_branch_expressions!(model::JuMP.Model, sys::SystemData, ::Type{IVR})
    return nothing
end

function add_branch_expressions!(model::JuMP.Model, sys::SystemData, ::Type{LPAC})
    return add_lpac_branch_expressions!(model, sys.branch)
end

function add_branch_expressions!(model::JuMP.Model, sys::SystemData, ::Type{SOCWR})
    return add_socwr_branch_expressions!(model, sys.branch)
end

function add_branch_expressions!(model::JuMP.Model, sys::SystemData, ::Type{SDPWRM})
    return add_socwr_branch_expressions!(model, sys.branch)
end

function add_branch_expressions!(model::JuMP.Model, sys::SystemData, ::Type{SOCBF})
    return nothing  # SOCBF: power flows are variables, not expressions
end

function add_branch_expressions!(model::JuMP.Model, sys::SystemData, ::Type{BFA})
    return nothing  # BFA: power flows are variables, not expressions
end

function add_branch_expressions!(model::JuMP.Model, sys::SystemData, ::Type{NFA})
    return add_nfa_branch_expressions!(model, sys.branch)
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

function add_branch_constraints!(model::JuMP.Model, sys::SystemData, ::Type{IVR})
    return add_ivr_branch_constraints!(model, sys.bus, sys.branch)
end

function add_branch_constraints!(model::JuMP.Model, sys::SystemData, ::Type{LPAC})
    return add_lpac_branch_constraints!(model, sys.branch)
end

function add_branch_constraints!(model::JuMP.Model, sys::SystemData, ::Type{SOCWR})
    return add_socwr_branch_constraints!(model, sys.branch)
end

function add_branch_constraints!(model::JuMP.Model, sys::SystemData, ::Type{SDPWRM})
    return add_sdpwrm_branch_constraints!(model, sys.bus, sys.branch)
end

function add_branch_constraints!(model::JuMP.Model, sys::SystemData, ::Type{SOCBF})
    return add_socbf_branch_constraints!(model, sys.bus, sys.branch)
end

function add_branch_constraints!(model::JuMP.Model, sys::SystemData, ::Type{BFA})
    return add_bfa_branch_constraints!(model, sys.bus, sys.branch)
end

function add_branch_constraints!(model::JuMP.Model, sys::SystemData, ::Type{NFA})
    return nothing
end
