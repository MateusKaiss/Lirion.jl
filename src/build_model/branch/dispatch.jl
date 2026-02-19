# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================

function branches!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
)
    br = sys.branch
    create_dc_branch_expressions!(model, br)
    create_dc_branch_constraints!(model, br)
    return nothing
end

function branches!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
)
    br = sys.branch
    create_ac_branch_expressions!(model, br)
    create_ac_branch_constraints!(model, br)
    return nothing
end

function branches!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCP},
)
    bus = sys.bus
    br = sys.branch

    create_rt_variables_with_bounds!(model, bus, br)
    create_socp_branch_expressions!(model, br)
    create_socp_branch_constraints!(model, br)

    return nothing
end
