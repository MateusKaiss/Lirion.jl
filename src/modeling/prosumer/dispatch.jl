function prosumers!(
    model::JuMP.Model,
    sys::SystemData,
    model_type::Type{<:OPFModel},
)
    pros = sys.prosumer

    if pros.number_of_prosumers == 0
        @info "No prosumers found. Skipping prosumer modeling."
        return nothing
    end

    add_prosumer_variables!(model, sys, model_type)
    add_prosumer_constraints!(model, sys, model_type)

    return nothing
end

# =======================================================================
# VARIABLES
# =======================================================================

function add_prosumer_variables!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{<:OPFModel},
)
    return add_prosumer_variables!(model, sys.prosumer)
end

# =======================================================================
# CONSTRAINTS
# =======================================================================

function add_prosumer_constraints!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{<:OPFModel},
)
    return add_prosumer_constraints!(model, sys.prosumer)
end
