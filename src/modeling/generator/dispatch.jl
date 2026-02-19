
# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================

function generators!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
)
    gen = sys.gen
    if gen.n == 0
        @info "No generators in SystemData (DC). Skipping generator modeling."
        return nothing
    end

    create_active_generation_variables!(model, gen)
    create_curtailment_variables!(model, gen)
    create_curtailment_constraints!(model, gen)
    return nothing
end

function generators!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
)
    gen = sys.gen
    if gen.n == 0
        @info "No generators in SystemData (AC). Skipping generator modeling."
        return nothing
    end

    create_active_generation_variables!(model, gen)
    create_reactive_generation_variables!(model, gen)
    create_curtailment_variables!(model, gen)
    create_curtailment_constraints!(model, gen)
    return nothing
end

function generators!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCP},
)
    gen = sys.gen
    if gen.n == 0
        @info "No generators in SystemData (SOCP). Skipping generator modeling."
        return nothing
    end

    create_active_generation_variables!(model, gen)
    create_reactive_generation_variables!(model, gen)
    create_curtailment_variables!(model, gen)
    create_curtailment_constraints!(model, gen)
    return nothing
end
