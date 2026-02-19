
# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================
function prosumer!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{T},
) where {T <: OPFModel}
    pros = sys.prosumer

    if pros.number_of_prosumers == 0
        @info "No prosumers found. Skipping prosumer modeling."
        return nothing
    end

    create_prosumer_variables!(model, pros)
    create_prosumer_constraints!(model, pros)

    return nothing
end
