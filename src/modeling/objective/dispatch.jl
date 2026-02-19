
# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================
function objective!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
)
    obj = sys.objective
    compute_generation_cost(model, obj)

    @objective(model, Min,
        model[:generation_cost]
    )
end

function objective!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
)
    obj = sys.objective
    compute_generation_cost(model, obj)

    @objective(model, Min,
        model[:generation_cost]
    )
end

function objective!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCP},
)
    obj = sys.objective
    compute_generation_cost(model, obj)

    @objective(model, Min,
        model[:generation_cost]
    )
end
