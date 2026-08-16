# ============================================================
# OPTIMIZATION / SOLVE HELPERS
# ============================================================

function optimize_model!(
    model::JuMP.Model;
    label::AbstractString = "centralized OPF",
    allow_local::Bool = true,
)
    JuMP.optimize!(model)
    status = JuMP.termination_status(model)

    ok = status == MOI.OPTIMAL || (allow_local && status == MOI.LOCALLY_SOLVED)

    if ok
        @info "Objective value ($label) = $(JuMP.objective_value(model)) | status=$status"
    else
        @warn "Model did not solve successfully ($label). Status: $status"
    end

    return model
end

function solve_model(
    inputs::CentralizedInputs;
    label::AbstractString = "centralized OPF",
    allow_local::Bool = true,
)
    model = build_model(inputs)
    optimize_model!(model; label = label, allow_local = allow_local)
    return model
end
