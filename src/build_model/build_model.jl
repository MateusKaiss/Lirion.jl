# =========================
# BUILD ENTRYPOINT
# =========================
@inline _select_default_optimizer(::Type{DC}) = OSQP.Optimizer
@inline _select_default_optimizer(::Type{SOCP}) = SCS.Optimizer
@inline _select_default_optimizer(::Type{<:OPFModel}) = Ipopt.Optimizer

function build_model(
    sys::SystemData,
    model_type::Type{<:OPFModel};
    optimizer = nothing,
)
    opt = optimizer === nothing ? _select_default_optimizer(model_type) : optimizer
    times = get_model_times(sys)

    if length(times) == 1
        return _build_singleperiod_model(sys, model_type, opt)
    else
        return _build_multiperiod_model(sys, model_type, opt, times)
    end
end

function build_model(inputs::LirionInputs; optimizer = nothing)
    sys = inputs.system
    model_type = inputs.model_type::Type{<:OPFModel}

    opt =
        optimizer !== nothing ? optimizer :
        inputs.optimizer !== nothing ? inputs.optimizer :
        _select_default_optimizer(model_type)

    return build_model(sys, model_type; optimizer = opt)
end

# =========================
# CENTRALIZED (DC/AC/SOCP)
# =========================
function get_model_times(sys::SystemData)
    time_vecs = Vector{Vector{Int}}()

    if !isempty(sys.bus_ts.time)
        push!(time_vecs, sys.bus_ts.time)
    end
    if !isempty(sys.prosumer_ts.time)
        push!(time_vecs, sys.prosumer_ts.time)
    end
    if !isempty(sys.consumer_ts.time)
        push!(time_vecs, sys.consumer_ts.time)
    end
    if !isempty(sys.gencost_ts.time)
        push!(time_vecs, sys.gencost_ts.time)
    end

    if isempty(time_vecs)
        return collect(1:sys.T)
    end

    times = unique(vcat(time_vecs...))
    sort!(times)
    return times
end

# ============================================
# SINGLE-PERIOD (build only)
# ============================================
function _build_singleperiod_model(
    sys::SystemData,
    model_type::Type{<:OPFModel},
    optimizer,
)
    model = JuMP.Model(optimizer)

    model.ext[:time_periods] = 1
    model.ext[:times] = [1]

    process_subproblem(model, sys, model_type, SubproblemBuild)

    return model
end

# ============================================
# MULTI-PERIOD (build only)
# ============================================

function _build_multiperiod_model(
    sys::SystemData,
    model_type::Type{<:OPFModel},
    optimizer,
    times::AbstractVector,
)
    model = JuMP.Model(optimizer)

    model.ext[:time_periods] = length(times)
    model.ext[:times] = collect(times)

    process_subproblem(model, sys, model_type, SubproblemBuild)

    return model
end

# ============================================
# OPTIMIZATION HELPER
# ============================================
function optimize_model!(model::JuMP.Model; label::AbstractString = "centralized OPF")
    JuMP.optimize!(model)
    status = JuMP.termination_status(model)

    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
        @info "Objective value ($label) = $(JuMP.objective_value(model)) | status=$status"
    else
        @warn "Model did not solve successfully ($label). Status: $status"
    end

    return model
end

function solve_model(args...; kwargs...)
    model = build_model(args...; kwargs...)
    optimize_model!(model)
    return model
end

# ============================================
# SUBPROBLEM PIPELINE
# ============================================
@inline function process_subproblem(
    model::JuMP.Model,
    sys::SystemData,
    model_type::Type{<:OPFModel},
    action_type::Type{<:SubproblemAction},
)
    buses!(model, sys, model_type, action_type)
    branches!(model, sys, model_type, action_type)
    generators!(model, sys, model_type, action_type)
    batteries!(model, sys, model_type, action_type)
    prosumer!(model, sys, model_type, action_type)
    active_power_balance!(model, sys, model_type, action_type)
    reactive_power_balance!(model, sys, model_type, action_type)
    objective!(model, sys, model_type, action_type)
    return nothing
end
