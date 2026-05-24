# ============================================================
# Utils
# ============================================================

@inline _select_default_optimizer(::Type{DC}) = OSQP.Optimizer
@inline _select_default_optimizer(::Type{<:OPFModel}) = Ipopt.Optimizer
@inline _select_default_optimizer(::Type{ACR}) = Ipopt.Optimizer
@inline _select_default_optimizer(::Type{IVR}) = Ipopt.Optimizer
@inline _select_default_optimizer(::Type{LPAC}) = Ipopt.Optimizer

@inline _select_default_optimizer(::Type{SOCWR}) = Clarabel.Optimizer
@inline _select_default_optimizer(::Type{QC}) = Ipopt.Optimizer
@inline _select_default_optimizer(::Type{SDPWRM}) = Clarabel.Optimizer
@inline _select_default_optimizer(::Type{SOCBF}) = Clarabel.Optimizer
@inline _select_default_optimizer(::Type{BFA}) = Clarabel.Optimizer
@inline _select_default_optimizer(::Type{NFA}) = OSQP.Optimizer

@inline function _resolve_optimizer(
    user_opt,
    input_opt,
    model_type::Type{<:OPFModel},
)
    return something(user_opt, input_opt, _select_default_optimizer(model_type))
end

# ============================================================
# TIMES (CENTRALIZED)
# ============================================================

function get_model_times(sys::SystemData)
    time_vecs = Vector{Vector{Int}}()

    !isempty(sys.bus_ts.time) && push!(time_vecs, sys.bus_ts.time)
    !isempty(sys.gencost_ts.time) && push!(time_vecs, sys.gencost_ts.time)

    if isempty(time_vecs)
        return collect(1:sys.T)
    end

    times = unique(vcat(time_vecs...))
    sort!(times)
    return times
end

# ============================================================
# SUBPROBLEM PIPELINE
# ============================================================

@inline function assemble_opf!(
    model::JuMP.Model,
    sys::SystemData,
    model_type::Type{<:OPFModel},
)
    buses!(model, sys, model_type)
    branches!(model, sys, model_type)
    generators!(model, sys, model_type)
    storages!(model, sys, model_type)
    active_power_balance!(model, sys, model_type)
    reactive_power_balance!(model, sys, model_type)
    objective!(model, sys, model_type)
    return nothing
end

# ============================================================
# BUILDER (single / multi)
# ============================================================

function _build_singleperiod_model(sys::SystemData, model_type::Type{<:OPFModel}, optimizer)
    model = JuMP.Model(optimizer)
    model.ext[:time_periods] = 1
    model.ext[:times] = [1]
    assemble_opf!(model, sys, model_type)
    return model
end

function _build_multiperiod_model(
    sys::SystemData,
    model_type::Type{<:OPFModel},
    optimizer,
    times::AbstractVector,
)
    model = JuMP.Model(optimizer)
    model.ext[:time_periods] = length(times)
    model.ext[:times] = collect(times)
    assemble_opf!(model, sys, model_type)
    return model
end

function build_model(
    sys::SystemData,
    model_type::Type{<:OPFModel};
    optimizer = nothing,
)
    opt = something(optimizer, _select_default_optimizer(model_type))
    times = get_model_times(sys)

    if length(times) == 1
        return _build_singleperiod_model(sys, model_type, opt)
    end

    return _build_multiperiod_model(sys, model_type, opt, times)
end

# ============================================================
# BUILDER
# ============================================================

function build_model(inputs::CentralizedInputs; optimizer = nothing)
    sys = inputs.system
    model_type = inputs.model_type

    opt = _resolve_optimizer(optimizer, inputs.optimizer, model_type)

    return build_model(sys, model_type; optimizer = opt)
end
