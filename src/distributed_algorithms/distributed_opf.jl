# ------------------------------------------------------------
# Coupling -> Variant type
# ------------------------------------------------------------
@inline function _variant_type(::Val{c}) where {c}
    c === :state && return State
    c === :power && return Power
    c === :hybrid && return Hybrid
    return error("Unknown coupling: $(c). Expected :state, :power, or :hybrid.")
end

# ------------------------------------------------------------
# Algorithm display name 
# ------------------------------------------------------------
_algo_name(::ADMMConfig) = "ADMM"
_algo_name(::ALADINConfig) = "ALADIN"

# ------------------------------------------------------------
# Distributed runner — single generic entry point
# ------------------------------------------------------------
function run_distributed(inputs::DistributedInputs)::LirionResult
    cfg = inputs.model_type
    sys = inputs.system

    @info "Running distributed OPF ($(_algo_name(cfg)))"
    @info "  coupling   = $(cfg.coupling)"
    @info "  sectioning = $(cfg.sectioning)"
    @info "  base_model = $(cfg.base_model)"
    @info "  optimizer  = $(inputs.optimizer)"

    optimizer_factory = resolve_optimizer(cfg, inputs.optimizer)
    return run_algorithm(sys, cfg; optimizer_factory = optimizer_factory)
end
