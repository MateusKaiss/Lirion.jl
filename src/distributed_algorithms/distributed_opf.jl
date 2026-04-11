# ------------------------------------------------------------
# Coupling -> Variant type
# ------------------------------------------------------------
@inline function _variant_type(::Val{c}) where {c}
    c === :state && return State
    c === :power && return Power
    c === :hybrid && return Hybrid
    return error("Unknown ADMM coupling: $(c). Expected :state, :power, or :hybrid.")
end

# ------------------------------------------------------------
# Distributed runner
# ------------------------------------------------------------
function run_distributed(inputs::DistributedInputs)::LirionRunResult
    cfg = inputs.model_type
    sys = inputs.system

    @info "Running distributed OPF"
    @info "  coupling   = $(cfg.coupling)"
    @info "  sectioning = $(cfg.sectioning)"
    @info "  base_model = $(cfg.base_model)"
    @info "  optimizer  = $(inputs.optimizer)"

    optimizer_factory = resolve_optimizer(cfg, inputs.optimizer)

    output = run_admm(sys, cfg, admm_tag(cfg); optimizer_factory = optimizer_factory)

    r = output.result
    base_nt = r isa NamedTuple ? r : (result = r,)

    meta_nt = (
        config = cfg,
        system = sys,
        database = inputs.database_path,
    )

    return LirionRunResult(;
        result = merge(base_nt, meta_nt),
        model = output.model,
        system = output.system,
        diagnostics = output.diagnostics,
    )
end
