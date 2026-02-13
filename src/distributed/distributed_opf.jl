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
# Base model -> Tag type
# ------------------------------------------------------------
@inline function _base_tag_type(::Type{B}) where {B <: OPFModel}
    B === DC && return DCTag
    B === AC && return ACTag
    B === SOCP && return SOCPTag
    return error("Unknown base_model: $(B). Expected DC, AC, or SOCP.")
end

# ------------------------------------------------------------
# ADMM tag (type-level) used by run_admm dispatch
# ------------------------------------------------------------
@inline function admm_tag(cfg::DistributedOPFModel)
    V = _variant_type(Val(cfg.coupling))
    B = _base_tag_type(cfg.base_model)
    return ADMM{V, B}()
end

# ------------------------------------------------------------
# Distributed runner
# ------------------------------------------------------------
function run_distributed(inputs::LirionInputs)::LirionRunResult
    @assert inputs.model_type isa DistributedOPFModel
    @assert inputs.admm_system !== nothing "Distributed run requires inputs.admm_system"

    cfg = inputs.model_type::DistributedOPFModel
    sys = inputs.admm_system::ADMMSystem

    @info "Running distributed OPF"
    @info "  coupling   = $(cfg.coupling)"
    @info "  sectioning = $(cfg.sectioning)"
    @info "  base_model = $(cfg.base_model)"
    @info "  optimizer  = $(inputs.optimizer)"

    optimizer_factory = resolve_optimizer(cfg, inputs.optimizer)

    output = run_admm(sys, cfg, admm_tag(cfg); optimizer_factory = optimizer_factory)

    if output isa LirionRunResult
        r = output.result
        base_nt = r isa NamedTuple ? r : (result = r,)

        meta_nt = (
            config = cfg,
            system = sys,
            database = inputs.database_path,
        )

        enriched = merge(base_nt, meta_nt)

        return LirionRunResult(;
            result = enriched,
            model = output.model,
            system = output.system,
            diagnostics = output.diagnostics,
        )
    end

    return error(
        "Runner returned unsupported type: $(typeof(output)). Expected LirionRunResult.",
    )
end
