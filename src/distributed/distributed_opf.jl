@inline function _variant_type(::Val{m}) where {m}
    m === :admm && return Bus
    m === :admm_branch && return Branch
    m === :admm_hybrid && return Hybrid
    return error(
        "Unknown ADMM method: $(m). Expected :admm, :admm_branch, or :admm_hybrid.",
    )
end

@inline function _base_tag_type(::Type{B}) where {B <: OPFModel}
    B === DC && return DCTag
    B === AC && return ACTag
    B === SOCP && return SOCPTag
    return error("Unknown base_model: $(B). Expected DC, AC, or SOCP.")
end

@inline function admm_tag(cfg::DistributedOPFModel)
    V = _variant_type(Val(cfg.method))
    B = _base_tag_type(cfg.base_model)
    return ADMM{V, B}()
end

function run_distributed(inputs::LirionInputs)::LirionRunResult
    @assert inputs.model_type isa DistributedOPFModel
    @assert inputs.admm_system !== nothing "Distributed run requires inputs.admm_system"

    cfg = inputs.model_type::DistributedOPFModel
    sys = inputs.admm_system::ADMMSystem

    @info "Running distributed OPF"
    @info "  method     = $(cfg.method)"
    @info "  base_model = $(cfg.base_model)"

    output = run_admm(sys, cfg, admm_tag(cfg))

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
