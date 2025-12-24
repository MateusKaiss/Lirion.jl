@inline get_method_runner(config::DistributedOPFModel) =
    get_method_runner(config.base_model, Val(config.method))

# DC
@inline get_method_runner(::Type{DC}, ::Val{:admm}) = run_admm
@inline get_method_runner(::Type{DC}, ::Val{:admm_branch}) = run_admm_branch
@inline get_method_runner(::Type{DC}, ::Val{:admm_hybrid}) = run_admm_hybrid

# AC
@inline get_method_runner(::Type{AC}, ::Val{:admm}) = run_admm_ac
@inline get_method_runner(::Type{AC}, ::Val{:admm_branch}) = run_admm_ac_branch
@inline get_method_runner(::Type{AC}, ::Val{:admm_hybrid}) = run_admm_ac_hybrid

# SOCP
@inline get_method_runner(::Type{SOCP}, ::Val{:admm}) = run_admm_socp
@inline get_method_runner(::Type{SOCP}, ::Val{:admm_branch}) = run_admm_socp_branch
@inline get_method_runner(::Type{SOCP}, ::Val{:admm_hybrid}) = run_admm_socp_hybrid

function get_method_runner(base::Type{<:OPFModel}, ::Val{m}) where {m}
    return error("Distributed method $(m) with base model $(base) not implemented.")
end

function run_distributed(inputs::LirionInputs)::LirionRunResult
    @assert inputs.model_type isa DistributedOPFModel

    config = inputs.model_type::DistributedOPFModel
    admm_sys = inputs.admm_system::ADMMSystem

    @info "Running distributed OPF"
    @info "  method     = $(config.method)"
    @info "  base_model = $(config.base_model)"

    runner = get_method_runner(config)
    output = runner(admm_sys, config)

    if output isa LirionRunResult
        r = output.result
        base_nt = r isa NamedTuple ? r : (result = r,)
        meta_nt = (config = config, system = inputs.system, database = inputs.database_path)
        enriched = merge(base_nt, meta_nt)

        return LirionRunResult(;
            result = enriched,
            model = output.model,
            system = output.system === nothing ? inputs.system : output.system,
            diagnostics = output.diagnostics,
        )
    end

    return error(
        "Runner returned unsupported type: $(typeof(output)). Expected LirionRunResult.",
    )
end
