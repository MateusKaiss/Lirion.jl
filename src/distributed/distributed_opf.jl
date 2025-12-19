function get_method_runner(config::DistributedOPFModel)::Function
    method = config.method
    base = config.base_model

    if method === :admm
        if base === DC
            return run_admm
        elseif base === AC
            return run_admm_ac
        elseif base === SOCP
            return run_admm_socp
        end

    elseif method === :admm_branch
        if base === DC
            return run_admm_branch
        elseif base === AC
            return run_admm_ac_branch
        elseif base === SOCP
            return run_admm_socp_branch
        end

    elseif method === :admm_hybrid
        if base === DC
            return run_admm_hybrid
        elseif base === AC
            return run_admm_ac_hybrid
        elseif base === SOCP
            return run_admm_socp_hybrid
        end
    end

    return error("Distributed method $(method) with base model $(base) not implemented.")
end

"""
Executa o método distribuído e retorna SEMPRE `LirionRunResult`.

Compatibilidade:
- Se o runner retornar `LirionRunResult` (novo padrão), usamos `out.result`.
- Se retornar a tupla legada `(total_cost, partitions, iter_done, elapsed_s)`, convertemos.
"""
function run_distributed(inputs::LirionInputs)::LirionRunResult
    @assert inputs.model_type isa DistributedOPFModel

    config = inputs.model_type::DistributedOPFModel
    admm_sys = inputs.admm_system

    @info "Running distributed OPF"
    @info "  method     = $(config.method)"
    @info "  base_model = $(config.base_model)"

    runner = get_method_runner(config)
    out = runner(config, admm_sys)

    # -----------------------------
    # NOVO PADRÃO: retorna LirionRunResult
    # -----------------------------
    if out isa LirionRunResult
        # Garanta que o `system`/`database` apareçam no result final (opcional, mas útil)
        # (não mexe no diagnostics)
        r = out.result
        enriched = merge(
            r isa NamedTuple ? r : (result = r,),
            (config = config,
                system = inputs.system,
                database = inputs.database_path),
        )
        return LirionRunResult(;
            result = enriched,
            model = out.model,
            system = out.system === nothing ? inputs.system : out.system,
            diagnostics = out.diagnostics,
        )
    end

    # -----------------------------
    # LEGADO: retorna tupla (total_cost, partitions, iter_done, elapsed_s)
    # -----------------------------
    if out isa Tuple && length(out) == 4
        total_real_cost, partitions, iter_done, elapsed_s = out
        result = (
            total_cost = total_real_cost,
            partitions = partitions,
            iterations = iter_done,
            elapsed_s = elapsed_s,
            config = config,
            system = inputs.system,
            database = inputs.database_path,
        )
        return LirionRunResult(;
            result = result,
            model = nothing,
            system = inputs.system,
            diagnostics = NoDiagnostics(),
        )
    end

    return error(
        "Runner returned unsupported type: $(typeof(out)). Expected LirionRunResult or a 4-tuple.",
    )
end
