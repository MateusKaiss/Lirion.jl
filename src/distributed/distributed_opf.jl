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

function run_distributed(inputs::LirionInputs)
    @assert inputs.model_type isa DistributedOPFModel

    config = inputs.model_type::DistributedOPFModel
    admm_sys = inputs.admm_system

    @info "Running distributed OPF"
    @info "  method     = $(config.method)"
    @info "  base_model = $(config.base_model)"

    runner = get_method_runner(config)
    total_real_cost, partitions, iter_done, elapsed_s = runner(config, admm_sys)

    return (
        total_cost = total_real_cost,
        partitions = partitions,
        iterations = iter_done,
        elapsed_s = elapsed_s,
        config = config,
        system = inputs.system,
        database = inputs.database_path,
    )
end
