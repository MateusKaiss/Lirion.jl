function get_method_runner(model::DistributedOPFModel)
    if model.method == :admm && model.base_model == DC
        return run_admm
    elseif model.method == :admm && model.base_model == AC
        return run_admm_ac
    elseif model.method == :admm_branch && model.base_model == DC
        return run_admm_branch
    elseif model.method == :admm_hybrid && model.base_model == DC
        return run_admm_hybrid
    else
        error(
            "Distributed method $(model.method) with base $(model.base_model) not implemented.",
        )
    end
end