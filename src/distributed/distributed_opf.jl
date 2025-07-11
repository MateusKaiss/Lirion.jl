
function get_method_runner(model::DistributedOPFModel)
    if model.method == :admm
        return ADMM.run_admm
    else
        error("Método distribuído $(model.method) não implementado.")
    end
end