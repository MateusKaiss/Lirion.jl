
function get_method_runner(model::DistributedOPFModel)
    if model.method == :admm
        return run_admm
    else
        error("Distributed method $(model.method) not implemented.")
    end
end
