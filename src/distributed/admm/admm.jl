module ADMM

function run_admm(model_config, path)
    println("Executando ADMM")
    println("Base model: ", model_config.base_model)
    println("Iterações: ", model_config.max_iter)
    println("Rho: ", model_config.rho)
    println("Dados: ", path)
    return "finalizado"
end

end
