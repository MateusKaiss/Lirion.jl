module DistributedOPF

export DistributedOPFModel, get_method_runner

const PowerOptimization = Main.PowerOptimization

import ..PowerOptimization: OPFModel

struct DistributedOPFModel <: OPFModel
    base_model::Type{<:OPFModel}  # DC, AC, SOCP
    method::Symbol                #:admm, :benders, etc.
    max_iter::Int
    rho::Float64
end

function get_method_runner(model::DistributedOPFModel)
    if model.method == :admm
        return ADMM.run_admm
    else
        error("Método distribuído $(model.method) não implementado.")
    end
end

include("admm/admm.jl")
using .ADMM

end
