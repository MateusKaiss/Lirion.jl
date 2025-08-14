abstract type OPFModel end

struct DC   <: OPFModel end
struct AC   <: OPFModel end
struct SOCP <: OPFModel end

mutable struct DistributedOPFModel
    base_model::Type{<:OPFModel}
    method::Symbol
    max_iter::Int
    rho_theta::Float64   # ρθ 
    rho_v::Float64       # ρV 
    ε_primal::Float64
    ε_dual::Float64
end

DistributedOPFModel(
    base_model::Type{<:OPFModel},
    method::Symbol,
    max_iter::Int,
    rho::Float64,
    ε_primal::Float64,
    ε_dual::Float64,
) = DistributedOPFModel(base_model, method, max_iter, rho, rho, ε_primal, ε_dual)

abstract type SubproblemAction end

struct SubproblemBuild                    <: SubproblemAction end
struct SubproblemUpdate                   <: SubproblemAction end
struct SubproblemGetResult                <: SubproblemAction end
struct SubproblemNameElementsDimensions   <: SubproblemAction end
struct SubproblemWriteResult               <: SubproblemAction end