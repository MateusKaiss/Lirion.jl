abstract type OPFModel end

struct DC <: OPFModel end
struct AC <: OPFModel end
struct SOCP <: OPFModel end

struct DistributedOPFModel <: OPFModel
    base_model::Type{<:OPFModel}  # DC, AC, SOCP
    method::Symbol                #:admm, :benders, etc.
    max_iter::Int
    rho::Float64
end

abstract type SubproblemAction end

struct SubproblemBuild <: SubproblemAction end
struct SubproblemUpdate <: SubproblemAction end
struct SubproblemGetResult <: SubproblemAction end
struct SubproblemNameElementsDimensions <: SubproblemAction end
struct SubproblemWriteResult <: SubproblemAction end
