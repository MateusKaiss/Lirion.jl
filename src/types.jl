abstract type OPFModel end
struct DC <: OPFModel end
struct AC <: OPFModel end
struct SOCP <: OPFModel end

mutable struct DistributedOPFModel
    base_model::Type{<:OPFModel}
    method::Symbol
    max_iter::Int
    rho_theta::Float64    # ρθ  (consenso em ângulos / DC-bus e AC)
    rho_v::Float64        # ρV  (consenso em |V| / AC)
    rho_flow::Float64     # ρe  (consenso em fluxo / branch & hybrid)
    ε_primal::Float64
    ε_dual::Float64

    function DistributedOPFModel(
        base_model::Type{<:OPFModel},
        method::Symbol,
        max_iter::Integer,
        rho_theta::Real,
        rho_v::Real,
        rho_flow::Real,
        ε_primal::Real,
        ε_dual::Real,
    )
        max_iter > 0 || throw(ArgumentError("max_iter must be > 0"))
        rho_theta > 0 || throw(ArgumentError("rho_theta must be > 0"))
        rho_v > 0 || throw(ArgumentError("rho_v must be > 0"))
        rho_flow > 0 || throw(ArgumentError("rho_flow must be > 0"))
        ε_primal >= 0 || throw(ArgumentError("ε_primal must be ≥ 0"))
        ε_dual >= 0 || throw(ArgumentError("ε_dual must be ≥ 0"))

        return new(
            base_model,
            method,
            Int(max_iter),
            Float64(rho_theta),
            Float64(rho_v),
            Float64(rho_flow),
            Float64(ε_primal),
            Float64(ε_dual),
        )
    end
end

DistributedOPFModel(
    base_model::Type{<:OPFModel},
    method::Symbol,
    max_iter::Int,
    rho::Float64,
    ε_primal::Float64,
    ε_dual::Float64,
) = DistributedOPFModel(base_model, method, max_iter, rho, rho, rho, ε_primal, ε_dual)

DistributedOPFModel(
    base_model::Type{<:OPFModel},
    method::Symbol,
    max_iter::Int;
    rho::Float64 = 1_000.0,
    rho_theta::Union{Nothing, Float64} = nothing,
    rho_v::Union{Nothing, Float64} = nothing,
    rho_flow::Union{Nothing, Float64} = nothing,
    ε_primal::Float64 = 1e-7,
    ε_dual::Float64 = 1e-7,
) = DistributedOPFModel(
    base_model, method, max_iter,
    coalesce(rho_theta, rho),
    coalesce(rho_v, rho),
    coalesce(rho_flow, rho),
    ε_primal, ε_dual,
)

# ---------- SHIM  ----------
rho(c::DistributedOPFModel) = c.rho_theta
const ρ = rho

Base.getproperty(c::DistributedOPFModel, s::Symbol) =
    s === :rho ? c.rho_theta : getfield(c, s)

Base.setproperty!(c::DistributedOPFModel, s::Symbol, v) =
    s === :rho ? (c.rho_theta = v; c.rho_v = v; c.rho_flow = v) : Base.setfield!(c, s, v)

Base.propertynames(c::DistributedOPFModel, private::Bool = false) =
    (propertynames(typeof(c), private)..., :rho)

abstract type SubproblemAction end

struct SubproblemBuild <: SubproblemAction end
struct SubproblemUpdate <: SubproblemAction end
struct SubproblemGetResult <: SubproblemAction end
struct SubproblemNameElementsDimensions <: SubproblemAction end
struct SubproblemWriteResult <: SubproblemAction end