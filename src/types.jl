abstract type OPFModel end
struct DC <: OPFModel end
struct AC <: OPFModel end
struct SOCP <: OPFModel end

mutable struct DistributedOPFModel
    base_model::Type{<:OPFModel}
    method::Symbol
    max_iter::Int
    rho_theta::Float64     # rho_theta (consenso de ângulo)     [AC-bus; DC-bus]
    rho_v::Float64         # rhoV  (consenso de |V| ou V²) [AC/SOCP-bus]
    rho_flow_p::Float64    # rhoP  (consenso de P_flow / Re(W))
    rho_flow_q::Float64    # rhoQ  (consenso de Q_flow / Im(W))
    eps_primal::Float64
    eps_dual::Float64
    alpha_overrelax::Union{Nothing, Float64}
    warm_start::Bool
    adaptive_rho::Bool
    mu_incr::Float64
    mu_decr::Float64
    tau_incr::Float64
    tau_decr::Float64
    accel_nesterov::Bool
    optimizer::Any                  # pode ser nothing ou, ex., SCS.Optimizer
    prefer_socp_solver::Symbol      # :auto | :gurobi | :cosmo | :scs

    function DistributedOPFModel(
        base_model::Type{<:OPFModel},
        method::Symbol,
        max_iter::Integer,
        rho_theta::Real,
        rho_v::Real,
        rho_flow_p::Real,
        rho_flow_q::Real,
        eps_primal::Real,
        eps_dual::Real;
        alpha_overrelax::Union{Nothing, Real} = nothing,
        warm_start::Bool = false,
        adaptive_rho::Bool = true,
        mu_incr::Real = 3.0,
        mu_decr::Real = 3.0,
        tau_incr::Real = 1.05,
        tau_decr::Real = 0.98,
        accel_nesterov::Bool = false,
        optimizer::Any = nothing,
        prefer_socp_solver::Symbol = :auto,
    )
        max_iter > 0 || throw(ArgumentError("max_iter must be > 0"))
        rho_theta > 0 || throw(ArgumentError("rho_theta must be > 0"))
        rho_v > 0 || throw(ArgumentError("rho_v must be > 0"))
        rho_flow_p > 0 || throw(ArgumentError("rho_flow_p must be > 0"))
        rho_flow_q > 0 || throw(ArgumentError("rho_flow_q must be > 0"))
        eps_primal ≥ 0 || throw(ArgumentError("eps_primal must be ≥ 0"))
        eps_dual ≥ 0 || throw(ArgumentError("eps_dual must be ≥ 0"))
        mu_incr > 0 || throw(ArgumentError("mu_incr must be > 0"))
        mu_decr > 0 || throw(ArgumentError("mu_decr must be > 0"))
        tau_incr > 0 || throw(ArgumentError("tau_incr must be > 0"))
        tau_decr > 0 || throw(ArgumentError("tau_decr must be > 0"))

        α = isnothing(alpha_overrelax) ? nothing : Float64(alpha_overrelax)

        return new(
            base_model, method, Int(max_iter),
            Float64(rho_theta), Float64(rho_v),
            Float64(rho_flow_p), Float64(rho_flow_q),
            Float64(eps_primal), Float64(eps_dual), α,
            warm_start,
            adaptive_rho,
            Float64(mu_incr), Float64(mu_decr), Float64(tau_incr), Float64(tau_decr),
            accel_nesterov,
            optimizer,
            prefer_socp_solver,
        )
    end
end

# ----- construtor compacto com rho único -----
DistributedOPFModel(
    base_model::Type{<:OPFModel},
    method::Symbol,
    max_iter::Int,
    rho::Float64,
    eps_primal::Float64,
    eps_dual::Float64;
    alpha_overrelax::Union{Nothing, Float64} = nothing,
    warm_start::Bool = false,
    adaptive_rho::Bool = true,
    mu_incr::Float64 = 3.0,
    mu_decr::Float64 = 3.0,
    tau_incr::Float64 = 1.05,
    tau_decr::Float64 = 0.98,
    accel_nesterov::Bool = false,
    optimizer::Any = nothing,
    prefer_socp_solver::Symbol = :auto,
) = DistributedOPFModel(
    base_model, method, max_iter, rho, rho, rho, rho, eps_primal, eps_dual;
    alpha_overrelax = alpha_overrelax,
    warm_start = warm_start,
    adaptive_rho = adaptive_rho,
    mu_incr = mu_incr,
    mu_decr = mu_decr,
    tau_incr = tau_incr,
    tau_decr = tau_decr,
    accel_nesterov = accel_nesterov,
    optimizer = optimizer,
    prefer_socp_solver = prefer_socp_solver,
)

# ----- construtor keyword-friendly (permite sobrepor rhos específicos) -----
DistributedOPFModel(
    base_model::Type{<:OPFModel},
    method::Symbol,
    max_iter::Int;
    rho::Float64 = 1_000.0,
    rho_theta::Union{Nothing, Float64} = nothing,
    rho_v::Union{Nothing, Float64} = nothing,
    rho_flow::Union{Nothing, Float64} = nothing,
    rho_flow_p::Union{Nothing, Float64} = nothing,
    rho_flow_q::Union{Nothing, Float64} = nothing,
    eps_primal::Float64 = 1e-7,
    eps_dual::Float64 = 1e-7,
    alpha_overrelax::Union{Nothing, Float64} = nothing,
    warm_start::Bool = false,
    adaptive_rho::Bool = true,
    mu_incr::Float64 = 3.0,
    mu_decr::Float64 = 3.0,
    tau_incr::Float64 = 1.05,
    tau_decr::Float64 = 0.98,
    accel_nesterov::Bool = false,
    optimizer::Any = nothing,
    prefer_socp_solver::Symbol = :auto,
) = DistributedOPFModel(
    base_model, method, max_iter,
    coalesce(rho_theta, rho),
    coalesce(rho_v, rho),
    coalesce(rho_flow_p, rho_flow, rho),
    coalesce(rho_flow_q, rho_flow, rho),
    eps_primal, eps_dual;
    alpha_overrelax = alpha_overrelax,
    warm_start = warm_start,
    adaptive_rho = adaptive_rho,
    mu_incr = mu_incr,
    mu_decr = mu_decr,
    tau_incr = tau_incr,
    tau_decr = tau_decr,
    accel_nesterov = accel_nesterov,
    optimizer = optimizer,
    prefer_socp_solver = prefer_socp_solver,
)

@inline is_distributed_model(::Type{<:OPFModel}) = false
@inline is_distributed_model(::DistributedOPFModel) = true
@inline is_centralized_model(x) = !is_distributed_model(x)

# ---------- SHIM ----------
rho(c::DistributedOPFModel) = c.rho_theta
const rho = rho

Base.getproperty(c::DistributedOPFModel, s::Symbol) =
    s === :rho ? c.rho_theta :
    s === :rho_flow ? c.rho_flow_p :
    getfield(c, s)

Base.setproperty!(c::DistributedOPFModel, s::Symbol, v) =
    s === :rho ? (c.rho_theta = v;
    c.rho_v = v;
    c.rho_flow_p = v;
    c.rho_flow_q = v) :
    s === :rho_flow ? (c.rho_flow_p = v;
    c.rho_flow_q = v) :
    Base.setfield!(c, s, v)

Base.propertynames(c::DistributedOPFModel, private::Bool = false) =
    (propertynames(typeof(c), private)..., :rho, :rho_flow)

abstract type SubproblemAction end
struct SubproblemBuild <: SubproblemAction end
struct SubproblemUpdate <: SubproblemAction end
struct SubproblemGetResult <: SubproblemAction end
struct SubproblemNameElementsDimensions <: SubproblemAction end
struct SubproblemWriteResult <: SubproblemAction end
