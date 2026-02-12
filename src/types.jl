# ============================================================
# OPF BASE MODELS
# ============================================================

abstract type OPFModel end
struct DC <: OPFModel end
struct AC <: OPFModel end
struct SOCP <: OPFModel end

# ============================================================
# ADMM VARIANTS
# ============================================================

abstract type ADMMVariant end
struct Bus <: ADMMVariant end
struct Branch <: ADMMVariant end
struct Hybrid <: ADMMVariant end

# ============================================================
# OPF TAGS FOR DISTRIBUTED VARIANTS
# ============================================================

abstract type OPFTag end
struct DCTag <: OPFTag end
struct ACTag <: OPFTag end
struct SOCPTag <: OPFTag end

# ============================================================
# ADMM TYPE COMBINATOR
# ============================================================

struct ADMM{V <: ADMMVariant, B <: OPFTag} end

# ============================================================
# DISTRIBUTED OPF MODEL CONFIGURATION
# ============================================================

mutable struct DistributedOPFModel
    base_model::Type{<:OPFModel}
    method::Symbol
    max_iter::Int

    # bus/state rhos
    rho_theta::Float64
    rho_v::Float64           # (SOCP) você pode enxergar como rho(V²)

    # branch/flow rhos (separa semanticamente)
    rho_p::Float64
    rho_q::Float64
    rho_r::Float64
    rho_t::Float64

    eps_abs::Float64
    eps_rel::Float64

    adaptive_rho::Bool
    mu_incr::Float64
    mu_decr::Float64
    tau_incr::Float64
    tau_decr::Float64

    function DistributedOPFModel(
        base_model::Type{<:OPFModel},
        method::Symbol,
        max_iter::Integer, rho_theta::Real,
        rho_v::Real,
        rho_p::Real,
        rho_q::Real,
        rho_r::Real,
        rho_t::Real, eps_abs::Real,
        eps_rel::Real;
        adaptive_rho::Bool = true,
        mu_incr::Real = 3.0,
        mu_decr::Real = 3.0,
        tau_incr::Real = 1.05,
        tau_decr::Real = 0.98,
    )
        max_iter > 0 || throw(ArgumentError("max_iter must be > 0"))
        rho_theta > 0 || throw(ArgumentError("rho_theta must be > 0"))
        rho_v > 0 || throw(ArgumentError("rho_v must be > 0"))
        rho_p > 0 || throw(ArgumentError("rho_p must be > 0"))
        rho_q > 0 || throw(ArgumentError("rho_q must be > 0"))
        rho_r > 0 || throw(ArgumentError("rho_r must be > 0"))
        rho_t > 0 || throw(ArgumentError("rho_t must be > 0"))
        eps_abs ≥ 0 || throw(ArgumentError("eps_abs must be ≥ 0"))
        eps_rel ≥ 0 || throw(ArgumentError("eps_rel must be ≥ 0"))
        mu_incr > 0 || throw(ArgumentError("mu_incr must be > 0"))
        mu_decr > 0 || throw(ArgumentError("mu_decr must be > 0"))
        tau_incr > 0 || throw(ArgumentError("tau_incr must be > 0"))
        tau_decr > 0 || throw(ArgumentError("tau_decr must be > 0"))

        return new(
            base_model, method, Int(max_iter), Float64(rho_theta),
            Float64(rho_v), Float64(rho_p),
            Float64(rho_q),
            Float64(rho_r),
            Float64(rho_t), Float64(eps_abs),
            Float64(eps_rel), adaptive_rho,
            Float64(mu_incr),
            Float64(mu_decr),
            Float64(tau_incr),
            Float64(tau_decr),
        )
    end
end

# ============================================================
# DISTRIBUTED OPF MODEL CONVENIENCE CONSTRUCTORS
# ============================================================

# (A) Posicional (seta tudo igual)
DistributedOPFModel(
    base_model::Type{<:OPFModel},
    method::Symbol,
    max_iter::Int,
    rho::Float64,
    eps_abs::Float64,
    eps_rel::Float64;
    adaptive_rho::Bool = true,
    mu_incr::Float64 = 3.0,
    mu_decr::Float64 = 3.0,
    tau_incr::Float64 = 1.05,
    tau_decr::Float64 = 0.98,
) = DistributedOPFModel(
    base_model, method, max_iter,
    rho, rho, rho, rho, rho, rho,
    eps_abs, eps_rel;
    adaptive_rho = adaptive_rho,
    mu_incr = mu_incr,
    mu_decr = mu_decr,
    tau_incr = tau_incr,
    tau_decr = tau_decr,
)

# (B) kwargs flexível (rho + overrides)
# Precedência:
#   rho_theta/rho_v/rho_p/rho_q/rho_r/rho_t > rho_flow_pq / rho_flow_rt / rho_flow > rho
DistributedOPFModel(
    base_model::Type{<:OPFModel},
    method::Symbol,
    max_iter::Int;
    rho::Float64 = 1_000.0, rho_theta::Union{Nothing, Float64} = nothing,
    rho_v::Union{Nothing, Float64} = nothing, rho_flow::Union{Nothing, Float64} = nothing,      # fallback para p,q,r,t
    rho_flow_pq::Union{Nothing, Float64} = nothing,   # fallback para p,q
    rho_flow_rt::Union{Nothing, Float64} = nothing,   # fallback para r,t

    rho_p::Union{Nothing, Float64} = nothing,
    rho_q::Union{Nothing, Float64} = nothing,
    rho_r::Union{Nothing, Float64} = nothing,
    rho_t::Union{Nothing, Float64} = nothing,

    rho_flow_p::Union{Nothing, Float64} = nothing,    # alias -> rho_p
    rho_flow_q::Union{Nothing, Float64} = nothing,    # alias -> rho_q
    rho_flow_r::Union{Nothing, Float64} = nothing,    # alias -> rho_r
    rho_flow_t::Union{Nothing, Float64} = nothing,    # alias -> rho_t
    eps_abs::Float64 = 1e-7,
    eps_rel::Float64 = 1e-7, adaptive_rho::Bool = true,
    mu_incr::Float64 = 3.0,
    mu_decr::Float64 = 3.0,
    tau_incr::Float64 = 1.05,
    tau_decr::Float64 = 0.98,
) = DistributedOPFModel(
    base_model, method, max_iter, something(rho_theta, rho), something(rho_v, rho),

    # p/q (aceita rho_p ou rho_flow_p; e fallback em rho_flow_pq/rho_flow/rho)
    something(rho_p, rho_flow_p, rho_flow_pq, rho_flow, rho),
    something(rho_q, rho_flow_q, rho_flow_pq, rho_flow, rho),

    # r/t (aceita rho_r ou rho_flow_r; e fallback em rho_flow_rt/rho_flow/rho)
    something(rho_r, rho_flow_r, rho_flow_rt, rho_flow, rho),
    something(rho_t, rho_flow_t, rho_flow_rt, rho_flow, rho), eps_abs, eps_rel;
    adaptive_rho = adaptive_rho,
    mu_incr = mu_incr,
    mu_decr = mu_decr,
    tau_incr = tau_incr,
    tau_decr = tau_decr,
)

# ============================================================
# MODEL TYPE CHECKS
# ============================================================

@inline is_distributed_model(::Type{<:OPFModel}) = false
@inline is_distributed_model(::DistributedOPFModel) = true
@inline is_centralized_model(x) = !is_distributed_model(x)

# ============================================================
# RHO SHIMS (compat + rho_v2)
# ============================================================

# rho() continua sendo rho_theta por padrão (como você tinha)
rho(c::DistributedOPFModel) = c.rho_theta
const rho = rho

# getproperty: mantém nomes antigos E adiciona rho_v2
Base.getproperty(c::DistributedOPFModel, s::Symbol) =
    s === :rho ? c.rho_theta :
    s === :rho_flow ? c.rho_p :          # "flow genérico" -> p (e set vai setar p,q,r,t)
    s === :rho_flow_p ? c.rho_p :
    s === :rho_flow_q ? c.rho_q :
    s === :rho_flow_r ? c.rho_r :
    s === :rho_flow_t ? c.rho_t :
    s === :rho_v2 ? c.rho_v :          # <-- pedido: alias pra V²
    getfield(c, s)

# setproperty: mantém nomes antigos E adiciona set de rho_v2
Base.setproperty!(c::DistributedOPFModel, s::Symbol, v) =
    s === :rho ? (c.rho_theta = v;
    c.rho_v = v;
    c.rho_p = v;
    c.rho_q = v;
    c.rho_r = v;
    c.rho_t = v) :
    s === :rho_flow ? (c.rho_p = v;
    c.rho_q = v;
    c.rho_r = v;
    c.rho_t = v) :
    s === :rho_flow_p ? (c.rho_p = v) :
    s === :rho_flow_q ? (c.rho_q = v) :
    s === :rho_flow_r ? (c.rho_r = v) :
    s === :rho_flow_t ? (c.rho_t = v) :
    s === :rho_v2 ? (c.rho_v = v) :
    Base.setfield!(c, s, v)

Base.propertynames(c::DistributedOPFModel, private::Bool = false) =
    (propertynames(typeof(c), private)...,
        :rho, :rho_flow, :rho_flow_p, :rho_flow_q, :rho_flow_r, :rho_flow_t,
        :rho_v2)

# ============================================================
# RUN DIAGNOSTICS
# ============================================================

abstract type AbstractRunDiagnostics end
struct NoDiagnostics <: AbstractRunDiagnostics end

Base.@kwdef struct ADMMHistory <: AbstractRunDiagnostics
    r_primal::Vector{Float64} = Float64[]
    r_dual::Vector{Float64} = Float64[]
    eps_primal::Vector{Float64} = Float64[]
    eps_dual::Vector{Float64} = Float64[]
    rhos::Dict{Symbol, Vector{Float64}} = Dict{Symbol, Vector{Float64}}()
    heatmaps::Dict{Symbol, Any} = Dict{Symbol, Any}()
    meta::Dict{Symbol, Any} = Dict{Symbol, Any}()
end

# ============================================================
# RUN RESULT
# ============================================================

Base.@kwdef struct LirionRunResult
    result::Any = nothing
    model::Any = nothing
    system::Any = nothing
    diagnostics::AbstractRunDiagnostics = NoDiagnostics()
end
