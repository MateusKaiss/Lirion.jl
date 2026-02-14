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

struct State <: ADMMVariant end
struct Power <: ADMMVariant end
struct Hybrid <: ADMMVariant end

# ============================================================
# OPF TAGS FOR DISTRIBUTED VARIANTS (mantém como estava)
# ============================================================

abstract type OPFTag end
struct DCTag <: OPFTag end
struct ACTag <: OPFTag end
struct SOCPTag <: OPFTag end

# ============================================================
# ADMM TYPE COMBINATOR (mantém)
# ============================================================

struct ADMM{V <: ADMMVariant, B <: OPFTag} end

# ============================================================
# Algorithm spec
# ============================================================

struct ADMMAlgorithmSpec
    coupling::Symbol        # :state | :power | :hybrid
    sectioning::Symbol      # :tie_line
    max_iter::Int
    kwargs::NamedTuple
end

# ============================================================
# RHO RULES
# ============================================================

const _RHO_RULES = Dict{Tuple{DataType, Symbol}, NamedTuple}(
    # ----------------------------
    # DC
    # ----------------------------
    (DC, :state) => (
        allowed = Set([:rho_theta]),
        required = Set([:rho_theta]),
    ),
    (DC, :power) => (
        allowed = Set([:rho_p]),
        required = Set([:rho_p]),
    ),
    (DC, :hybrid) => (
        allowed = Set([:rho_theta, :rho_p]),
        required = Set([:rho_theta, :rho_p]),
    ),

    # ----------------------------
    # AC
    # ----------------------------
    (AC, :state) => (
        allowed = Set([:rho_theta, :rho_v]),
        required = Set([:rho_theta, :rho_v]),
    ),
    (AC, :power) => (
        allowed = Set([:rho_p, :rho_q]),
        required = Set([:rho_p, :rho_q]),
    ),
    (AC, :hybrid) => (
        allowed = Set([:rho_theta, :rho_v, :rho_p, :rho_q]),
        required = Set([:rho_theta, :rho_v, :rho_p, :rho_q]),
    ),

    # ----------------------------
    # SOCP (SEM rho_theta; usa rho_v2 e exige rho_r/rho_t no state)
    # ----------------------------
    (SOCP, :state) => (
        allowed = Set([:rho_v2, :rho_r, :rho_t]),
        required = Set([:rho_v2, :rho_r, :rho_t]),
    ),
    (SOCP, :power) => (
        allowed = Set([:rho_p, :rho_q]),
        required = Set([:rho_p, :rho_q]),
    ),
    (SOCP, :hybrid) => (
        allowed = Set([:rho_v2, :rho_p, :rho_q, :rho_r, :rho_t]),
        required = Set([:rho_v2, :rho_p, :rho_q, :rho_r, :rho_t]),
    ),
)

@inline function _rules_for(base_model::Type{<:OPFModel}, coupling::Symbol)
    key = (base_model, coupling)
    haskey(_RHO_RULES, key) ||
        error("No rho rules for base_model=$(base_model), coupling=$(coupling)")
    return _RHO_RULES[key]
end

@inline function _asfloat(x, name::Symbol)
    x isa Real || throw(ArgumentError("$name must be Real, got $(typeof(x))"))
    xf = Float64(x)
    xf > 0 || throw(ArgumentError("$name must be > 0, got $xf"))
    return xf
end

function _normalize_rhos_strict(
    base_model::Type{<:OPFModel},
    coupling::Symbol;
    # generic
    rho::Union{Nothing, Real} = nothing,
    rho_theta::Union{Nothing, Real} = nothing,
    rho_v::Union{Nothing, Real} = nothing,      # AC only (forbidden in SOCP)
    rho_v2::Union{Nothing, Real} = nothing,     # SOCP only (stored in field rho_v)
    rho_p::Union{Nothing, Real} = nothing,
    rho_q::Union{Nothing, Real} = nothing,
    rho_r::Union{Nothing, Real} = nothing,
    rho_t::Union{Nothing, Real} = nothing,
)
    rules = _rules_for(base_model, coupling)
    allowed = rules.allowed
    required = rules.required

    # ----------------------------
    # Hard SOCP bans
    # ----------------------------
    if base_model === SOCP
        rho_theta === nothing ||
            error("SOCP does not accept rho_theta. Use rho_v2, rho_r, rho_t (or rho=...).")
        rho_v === nothing ||
            error("SOCP does not accept rho_v. Use rho_v2 (or rho=...).")
    end

    provided = Set{Symbol}()
    rho_theta !== nothing && push!(provided, :rho_theta)
    rho_v !== nothing && push!(provided, :rho_v)
    rho_v2 !== nothing && push!(provided, :rho_v2)
    rho_p !== nothing && push!(provided, :rho_p)
    rho_q !== nothing && push!(provided, :rho_q)
    rho_r !== nothing && push!(provided, :rho_r)
    rho_t !== nothing && push!(provided, :rho_t)

    for sym in provided
        sym in allowed || error(
            "Invalid rho parameter $sym for base_model=$(base_model), coupling=$(coupling). " *
            "Allowed: $(sort!(collect(allowed))). (You can also pass generic rho=...)",
        )
    end

    if rho === nothing
        missing = setdiff(required, provided)
        isempty(missing) || error(
            "Missing required rho parameters for base_model=$(base_model), coupling=$(coupling): " *
            "$(sort!(collect(missing))). Either provide them explicitly, or pass generic rho=... to fill required.",
        )
    end

    ρ_generic = rho === nothing ? nothing : _asfloat(rho, :rho)

    if ρ_generic !== nothing
        if :rho_theta in required
            rho_theta = something(rho_theta, ρ_generic)
        end
        if :rho_v in required
            rho_v = something(rho_v, ρ_generic)
        end
        if :rho_v2 in required
            rho_v2 = something(rho_v2, ρ_generic)
        end
        if :rho_p in required
            rho_p = something(rho_p, ρ_generic)
        end
        if :rho_q in required
            rho_q = something(rho_q, ρ_generic)
        end
        if :rho_r in required
            rho_r = something(rho_r, ρ_generic)
        end
        if :rho_t in required
            rho_t = something(rho_t, ρ_generic)
        end
    end

    function _get_required_value(sym::Symbol)
        sym === :rho_theta && return _asfloat(rho_theta, :rho_theta)
        sym === :rho_v && return _asfloat(rho_v, :rho_v)
        sym === :rho_v2 && return _asfloat(rho_v2, :rho_v2)
        sym === :rho_p && return _asfloat(rho_p, :rho_p)
        sym === :rho_q && return _asfloat(rho_q, :rho_q)
        sym === :rho_r && return _asfloat(rho_r, :rho_r)
        sym === :rho_t && return _asfloat(rho_t, :rho_t)
        return error("unexpected required sym: $sym")
    end

    req_sorted = sort!(collect(required))
    anchor_sym = req_sorted[1]
    anchor_val = _get_required_value(anchor_sym)

    ρθ =
        (base_model === SOCP) ? anchor_val :
        ((:rho_theta in required) ? _asfloat(rho_theta, :rho_theta) : anchor_val)

    ρv = if base_model === AC
        (:rho_v in required) ? _asfloat(rho_v, :rho_v) : anchor_val
    elseif base_model === SOCP
        (:rho_v2 in required) ? _asfloat(rho_v2, :rho_v2) : anchor_val
    else
        anchor_val
    end

    ρp = (:rho_p in required) ? _asfloat(rho_p, :rho_p) : anchor_val
    ρq = (:rho_q in required) ? _asfloat(rho_q, :rho_q) : anchor_val
    ρr = (:rho_r in required) ? _asfloat(rho_r, :rho_r) : anchor_val
    ρt = (:rho_t in required) ? _asfloat(rho_t, :rho_t) : anchor_val

    return (rho_theta = ρθ, rho_v = ρv, rho_p = ρp, rho_q = ρq, rho_r = ρr, rho_t = ρt)
end

# ============================================================
# DISTRIBUTED OPF MODEL CONFIGURATION
# ============================================================

mutable struct DistributedOPFModel
    base_model::Type{<:OPFModel}
    coupling::Symbol
    sectioning::Symbol
    max_iter::Int

    # state rhos
    rho_theta::Float64
    rho_v::Float64

    # power/flow rhos
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
        coupling::Symbol,
        sectioning::Symbol,
        max_iter::Integer,
        rho_theta::Real,
        rho_v::Real,
        rho_p::Real,
        rho_q::Real,
        rho_r::Real,
        rho_t::Real,
        eps_abs::Real,
        eps_rel::Real;
        adaptive_rho::Bool = true,
        mu_incr::Real = 3.0,
        mu_decr::Real = 3.0,
        tau_incr::Real = 1.05,
        tau_decr::Real = 0.98,
    )
        coupling in (:state, :power, :hybrid) ||
            throw(
                ArgumentError(
                    "Invalid coupling=$coupling. Allowed: :state, :power, :hybrid",
                ),
            )
        sectioning === :tie_line ||
            throw(
                ArgumentError(
                    "Invalid sectioning=$sectioning. Currently supported: :tie_line",
                ),
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
            base_model,
            coupling,
            sectioning,
            Int(max_iter),
            Float64(rho_theta),
            Float64(rho_v),
            Float64(rho_p),
            Float64(rho_q),
            Float64(rho_r),
            Float64(rho_t),
            Float64(eps_abs),
            Float64(eps_rel),
            adaptive_rho,
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

DistributedOPFModel(
    base_model::Type{<:OPFModel},
    coupling::Symbol,
    sectioning::Symbol,
    max_iter::Int,
    rho::Real,
    eps_abs::Real,
    eps_rel::Real;
    adaptive_rho::Bool = true,
    mu_incr::Real = 3.0,
    mu_decr::Real = 3.0,
    tau_incr::Real = 1.05,
    tau_decr::Real = 0.98,
) = begin
    rhos = _normalize_rhos_strict(base_model, coupling; rho = rho)
    DistributedOPFModel(
        base_model, coupling, sectioning, max_iter,
        rhos.rho_theta, rhos.rho_v, rhos.rho_p, rhos.rho_q, rhos.rho_r, rhos.rho_t,
        eps_abs, eps_rel;
        adaptive_rho = adaptive_rho,
        mu_incr = mu_incr,
        mu_decr = mu_decr,
        tau_incr = tau_incr,
        tau_decr = tau_decr,
    )
end

DistributedOPFModel(
    base_model::Type{<:OPFModel},
    coupling::Symbol,
    sectioning::Symbol,
    max_iter::Int;
    rho::Union{Nothing, Real} = nothing,
    rho_theta::Union{Nothing, Real} = nothing,
    rho_v::Union{Nothing, Real} = nothing,
    rho_v2::Union{Nothing, Real} = nothing,
    rho_p::Union{Nothing, Real} = nothing,
    rho_q::Union{Nothing, Real} = nothing,
    rho_r::Union{Nothing, Real} = nothing,
    rho_t::Union{Nothing, Real} = nothing,
    eps_abs::Real = 1e-7,
    eps_rel::Real = 1e-7,
    adaptive_rho::Bool = true,
    mu_incr::Real = 3.0,
    mu_decr::Real = 3.0,
    tau_incr::Real = 1.05,
    tau_decr::Real = 0.98,
) = begin
    rhos = _normalize_rhos_strict(
        base_model, coupling;
        rho = rho,
        rho_theta = rho_theta,
        rho_v = rho_v,
        rho_v2 = rho_v2,
        rho_p = rho_p,
        rho_q = rho_q,
        rho_r = rho_r,
        rho_t = rho_t,
    )

    DistributedOPFModel(
        base_model, coupling, sectioning, max_iter,
        rhos.rho_theta, rhos.rho_v, rhos.rho_p, rhos.rho_q, rhos.rho_r, rhos.rho_t,
        eps_abs, eps_rel;
        adaptive_rho = adaptive_rho,
        mu_incr = mu_incr,
        mu_decr = mu_decr,
        tau_incr = tau_incr,
        tau_decr = tau_decr,
    )
end

# ============================================================
# MODEL TYPE CHECKS
# ============================================================

@inline is_distributed_model(::Type{<:OPFModel}) = false
@inline is_distributed_model(::DistributedOPFModel) = true
@inline is_centralized_model(x) = !is_distributed_model(x)

# ============================================================
# shims
# ============================================================

rho(c::DistributedOPFModel) = c.rho_theta

Base.getproperty(c::DistributedOPFModel, s::Symbol) =
    s === :rho ? c.rho_theta :
    s === :rho_v2 ? c.rho_v :
    getfield(c, s)

function Base.setproperty!(c::DistributedOPFModel, s::Symbol, v)
    if s === :rho
        v isa Real || throw(ArgumentError("rho must be Real, got $(typeof(v))"))
        vf = Float64(v)
        vf > 0 || throw(ArgumentError("rho must be > 0, got $vf"))

        rules = _rules_for(c.base_model, c.coupling)
        req = rules.required

        if :rho_theta in req
            c.rho_theta = vf
        end
        if :rho_v in req
            c.rho_v = vf
        end
        if :rho_v2 in req
            c.rho_v = vf
        end
        if :rho_p in req
            c.rho_p = vf
        end
        if :rho_q in req
            c.rho_q = vf
        end
        if :rho_r in req
            c.rho_r = vf
        end
        if :rho_t in req
            c.rho_t = vf
        end

        return vf
    elseif s === :rho_v2
        c.rho_v = Float64(v)
        return c.rho_v
    else
        Base.setfield!(c, s, v)
        return v
    end
end

Base.propertynames(c::DistributedOPFModel, private::Bool = false) =
    (propertynames(typeof(c), private)..., :rho, :rho_v2)

# ============================================================
# Diagnostics
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

Base.@kwdef struct LirionRunResult
    result::Any = nothing
    model::Any = nothing
    system::Any = nothing
    diagnostics::AbstractRunDiagnostics = NoDiagnostics()
end
