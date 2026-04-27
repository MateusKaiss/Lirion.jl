# ============================================================
# OPF BASE MODELS
# ============================================================

abstract type OPFModel end

struct DC <: OPFModel end
struct AC <: OPFModel end     # AC polar (V, θ)
struct ACR <: OPFModel end    # AC rectangular (vr, vi)
struct IVR <: OPFModel end    # Current-Voltage rectangular (vr, vi, cr, ci)
struct LPAC <: OPFModel end   # Linear Power AC (Coffrin & Van Hentenryck 2014)
struct SOCWR <: OPFModel end  # Jabr SOCP (W-rectangular)
struct QC <: OPFModel end     # Quadratic-Convex relaxation (Hijazi-Coffrin-VH 2017)
struct SDPWRM <: OPFModel end # Bai-Wei-Fujisawa SDP (W-rectangular matrix)
struct SOCBF <: OPFModel end  # Branch Flow SOCP
struct BFA <: OPFModel end    # Linear Branch Flow Approx (LinDistFlow)
struct NFA <: OPFModel end    # Network Flow Approximation (LP, active-only)

# ============================================================
# ADMM VARIANTS
# ============================================================

abstract type ADMMVariant end

struct State <: ADMMVariant end
struct Power <: ADMMVariant end
struct Hybrid <: ADMMVariant end

# ============================================================
# DISTRIBUTED ALGORITHM HIERARCHY
# ============================================================
abstract type DistributedAlgorithm end

struct ADMM{V <: ADMMVariant, M <: OPFModel} <: DistributedAlgorithm end
struct ALADIN{V <: ADMMVariant, M <: OPFModel} <: DistributedAlgorithm end

# ============================================================
# SECTIONING / DECOMPOSITION STRATEGIES
# ============================================================

abstract type Sectioning end
struct TieLine <: Sectioning end
struct DummyBus <: Sectioning end
struct SplitBus <: Sectioning end
struct Metis <: Sectioning end

# ============================================================
# Algorithm Specifications
# ============================================================

abstract type AbstractAlgorithm end
abstract type AbstractDistributedAlgorithm <: AbstractAlgorithm end

struct CentralizedAlgorithm <: AbstractAlgorithm end

struct ADMMAlgorithm <: AbstractDistributedAlgorithm
    coupling::Symbol
    sectioning::Sectioning
    max_iter::Int
    kwargs::NamedTuple
end

struct ALADINAlgorithm <: AbstractDistributedAlgorithm
    coupling::Symbol
    sectioning::Sectioning
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
    # SOCWR
    # ----------------------------
    (SOCWR, :state) => (
        allowed = Set([:rho_v2, :rho_r, :rho_t]),
        required = Set([:rho_v2, :rho_r, :rho_t]),
    ),
    (SOCWR, :power) => (
        allowed = Set([:rho_p, :rho_q]),
        required = Set([:rho_p, :rho_q]),
    ),
    (SOCWR, :hybrid) => (
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
    rho::Union{Nothing, Real} = nothing,
    rho_theta::Union{Nothing, Real} = nothing,
    rho_v::Union{Nothing, Real} = nothing,
    rho_v2::Union{Nothing, Real} = nothing,
    rho_p::Union{Nothing, Real} = nothing,
    rho_q::Union{Nothing, Real} = nothing,
    rho_r::Union{Nothing, Real} = nothing,
    rho_t::Union{Nothing, Real} = nothing,
)
    rules = _rules_for(base_model, coupling)
    allowed = rules.allowed
    required = rules.required

    if base_model === SOCWR
        rho_theta === nothing ||
            error("SOCWR does not accept rho_theta. Use rho_v2, rho_r, rho_t (or rho=...).")
        rho_v === nothing ||
            error("SOCWR does not accept rho_v. Use rho_v2 (or rho=...).")
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

    rho_generic = rho === nothing ? nothing : _asfloat(rho, :rho)

    if rho_generic !== nothing
        (:rho_theta in required) && (rho_theta = something(rho_theta, rho_generic))
        (:rho_v in required) && (rho_v = something(rho_v, rho_generic))
        (:rho_v2 in required) && (rho_v2 = something(rho_v2, rho_generic))
        (:rho_p in required) && (rho_p = something(rho_p, rho_generic))
        (:rho_q in required) && (rho_q = something(rho_q, rho_generic))
        (:rho_r in required) && (rho_r = something(rho_r, rho_generic))
        (:rho_t in required) && (rho_t = something(rho_t, rho_generic))
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
    anchor_val = _get_required_value(req_sorted[1])

    rho_theta_val =
        (base_model === SOCWR) ? anchor_val :
        ((:rho_theta in required) ? _asfloat(rho_theta, :rho_theta) : anchor_val)

    rho_v_val = if base_model === AC
        (:rho_v in required) ? _asfloat(rho_v, :rho_v) : anchor_val
    elseif base_model === SOCWR
        (:rho_v2 in required) ? _asfloat(rho_v2, :rho_v2) : anchor_val
    else
        anchor_val
    end

    rho_p_val = (:rho_p in required) ? _asfloat(rho_p, :rho_p) : anchor_val
    rho_q_val = (:rho_q in required) ? _asfloat(rho_q, :rho_q) : anchor_val
    rho_r_val = (:rho_r in required) ? _asfloat(rho_r, :rho_r) : anchor_val
    rho_t_val = (:rho_t in required) ? _asfloat(rho_t, :rho_t) : anchor_val

    return (
        rho_theta = rho_theta_val,
        rho_v = rho_v_val,
        rho_p = rho_p_val,
        rho_q = rho_q_val,
        rho_r = rho_r_val,
        rho_t = rho_t_val,
    )
end

# ============================================================
# DISTRIBUTED OPF MODEL CONFIGURATION
# ============================================================

abstract type AbstractDistributedConfig end

mutable struct ADMMConfig <: AbstractDistributedConfig
    base_model::Type{<:OPFModel}
    coupling::Symbol
    sectioning::Sectioning
    max_iter::Int

    rho_theta::Float64
    rho_v::Float64
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
    adapt_warmup::Int
    rho_floor::Float64
    adapt_freeze::Float64

    function ADMMConfig(
        base_model::Type{<:OPFModel},
        coupling::Symbol,
        sectioning::Sectioning,
        max_iter::Integer,
        rho_theta::Real,
        rho_v::Real,
        rho_p::Real,
        rho_q::Real,
        rho_r::Real,
        rho_t::Real,
        eps_abs::Real,
        eps_rel::Real;
        adaptive_rho::Bool = false,
        mu_incr::Real = 30.0,
        mu_decr::Real = 10.0,
        tau_incr::Real = 1.2,
        tau_decr::Real = 0.5,
        adapt_warmup::Integer = 5,
        rho_floor::Real = 0.0,
        adapt_freeze::Real = 100.0,
    )
        base_model in (DC, AC, SOCWR) ||
            throw(ArgumentError("Invalid base_model=$base_model. Allowed: DC, AC, SOCWR"))

        coupling in (:state, :power, :hybrid) ||
            throw(
                ArgumentError(
                    "Invalid coupling=$coupling. Allowed: :state, :power, :hybrid",
                ),
            )

        max_iter > 0 || throw(ArgumentError("max_iter must be > 0"))

        rho_theta > 0 || throw(ArgumentError("rho_theta must be > 0"))
        rho_v > 0 || throw(ArgumentError("rho_v must be > 0"))
        rho_p > 0 || throw(ArgumentError("rho_p must be > 0"))
        rho_q > 0 || throw(ArgumentError("rho_q must be > 0"))
        rho_r > 0 || throw(ArgumentError("rho_r must be > 0"))
        rho_t > 0 || throw(ArgumentError("rho_t must be > 0"))

        eps_abs >= 0 || throw(ArgumentError("eps_abs must be ≥ 0"))
        eps_rel >= 0 || throw(ArgumentError("eps_rel must be ≥ 0"))

        mu_incr > 0 || throw(ArgumentError("mu_incr must be > 0"))
        mu_decr > 0 || throw(ArgumentError("mu_decr must be > 0"))
        tau_incr > 0 || throw(ArgumentError("tau_incr must be > 0"))
        tau_decr > 0 || throw(ArgumentError("tau_decr must be > 0"))
        adapt_warmup >= 0 || throw(ArgumentError("adapt_warmup must be ≥ 0"))
        rho_floor >= 0 || throw(ArgumentError("rho_floor must be ≥ 0"))
        adapt_freeze >= 0 || throw(ArgumentError("adapt_freeze must be ≥ 0"))

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
            Int(adapt_warmup),
            Float64(rho_floor),
            Float64(adapt_freeze),
        )
    end
end

# ============================================================
# DISTRIBUTED OPF MODEL CONVENIENCE CONSTRUCTORS
# ============================================================

ADMMConfig(
    base_model::Type{<:OPFModel},
    coupling::Symbol,
    sectioning::Union{Sectioning, Symbol},
    max_iter::Int,
    rho::Real,
    eps_abs::Real,
    eps_rel::Real;
    adaptive_rho::Bool = false,
    mu_incr::Real = 30.0,
    mu_decr::Real = 10.0,
    tau_incr::Real = 1.2,
    tau_decr::Real = 0.5,
    adapt_warmup::Integer = 5,
    rho_floor::Real = 0.0,
    adapt_freeze::Real = 100.0,
) = begin
    sec = _sectioning(sectioning)
    rhos = _normalize_rhos_strict(base_model, coupling; rho = rho)

    ADMMConfig(
        base_model,
        coupling,
        sec,
        max_iter,
        rhos.rho_theta,
        rhos.rho_v,
        rhos.rho_p,
        rhos.rho_q,
        rhos.rho_r,
        rhos.rho_t,
        eps_abs,
        eps_rel;
        adaptive_rho = adaptive_rho,
        mu_incr = mu_incr,
        mu_decr = mu_decr,
        tau_incr = tau_incr,
        tau_decr = tau_decr,
        adapt_warmup = adapt_warmup,
        rho_floor = rho_floor,
        adapt_freeze = adapt_freeze,
    )
end

ADMMConfig(
    base_model::Type{<:OPFModel},
    coupling::Symbol,
    sectioning::Union{Sectioning, Symbol},
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
    adaptive_rho::Bool = false,
    mu_incr::Real = 30.0,
    mu_decr::Real = 10.0,
    tau_incr::Real = 1.2,
    tau_decr::Real = 0.5,
    adapt_warmup::Integer = 5,
    rho_floor::Real = 0.0,
    adapt_freeze::Real = 100.0,
) = begin
    sec = _sectioning(sectioning)

    rhos = _normalize_rhos_strict(
        base_model,
        coupling;
        rho = rho,
        rho_theta = rho_theta,
        rho_v = rho_v,
        rho_v2 = rho_v2,
        rho_p = rho_p,
        rho_q = rho_q,
        rho_r = rho_r,
        rho_t = rho_t,
    )

    ADMMConfig(
        base_model,
        coupling,
        sec,
        max_iter,
        rhos.rho_theta,
        rhos.rho_v,
        rhos.rho_p,
        rhos.rho_q,
        rhos.rho_r,
        rhos.rho_t,
        eps_abs,
        eps_rel;
        adaptive_rho = adaptive_rho,
        mu_incr = mu_incr,
        mu_decr = mu_decr,
        tau_incr = tau_incr,
        tau_decr = tau_decr,
        adapt_warmup = adapt_warmup,
        rho_floor = rho_floor,
        adapt_freeze = adapt_freeze,
    )
end

# ============================================================
# ALADIN CONFIGURATION
# ============================================================

mutable struct ALADINConfig <: AbstractDistributedConfig
    base_model::Type{<:OPFModel}
    coupling::Symbol
    sectioning::Sectioning
    max_iter::Int

    rho::Float64
    mu::Float64

    alpha1::Float64
    alpha2::Float64
    alpha3::Float64

    sigma_theta_shared::Float64
    sigma_pg::Float64
    sigma_other::Float64
    sigma_pflow::Float64
    sigma_qflow::Float64
    eps_abs::Float64
    eps_rel::Float64
    active_set_tol::Float64
    adaptive::Bool
    rho_growth::Float64
    mu_growth::Float64
    rho_upper::Float64
    mu_upper::Float64
    proximal_mode::Symbol

    function ALADINConfig(
        base_model::Type{<:OPFModel},
        coupling::Symbol,
        sectioning::Sectioning,
        max_iter::Integer;
        rho::Real = 1e3,
        mu::Real = 1e4,
        alpha1::Real = 1.0,
        alpha2::Real = 1.0,
        alpha3::Real = 1.0,
        sigma_theta_shared::Real = 10.0,
        sigma_pg::Real = 1.0,
        sigma_other::Real = 1.0,
        sigma_pflow::Real = 10.0,
        sigma_qflow::Real = 10.0,
        eps_abs::Real = 1e-4,
        eps_rel::Real = 1e-4,
        active_set_tol::Real = 1e-6,
        adaptive::Bool = false,
        rho_growth::Real = 1.3,
        mu_growth::Real = 2.0,
        rho_upper::Real = 1e6,
        mu_upper::Real = 1e6,
        proximal_mode::Symbol = :primal,
    )
        base_model in (DC, AC, SOCWR) ||
            throw(
                ArgumentError(
                    "ALADIN is currently implemented for DC, AC, and SOCWR. Got: $(base_model)",
                ),
            )
        coupling in (:state, :power, :hybrid) ||
            throw(
                ArgumentError(
                    "Invalid coupling=$(coupling). Allowed: :state, :power, :hybrid",
                ),
            )

        base_model === SOCWR && sectioning isa DummyBus &&
            throw(
                ArgumentError(
                    "ALADIN is not supported for SOCWR + DummyBus sectioning " *
                    "(coordinator QP converges to a relaxed solution). " *
                    "Use ADMM for this combination, or switch sectioning to " *
                    ":tie_line / :split_bus for ALADIN on SOCWR.",
                ),
            )
        max_iter > 0 || throw(ArgumentError("max_iter must be > 0"))
        rho > 0 || throw(ArgumentError("rho must be > 0"))
        mu > 0 || throw(ArgumentError("mu must be > 0"))
        alpha1 > 0 || throw(ArgumentError("alpha1 must be > 0"))
        alpha2 > 0 || throw(ArgumentError("alpha2 must be > 0"))
        alpha3 > 0 || throw(ArgumentError("alpha3 must be > 0"))
        eps_abs >= 0 || throw(ArgumentError("eps_abs must be ≥ 0"))
        eps_rel >= 0 || throw(ArgumentError("eps_rel must be ≥ 0"))
        active_set_tol > 0 || throw(ArgumentError("active_set_tol must be > 0"))
        rho_growth >= 1.0 ||
            throw(ArgumentError("rho_growth must be ≥ 1 (monotonic growth)"))
        mu_growth >= 1.0 ||
            throw(ArgumentError("mu_growth must be ≥ 1 (monotonic growth)"))
        rho_upper >= rho ||
            throw(ArgumentError("rho_upper must be ≥ rho"))
        mu_upper >= mu ||
            throw(ArgumentError("mu_upper must be ≥ mu"))
        sigma_pflow > 0 || throw(ArgumentError("sigma_pflow must be > 0"))
        sigma_qflow > 0 || throw(ArgumentError("sigma_qflow must be > 0"))
        proximal_mode in (:primal, :pflow) ||
            throw(
                ArgumentError(
                    "Invalid proximal_mode=$proximal_mode. Allowed: :primal, :pflow",
                ),
            )
        proximal_mode === :pflow && coupling !== :power &&
            throw(
                ArgumentError(
                    "proximal_mode=:pflow is only supported with coupling=:power " *
                    "(got coupling=$(coupling)). Use coupling=:power or " *
                    "proximal_mode=:primal.",
                ),
            )

        return new(
            base_model,
            coupling,
            sectioning,
            Int(max_iter),
            Float64(rho),
            Float64(mu),
            Float64(alpha1),
            Float64(alpha2),
            Float64(alpha3),
            Float64(sigma_theta_shared),
            Float64(sigma_pg),
            Float64(sigma_other),
            Float64(sigma_pflow),
            Float64(sigma_qflow),
            Float64(eps_abs),
            Float64(eps_rel),
            Float64(active_set_tol),
            adaptive,
            Float64(rho_growth),
            Float64(mu_growth),
            Float64(rho_upper),
            Float64(mu_upper),
            proximal_mode,
        )
    end
end

# ============================================================
# MODEL TYPE CHECKS
# ============================================================

@inline is_distributed_model(::Type{<:OPFModel}) = false
@inline is_distributed_model(::AbstractDistributedConfig) = true
@inline is_centralized_model(x) = !is_distributed_model(x)

# ============================================================
# Shims
# ============================================================

rho(c::ADMMConfig) = c.rho_theta

Base.getproperty(c::ADMMConfig, s::Symbol) =
    s === :rho ? c.rho_theta :
    s === :rho_v2 ? c.rho_v :
    getfield(c, s)

function Base.setproperty!(c::ADMMConfig, s::Symbol, v)
    if s === :rho
        v isa Real || throw(ArgumentError("rho must be Real, got $(typeof(v))"))
        vf = Float64(v)
        vf > 0 || throw(ArgumentError("rho must be > 0, got $vf"))

        rules = _rules_for(c.base_model, c.coupling)
        req = rules.required

        (:rho_theta in req) && (c.rho_theta = vf)
        (:rho_v in req) && (c.rho_v = vf)
        (:rho_v2 in req) && (c.rho_v = vf)
        (:rho_p in req) && (c.rho_p = vf)
        (:rho_q in req) && (c.rho_q = vf)
        (:rho_r in req) && (c.rho_r = vf)
        (:rho_t in req) && (c.rho_t = vf)

        return vf
    elseif s === :rho_v2
        v isa Real || throw(ArgumentError("rho_v2 must be Real, got $(typeof(v))"))
        vf = Float64(v)
        vf > 0 || throw(ArgumentError("rho_v2 must be > 0, got $vf"))
        c.rho_v = vf
        return c.rho_v
    else
        Base.setfield!(c, s, v)
        return v
    end
end

Base.propertynames(c::ADMMConfig, private::Bool = false) =
    (fieldnames(typeof(c))..., :rho, :rho_v2)

# ============================================================
# Diagnostics
# ============================================================

abstract type AbstractRunDiagnostics end
struct NoDiagnostics <: AbstractRunDiagnostics end

# Root of the distributed-workspace hierarchy. Both ADMM and ALADIN
# workspaces descend from this type so that generic code can dispatch
# on a single "distributed workspace" concept without having to list
# the two families separately.
abstract type AbstractDistributedWorkspace end

abstract type AbstractADMMWorkspace <: AbstractDistributedWorkspace end

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
# Inputs Format
# ============================================================

abstract type InputFormat end
struct MatpowerM <: InputFormat end
struct CsvFolder <: InputFormat end
