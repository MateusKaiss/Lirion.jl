# ============================================================
# ADMM
# ============================================================

@inline function _push_rho!(hist::ADMMHistory, key::Symbol, val::Real)
    v = get!(hist.rhos, key, Float64[])
    push!(v, Float64(val))
    return nothing
end

@inline function _record_iter!(ws::AbstractADMMWorkspace)
    push!(ws.hist.r_primal, ws.r_primal)
    push!(ws.hist.r_dual, ws.r_dual)
    push!(ws.hist.eps_primal, ws.eps_primal)
    push!(ws.hist.eps_dual, ws.eps_dual)
    for (k, v) in pairs(_get_rhos(ws))
        _push_rho!(ws.hist, k, v)
    end
    return nothing
end

# ============================================================
# ALADIN
# ============================================================

const _PSD_DELTA_LAST = Base.RefValue{Float64}(0.0)
const _PSD_DELTA_LOCK = ReentrantLock()
const _SD_PER_ITER = Ref(false)

const _ALMOST_LOCK = ReentrantLock()
const _ALMOST_ZONE_SOLVES = Ref(0)
const _ALMOST_QP_SOLVES = Ref(0)

@inline _is_almost_status(ts) =
    ts == MOI.ALMOST_OPTIMAL || ts == MOI.ALMOST_LOCALLY_SOLVED

function _note_almost_zone_solve!(model::JuMP.Model)
    _is_almost_status(JuMP.termination_status(model)) || return nothing
    lock(_ALMOST_LOCK) do
        _ALMOST_ZONE_SOLVES[] += 1
    end
    return nothing
end

function _note_almost_qp_solve!(model::JuMP.Model)
    _is_almost_status(JuMP.termination_status(model)) || return nothing
    lock(_ALMOST_LOCK) do
        _ALMOST_QP_SOLVES[] += 1
    end
    return nothing
end

function _reset_solve_counters!()
    lock(_ALMOST_LOCK) do
        _ALMOST_ZONE_SOLVES[] = 0
        _ALMOST_QP_SOLVES[] = 0
    end
    return nothing
end

function _take_solve_counters!()
    return lock(_ALMOST_LOCK) do
        z = _ALMOST_ZONE_SOLVES[]
        q = _ALMOST_QP_SOLVES[]
        _ALMOST_ZONE_SOLVES[] = 0
        _ALMOST_QP_SOLVES[] = 0
        (z, q)
    end
end

function _flush_solve_counters!(ws::AbstractALADINWorkspace)
    z, q = _take_solve_counters!()
    zv = get!(ws.hist.meta, :almost_zone_solves, Int[])::Vector{Int}
    qv = get!(ws.hist.meta, :almost_qp_solves, Int[])::Vector{Int}
    isempty(zv) ? push!(zv, z) : (zv[end] += z)
    isempty(qv) ? push!(qv, q) : (qv[end] += q)
    return nothing
end

function _sd_now(ws)::Float64
    lam1 = 0.0; nlam = 0
    for st in values(ws.vars)
        lam1 += sum(abs, st.lambda); nlam += length(st.lambda)
    end
    kap1 = 0.0; m_loc = 0; n_vars = 0
    for zone in ws.zones
        model = zone.local_model
        n_vars += zone.n_vars
        for (F, S) in JuMP.list_of_constraint_types(model)
            F === JuMP.VariableRef && continue
            for c in JuMP.all_constraints(model, F, S)
                m_loc += 1
                d = try; JuMP.dual(c); catch; nothing; end
                d === nothing && continue
                kap1 += d isa Number ? abs(d) : sum(abs, d)
            end
        end
        for v in zone.all_vars
            try
                JuMP.has_lower_bound(v) && (kap1 += abs(JuMP.dual(JuMP.LowerBoundRef(v))))
                JuMP.has_upper_bound(v) && (kap1 += abs(JuMP.dual(JuMP.UpperBoundRef(v))))
            catch
            end
        end
    end
    s_max = 100.0
    return max(s_max, (lam1 + kap1) / max(m_loc + nlam + n_vars, 1)) / s_max
end

Base.@kwdef struct ALADINHistory <: AbstractRunDiagnostics
    r_primal::Vector{Float64} = Float64[]
    r_step_primal::Vector{Float64} = Float64[]
    r_dual::Vector{Float64} = Float64[]
    eps_primal::Vector{Float64} = Float64[]
    eps_step_primal::Vector{Float64} = Float64[]
    eps_dual::Vector{Float64} = Float64[]
    lambda_step::Vector{Float64} = Float64[]
    rhos::Dict{Symbol, Vector{Float64}} = Dict{Symbol, Vector{Float64}}()
    mus::Vector{Float64} = Float64[]
    meta::Dict{Symbol, Any} = Dict{Symbol, Any}()
end

@inline function _push_rho_aladin!(hist::ALADINHistory, key::Symbol, val::Real)
    v = get!(hist.rhos, key, Float64[])
    push!(v, Float64(val))
    return nothing
end

@inline function _record_iter!(ws::AbstractALADINWorkspace)
    push!(ws.hist.r_primal, ws.r_primal)
    push!(ws.hist.r_step_primal, ws.r_step_primal)
    push!(ws.hist.r_dual, ws.r_dual)
    push!(ws.hist.eps_primal, ws.eps_primal)
    push!(ws.hist.eps_step_primal, ws.eps_step_primal)
    push!(ws.hist.eps_dual, ws.eps_dual)
    push!(ws.hist.lambda_step, ws.lambda_step)
    push!(ws.hist.mus, ws.mu[])
    for (k, v) in pairs(_get_rhos(ws))
        _push_rho_aladin!(ws.hist, k, v)
    end
    # δ-PSD applied this iteration (max over zones; take-and-reset)
    delta = lock(_PSD_DELTA_LOCK) do
        d = _PSD_DELTA_LAST[]
        _PSD_DELTA_LAST[] = 0.0
        d
    end
    push!(get!(ws.hist.meta, :psd_delta, Float64[])::Vector{Float64}, delta)
    # ALMOST_*-accepted zone / coordinator-QP solves this iteration (take-and-reset)
    z_almost, q_almost = _take_solve_counters!()
    push!(get!(ws.hist.meta, :almost_zone_solves, Int[])::Vector{Int}, z_almost)
    push!(get!(ws.hist.meta, :almost_qp_solves, Int[])::Vector{Int}, q_almost)
    # s_d per-iteration (opt-in; serial context — safe for JuMP.dual)
    if _SD_PER_ITER[]
        push!(get!(ws.hist.meta, :sd, Float64[])::Vector{Float64}, _sd_now(ws))
    end
    return nothing
end
