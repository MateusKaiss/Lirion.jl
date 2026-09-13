# ============================================================
# ALADIN Algorithm 3 — globalization branch (c): dual-only step
# ============================================================

@inline function _mark_local_models_probed!(ws::AbstractALADINWorkspace)
    ws.hist.meta[:local_models_probed] = true
    return nothing
end

function _vrho_resolve_sum!(ws::AbstractALADINWorkspace)
    _mark_local_models_probed!(ws)
    parts = ws.partitions.parts
    vals = zeros(Float64, length(parts))
    Threads.@threads for i in eachindex(parts)
        part = parts[i]
        update_objective!(part, ws)
        JuMP.optimize!(part.local_model)
        _require_optimal!(part.local_model; iter = ws.iter, zone_id = part.id)
        _note_almost_zone_solve!(part.local_model)
        vals[i] = JuMP.objective_value(part.local_model)
    end
    return sum(vals)
end

function _golden_max(f, a::Float64, b::Float64; iters::Int = 6)
    invphi = 0.6180339887498949            # 1/φ
    c = b - invphi * (b - a)
    d = a + invphi * (b - a)
    fc = f(c)
    fd = f(d)
    @inbounds for _ in 1:iters
        if fc < fd
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd = f(d)
        else
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc = f(c)
        end
    end
    return fc > fd ? c : d
end

function _branch_c_dual_linesearch!(ws::AbstractALADINWorkspace)
    states = collect(values(ws.vars))
    saved = [copy(s.lambda) for s in states]
    dir = [s.nu .- s.lambda for s in states]

    function setlambda!(alpha)
        @inbounds for (s, s0, d) in zip(states, saved, dir)
            @. s.lambda = s0 + alpha * d
        end
        return nothing
    end

    vrho(alpha) = (setlambda!(alpha); _vrho_resolve_sum!(ws))

    alpha3 = _golden_max(vrho, 0.0, 1.0; iters = 6)
    @inbounds for (s, s0) in zip(states, saved)
        copyto!(s.lambda, s0)
    end

    ws.alpha1 = 0.0
    ws.alpha2 = 0.0
    ws.alpha3 = alpha3
    return nothing
end

# ============================================================
# ALADIN Algorithm 3 — branch (c) for :alg3_original (Wolfe ascent search)
# ============================================================

function _vrho_resolve!(ws::AbstractALADINWorkspace)
    _mark_local_models_probed!(ws)
    parts = ws.partitions.parts
    vals = zeros(Float64, length(parts))
    Threads.@threads for i in eachindex(parts)
        part = parts[i]
        update_objective!(part, ws)
        JuMP.optimize!(part.local_model)
        _require_optimal!(part.local_model; iter = ws.iter, zone_id = part.id)
        _note_almost_zone_solve!(part.local_model)
        vals[i] = JuMP.objective_value(part.local_model)
        extract_solve_values!(ws.zones[i], part.local_model)
    end
    return sum(vals)
end

function _vrho_cached_sum(ws::AbstractALADINWorkspace)
    s = 0.0
    for part in ws.partitions.parts
        s += JuMP.objective_value(part.local_model)
    end
    return s
end

_dual_slope(ws::AbstractALADINWorkspace, dirs) = error(
    "globalization=:alg3_original needs `_dual_slope(::$(typeof(ws)), dirs)` " *
    "(Danskin slope of V_ρ along λ_QP − λ). Implemented for AC, DC and SOCWR " *
    "in aladin/{ac,dc,socwr}/line_search.jl.",
)

_dual_residual(ws::AbstractALADINWorkspace) = error(
    "globalization=:alg3_original needs `_dual_residual(::$(typeof(ws)))` " *
    "(Danskin gradient Σ Aᵢyᵢ − b, keyed like ws.vars). Implemented for AC, " *
    "DC and SOCWR in aladin/{ac,dc,socwr}/line_search.jl.",
)

const _ALG3_FALLBACK_GATE = Ref(1e-8)

# Outcome relabelling for a step found along the Danskin gradient.
@inline function _branch_c_grad_status(status::Symbol)
    status === :full_step && return :full_step_grad
    status === :wolfe && return :wolfe_grad
    status === :armijo_only && return :armijo_grad
    status === :ascent_only && return :ascent_grad
    return status                                   # :none passes through
end

function _vrho_resolve_trial!(ws::AbstractALADINWorkspace)
    try
        return _vrho_resolve!(ws)
    catch err
        ws.hist.meta[:branch_c_trial_failures] =
            Int(get(ws.hist.meta, :branch_c_trial_failures, 0)) + 1
        @warn "ALADIN :alg3_original branch (c) — a zone NLP failed inside a " *
              "trial V_ρ evaluation (iter $(ws.iter)); the trial step is " *
              "rejected (φ = −∞). Counted in hist.meta[:branch_c_trial_failures]." exception = err maxlog = 1
        return -Inf
    end
end

function _branch_c_dual_linesearch_wolfe!(
    ws::AbstractALADINWorkspace;
    c1::Float64 = 1e-4,
    c2::Float64 = 0.1,
    max_evals::Int = 10,
    alpha_min::Float64 = 1e-4,
    phi_tol_rel::Float64 = 1e-8,
)
    states = collect(values(ws.vars))
    saved = [copy(s.lambda) for s in states]
    dirs = map(s -> s.nu .- s.lambda, ws.vars)     # NamedTuple, keys of ws.vars
    dirs_v = collect(values(dirs))
    saved_y = [copy(z.y) for z in ws.zones]

    # λ(α) = λ_saved + α·dir — absolute from the saved λ for ANY direction.
    function setlambda!(alpha, dir_v)
        @inbounds for (s, s0, d) in zip(states, saved, dir_v)
            @. s.lambda = s0 + alpha * d
        end
        return nothing
    end
    function restore!()
        @inbounds for (s, s0) in zip(states, saved)
            copyto!(s.lambda, s0)
        end
        @inbounds for (z, y0) in zip(ws.zones, saved_y)
            copyto!(z.y, y0)
        end
        return nothing
    end

    n_evals0 = 0
    if ws.iter == 1
        setlambda!(0.0, dirs_v)
        phi0 = _vrho_resolve!(ws)
        n_evals0 += 1
    else
        phi0 = _vrho_cached_sum(ws)
    end
    slope0_d = _dual_slope(ws, dirs)

    grad = _dual_residual(ws)                      # NamedTuple, keys of ws.vars
    mu = ws.mu[]
    gdir = map(g -> mu .* g, grad)                 # d_g = μ g, keys of ws.vars
    gdir_v = collect(values(gdir))
    slope0_g = _dual_slope(ws, gdir)               # = μ‖g‖² ≥ 0
    g_inf = 0.0                                    # ‖g‖∞ of the RAW g (gate)
    for g in values(grad)
        @inbounds for v in g
            av = abs(v)
            av > g_inf && (g_inf = av)
        end
    end

    tau = phi_tol_rel * max(1.0, abs(phi0))

    function run_search(dir, dir_v, slope0)
        n_evals = n_evals0
        best_alpha = 0.0
        best_phi = -Inf
        function evalphi(alpha)
            setlambda!(alpha, dir_v)
            phi = _vrho_resolve_trial!(ws)         # −∞ on a failed zone solve
            slope = isfinite(phi) ? _dual_slope(ws, dir) : NaN
            n_evals += 1
            if alpha > 0.0 && phi > best_phi
                best_phi = phi
                best_alpha = alpha
            end
            return phi, slope
        end
        armijo(a, p) = p >= phi0 + c1 * a * slope0 - tau

        alpha = 0.0
        status = :none

        phi1, slope1 = evalphi(1.0)
        if isfinite(phi1) && armijo(1.0, phi1) && slope1 >= -c2 * slope0
            alpha = 1.0
            status = :full_step
        else
            if isfinite(phi1) && armijo(1.0, phi1)
                # overshoot: φ(1) is good but φ' strongly negative — the
                # maximiser is interior; bracket runs from 1 back toward 0.
                lo, phi_lo, slope_lo = 1.0, phi1, slope1
                hi, phi_hi, slope_hi = 0.0, phi0, slope0
            else
                lo, phi_lo, slope_lo = 0.0, phi0, slope0
                hi, phi_hi, slope_hi = 1.0, phi1, slope1
            end
            while n_evals < max_evals
                w = abs(hi - lo)
                w < alpha_min && break
                a = if isfinite(slope_hi) && slope_lo * slope_hi < 0.0
                    (lo * slope_hi - hi * slope_lo) / (slope_hi - slope_lo)
                else
                    0.5 * (lo + hi)
                end
                a = clamp(a, min(lo, hi) + 0.1 * w, max(lo, hi) - 0.1 * w)
                p, s = evalphi(a)
                if !(isfinite(p) && armijo(a, p)) || p <= phi_lo
                    hi, phi_hi, slope_hi = a, p, s
                else
                    if abs(s) <= c2 * slope0
                        alpha = a
                        status = :wolfe
                        break
                    end
                    if s * (hi - lo) <= 0.0
                        hi, phi_hi, slope_hi = lo, phi_lo, slope_lo
                    end
                    lo, phi_lo, slope_lo = a, p, s
                end
            end
            if status === :none && lo > 0.0
                alpha = lo
                status = :armijo_only
            end
        end

        if status === :none && best_alpha > 0.0 && best_phi > phi0
            alpha = best_alpha
            status = :ascent_only
        end
        return alpha, status
    end

    function run_backtrack(dir, dir_v, slope0; beta::Float64 = 0.1)
        n_evals = n_evals0
        best_alpha = 0.0
        best_phi = -Inf
        armijo(a, p) = p >= phi0 + c1 * a * slope0 - tau

        a = 1.0
        while n_evals < max_evals
            setlambda!(a, dir_v)
            p = _vrho_resolve_trial!(ws)           # −∞ on a failed zone solve
            n_evals += 1
            if isfinite(p) && p > best_phi
                best_phi = p
                best_alpha = a
            end
            if isfinite(p) && armijo(a, p) && p > phi0
                return a, (a == 1.0 ? :full_step : :wolfe)
            end
            a *= beta
        end
        if best_alpha > 0.0 && best_phi > phi0
            return best_alpha, :ascent_only
        end
        return 0.0, :none
    end

    alpha = 0.0
    status = :none

    if isfinite(slope0_d) && slope0_d > 0.0
        alpha, status = run_search(dirs, dirs_v, slope0_d)
    end

    grad_step = false
    if status === :none && isfinite(slope0_g) && slope0_g > 0.0 &&
       g_inf > _ALG3_FALLBACK_GATE[]
        alpha, status = run_backtrack(gdir, gdir_v, slope0_g)
        grad_step = status !== :none
        status = _branch_c_grad_status(status)
    end

    if status === :none
        alpha = 0.0
        status = :no_ascent
    end

    restore!()
    if grad_step
        @inbounds for (s, s0, dg) in zip(states, saved, gdir_v)
            @. s.nu = s0 + dg
        end
    end
    ws.alpha1 = 0.0
    ws.alpha2 = 0.0
    ws.alpha3 = alpha
    return status
end

# ============================================================
# Post-loop: put the local models back on the accepted iterate
# ============================================================
function _restore_local_models!(ws::AbstractALADINWorkspace)
    (ws.alpha1 == 0.0 && ws.alpha2 == 0.0) || return :probe   # x moved: no restore point
    rho_key = hasproperty(ws, :rho) ? :rho : :rho_bus
    rho_ref = getproperty(ws, rho_key)::Base.RefValue{Float64}
    rho_now = rho_ref[]
    rho_rec = get(ws.hist.rhos, rho_key, Float64[])
    rho_acc = isempty(rho_rec) ? rho_now : rho_rec[end]
    rho_swap = rho_acc != rho_now

    states = collect(values(ws.vars))
    lam_now = [copy(s.lambda) for s in states]
    @inbounds for s in states
        copyto!(s.lambda, s.lambda_prev)
    end
    if rho_swap
        rho_ref[] = rho_acc
        refresh_rho!(ws)
    end

    outcome = :restored
    try
        parts = ws.partitions.parts
        Threads.@threads for i in eachindex(parts)
            part = parts[i]
            update_objective!(part, ws)
            JuMP.optimize!(part.local_model)
            _require_optimal!(part.local_model; iter = ws.iter, zone_id = part.id)
            _note_almost_zone_solve!(part.local_model)
        end
    catch err
        outcome = :probe
        @warn "ALADIN — could not re-solve the local models at the accepted " *
              "iterate after the final branch-(c) search; per-zone solutions " *
              "are at a PROBE point (hist.meta[:local_models] = :probe)." exception = err
    finally
        @inbounds for (s, l) in zip(states, lam_now)
            copyto!(s.lambda, l)
        end
        if rho_swap
            rho_ref[] = rho_now
            refresh_rho!(ws)
        end
    end
    return outcome
end

function _finalize_local_models!(ws::AbstractALADINWorkspace)
    probed = get(ws.hist.meta, :local_models_probed, false) === true
    outcome = probed ? _restore_local_models!(ws) : :accepted
    ws.hist.meta[:local_models] = outcome
    return outcome
end
