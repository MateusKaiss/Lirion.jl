# ============================================================
# Core loop
# ============================================================

# Certified from iteration 2 on: at iteration 1 r_step_primal ≡ 0 by construction.
@inline function stop_if_converged(ws::AbstractALADINWorkspace, params::ALADINParams)
    ws.iter >= 2 || return false
    return (ws.r_primal < ws.eps_primal) && (ws.r_step_primal < ws.eps_step_primal)
end

@inline _line_search_step!(::AbstractALADINWorkspace) = nothing

# Terminal status (ws.hist.meta[:terminal_status]):
#   :converged        stop_if_converged fired
#   :stalled          max_iter reached with α = (0, 0, 0) (iterate frozen)
#   :iteration_limit  max_iter reached while the iterate was still moving
@inline _aladin_terminal_status(converged::Bool, zero_step::Bool) =
    converged ? :converged : (zero_step ? :stalled : :iteration_limit)

# Early exit after this many consecutive frozen iterations (α = (0,0,0) with
# ρ, μ unchanged); the terminal status stays :stalled.
const _STALL_EXIT_STEPS = Ref(10)

function run_algorithm!(ws::AbstractALADINWorkspace, params::ALADINParams)
    tag = _print_tag(ws)
    alpha1_cfg = ws.alpha1    # cfg.alpha1 (HFD step 5 α₁)
    alpha3_cfg = ws.alpha3
    _reset_solve_counters!()
    converged = false
    zero_step = false        # last globalization returned α = (0, 0, 0)
    n_zero_steps = 0
    frozen_streak = 0        # consecutive iterations with α = (0,0,0) AND ρ, μ unchanged

    for iter in 1:params.max_iter
        ws.iter = iter
        _print_iter_header!(tag, iter)

        local_solve!(ws)
        ws.hist.meta[:local_models_probed] = false   # models = accepted solve
        solve_coordinator_qp!(ws)
        compute_residuals!(ws)
        compute_tolerances!(ws, params.eps_abs, params.eps_rel)
        _record_iter!(ws)
        _print_aladin_residuals!(ws)

        if stop_if_converged(ws, params)
            _print_converged!()
            converged = true
            break
        end

        ws.alpha1 = alpha1_cfg
        ws.alpha3 = alpha3_cfg
        _line_search_step!(ws)
        zero_step = ws.alpha1 == 0.0 && ws.alpha2 == 0.0 && ws.alpha3 == 0.0
        zero_step && (n_zero_steps += 1)
        primal_update!(ws)
        dual_update!(ws)

        rhos_before = _get_rhos(ws)
        mu_before = ws.mu[]
        adapt_rho_mu!(ws, params)
        frozen = zero_step && _get_rhos(ws) == rhos_before && ws.mu[] == mu_before
        frozen_streak = frozen ? frozen_streak + 1 : 0
        if frozen_streak >= _STALL_EXIT_STEPS[] && iter < params.max_iter
            ws.hist.meta[:early_stall_exit] = iter
            break
        end
    end
    # re-solve the local models at the accepted iterate if branch (c) probed them
    _finalize_local_models!(ws)
    _flush_solve_counters!(ws)
    ws.hist.meta[:zero_step_iters] = n_zero_steps
    ws.hist.meta[:terminal_status] = _aladin_terminal_status(converged, zero_step)
    return nothing
end

# ============================================================
# Entry point
# ============================================================

function run_algorithm(
    sys::DistributedSystem,
    config::ALADINConfig;
    optimizer_factory = nothing,
)
    t0_solver = time_ns()
    tag = ALADIN{_variant_type(Val(config.coupling)), config.base_model}()

    params = build_params(config)
    ws = build_workspace(sys, config, optimizer_factory, tag)
    initial_rhos = _get_rhos(ws)
    run_algorithm!(ws, params)

    tag_str = _print_tag(ws)
    final_rhos = _get_rhos(ws)
    rho_lines = [
        "$(k) initial: $(initial_rhos[k]) | $(k) final: $(final_rhos[k])"
        for k in keys(final_rhos)
    ]
    push!(rho_lines, "mu: $(ws.mu[])")

    if _FAST_USE_LOGS
        _print_summary!(
            tag_str, "ALADIN",
            _finalize(ws.partitions, t0_solver, ws.iter);
            rho_lines = rho_lines,
        )
    end

    return ws.partitions, ws, ws.iter
end

# ============================================================
# Runtime manifest (reproducibility)
# ============================================================

# Every module-level Ref/ENV toggle that can change an ALADIN trajectory, read at call time.
function aladin_runtime_manifest()
    return Dict{String, Any}(
        # coordinator-QP Hessian: prox-in-B (ρΣ added to Bᵢ) toggles
        "AC_QP_NO_PROX" => _AC_QP_NO_PROX[],
        "SOCWR_QP_PROX_SCALE" => Float64(_SOCWR_QP_PROX_SCALE[]),
        # SOCWR active-set activity thresholds (ported knob names)
        "SOCWR_ACTIVE_DUAL_TOL" => Float64(_SOCWR_ACTIVE_DUAL_TOL[]),
        "SOCWR_ACTIVE_CONE_TOL" => Float64(_SOCWR_ACTIVE_CONE_TOL[]),
        # HFD modified gradient gᵢ = ∇fᵢ + (Cᵢ* − Cᵢ)ᵀκᵢ (AC)
        "MODIFIED_GRADIENT" => _MODIFIED_GRADIENT[],
        # Algorithm 3 / :alg3_original constants
        "ALG3_PHI_LOCAL_VIOL" => _ALG3_PHI_LOCAL_VIOL[],
        "ALG3_ORIG_SAFETY" => Float64(_ALG3_ORIG_SAFETY[]),
        "ALG3_ORIG_FLOOR" => Float64(_ALG3_ORIG_FLOOR[]),
        "ALG3O_LAMBAR_HYST" => _ALG3O_LAMBAR_HYST,
        "ALG3O_LAMBAR_FACTOR" => _ALG3O_LAMBAR_FACTOR,
        "ALG3O_LAMBAR_MIN" => _ALG3O_LAMBAR_MIN,
        "ALG3_ORIGINAL_LAM_BAR_FLOOR_DC" => _ALG3_ORIGINAL_LAM_BAR_FLOOR,
        # :alg3_original branch-(c) fallback gate (core/branch_c.jl, ‖g‖∞ floor)
        "ALG3_FALLBACK_GATE" => Float64(_ALG3_FALLBACK_GATE[]),
        # frozen-run early exit threshold (consecutive α=(0,0,0) with ρ,μ fixed)
        "STALL_EXIT_STEPS" => Int(_STALL_EXIT_STEPS[]),
        # active-set / regularisation toggles
        "BOUND_DETECT_MODE" => String(_BOUND_DETECT_MODE[]),
        "GAUGE_REG_MODE" => String(_GAUGE_REG[][1]),
        "GAUGE_REG_GAMMA" => Float64(_GAUGE_REG[][2]),
        "LIRION_GAUGEFIX" => get(ENV, "LIRION_GAUGEFIX", "unset"),
        "ACTIVE_DEP_TOL" => Float64(_ACTIVE_DEP_TOL[]),
        "ACTIVE_DUAL_TOL" => Float64(_ACTIVE_DUAL_TOL[]),
        "QP_HESS_REG" => Float64(_QP_HESS_REG[]),
        # AC/SOCWRNL coordinator-Hessian PSD safeguard: shift (default) | flip | clip
        "AC_PSD_MODE" => String(_PSD_MODE[]),
        "SD_PER_ITER" => _SD_PER_ITER[],
        # SOCWRNL exact-Hessian prox-in-B opt-out (socwr/sensitivities.jl)
        "SOCWRNL_QP_NO_PROX" => _SOCWRNL_QP_NO_PROX[],
        # modified-gradient weak-row floor (ac/sensitivities.jl)
        "MODGRAD_DUAL_FLOOR" => _MODGRAD_DUAL_FLOOR,
        # box-bound activity test (core/sensitivities.jl detect_active_bounds!, all models)
        "BOUND_DETECT_TOL" => _BOUND_DETECT_TOL,
        "BOUND_DETECT_DUAL_TOL" => _BOUND_DETECT_DUAL_TOL,
        # DC / conic-SOCWR row activity thresholds (kwarg constants)
        "DC_ACTIVE_DUAL_TOL" => _DC_ACTIVE_DUAL_TOL,
        "DC_ACTIVE_DEP_TOL" => _DC_ACTIVE_DEP_TOL,
        "SOCWR_ACTIVE_DEP_TOL" => _SOCWR_ACTIVE_DEP_TOL,
        "LIRION_HSLLIB" => get(ENV, "LIRION_HSLLIB", "unset"),
        "LIRION_STRICT_SOLVERS" => get(ENV, "LIRION_STRICT_SOLVERS", "unset"),
        "JULIA_NUM_THREADS" => Threads.nthreads(),
    )
end

# Per-run record of what was actually applied (globalization, λ̄/κ̄, terminal
# status, solve counts, effective Σ). `cfg` is optional; when given, its
# configured values are recorded verbatim.
function aladin_run_manifest(ws::AbstractALADINWorkspace, cfg = nothing)
    meta = ws.hist.meta
    lam_bar = hasproperty(ws, :lam_bar) ? Float64(ws.lam_bar) : NaN
    kap_bar = hasproperty(ws, :kap_bar) ? Float64(ws.kap_bar) : lam_bar
    fixed = hasproperty(ws, :lam_bar_fixed) ? Bool(ws.lam_bar_fixed) :
            (hasproperty(ws, :alg3_lam_bar) ? ws.alg3_lam_bar > 0.0 : false)
    glob = cfg isa ALADINConfig ? cfg.globalization :
           (hasproperty(ws, :globalization) ? ws.globalization : :unknown)
    lam_cfg = cfg isa ALADINConfig ? Float64(cfg.alg3_lam_bar) :
              (fixed ? lam_bar : 0.0)
    pmode = hasproperty(ws, :proximal_mode) ? String(ws.proximal_mode) : "primal"
    return Dict{String, Any}(
        "GLOBALIZATION" => String(glob),
        "ALG3_LAM_BAR_CFG" => lam_cfg,
        "LAM_BAR_APPLIED" => lam_bar,
        "KAP_BAR_APPLIED" => kap_bar,
        "LAM_BAR_FIXED" => fixed,
        "LAM_BAR_VIOLATIONS" => _lam_bar_violations(meta),
        "TERMINAL_STATUS" => String(get(meta, :terminal_status, :unknown)),
        "ZERO_STEP_ITERS" => Int(get(meta, :zero_step_iters, 0)),
        # iteration at which a frozen run was cut (0 = ran to stop/max_iter)
        "EARLY_STALL_EXIT" => Int(get(meta, :early_stall_exit, 0)),
        # δ-PSD shift of the AC coordinator Hessian (hist.meta[:psd_delta]); NaN/0 when not recorded
        "PSD_DELTA_FRAC" => _psd_delta_frac(meta),
        "PSD_DELTA_MAX" => _psd_delta_max(meta),
        "PSD_DELTA_ITERS" => _psd_delta_iters(meta),
        # :alg3_original branch-(c) outcome counts: visits, gradient fallback steps, :no_ascent freezes
        "N_BRANCH_C" => _alg3o_branch_c_counts(meta)[1],
        "N_GRAD_STEPS" => _alg3o_branch_c_counts(meta)[2],
        "N_NO_ASCENT" => _alg3o_branch_c_counts(meta)[3],
        # branch-(c) visits that returned on eval-budget exhaustion
        "N_BUDGET_C" => _alg3o_branch_c_counts(meta)[4],
        # row-activity / LICQ thresholds used by this workspace's active set
        "ACTIVE_ROW_DUAL_TOL_EFFECTIVE" => _active_row_tols(ws)[1],
        "ACTIVE_DEP_TOL_EFFECTIVE" => _active_row_tols(ws)[2],
        # κ̄ (= λ̄) majorisation violations against the local duals
        "KAP_BAR_VIOLATIONS" => _kap_bar_violations(meta),
        "ALMOST_ZONE_SOLVES" => sum(get(meta, :almost_zone_solves, Int[]); init = 0),
        "ALMOST_QP_SOLVES" => sum(get(meta, :almost_qp_solves, Int[]); init = 0),
        "LOCAL_MODELS" => String(get(meta, :local_models, :unknown)),
        "PROXIMAL_MODE" => pmode,
        "SIGMA_EFFECTIVE" => _sigma_effective(ws),
    )
end

# Fraction of recorded iterations needing a δ > 0 PSD shift; NaN when nothing recorded.
function _psd_delta_frac(meta::Dict{Symbol, Any})
    v = get(meta, :psd_delta, nothing)
    (v isa AbstractVector && !isempty(v)) || return NaN
    return count(>(0.0), v) / length(v)
end
function _psd_delta_max(meta::Dict{Symbol, Any})
    v = get(meta, :psd_delta, nothing)
    (v isa AbstractVector && !isempty(v)) || return NaN
    return Float64(maximum(v))
end
function _psd_delta_iters(meta::Dict{Symbol, Any})
    v = get(meta, :psd_delta, nothing)
    (v isa AbstractVector && !isempty(v)) || return 0
    return count(>(0.0), v)
end

# (n_branch_c, n_grad_steps, n_no_ascent, n_budget) of :alg3_original; zeros otherwise.
function _alg3o_branch_c_counts(meta::Dict{Symbol, Any})
    rec = get(meta, :alg3_original, nothing)
    if rec isa AbstractVector
        n_c = 0; n_g = 0; n_na = 0; n_b = 0
        for e in rec
            e isa Tuple && length(e) >= 4 || continue
            e[1] === :c || continue
            n_c += 1
            st = e[4]
            st isa Symbol || continue
            endswith(String(st), "_grad") && (n_g += 1)
            st === :no_ascent && (n_na += 1)
            (st === :armijo_only || st === :ascent_only || st === :armijo_grad || st === :ascent_grad) && (n_b += 1)
        end
        return (n_c, n_g, n_na, n_b)
    elseif rec isa Dict
        g(k) = Int(get(rec, k, 0))
        return (g(:c), g(:c_full_grad) + g(:c_wolfe_grad) + g(:c_budget_grad), g(:c_no_ascent), g(:c_budget) + g(:c_budget_grad))
    end
    return (0, 0, 0, 0)
end

# (row dual_tol, dep_tol) used by this workspace type's active-set builder.
function _active_row_tols(ws::AbstractALADINWorkspace)
    ws isa _DCPhase2Workspace && return (_DC_ACTIVE_DUAL_TOL, _DC_ACTIVE_DEP_TOL)
    ws isa _ACPhase2Workspace && return (_ACTIVE_DUAL_TOL[], _ACTIVE_DEP_TOL[])
    exact = !isempty(ws.zones) && ws.zones[1].exact_hessian
    exact && return (_ACTIVE_DUAL_TOL[], _ACTIVE_DEP_TOL[])
    return (_SOCWR_ACTIVE_DUAL_TOL[], _SOCWR_ACTIVE_DEP_TOL)
end

# max|κ_local| > fixed κ̄ (= λ̄) events.
function _kap_bar_violations(meta::Dict{Symbol, Any})
    n = Int(get(meta, :kap_bar_violations, 0))
    tally = get(meta, :alg3_original, nothing)
    tally isa Dict && (n += Int(get(tally, :kap_bar_violations, 0)))
    return n
end

# |λ|∞ > fixed λ̄ events (:alg3_original with cfg.alg3_lam_bar > 0).
function _lam_bar_violations(meta::Dict{Symbol, Any})
    n = Int(get(meta, :lam_bar_violations, 0))
    tally = get(meta, :alg3_original, nothing)
    tally isa Dict && (n += Int(get(tally, :lam_bar_violations, 0)))
    return n
end

# Per-class Σ weights the run actually used (every `sigma_*` field of the workspace).
function _sigma_effective(ws::AbstractALADINWorkspace)
    out = Dict{String, Float64}()
    for f in fieldnames(typeof(ws))
        startswith(String(f), "sigma_") || continue
        out[String(f)] = Float64(getfield(ws, f))
    end
    return out
end
