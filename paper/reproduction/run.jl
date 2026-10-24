# Re-runs every Table I cell from cases/<case>/gold_params.toml, saves the
# residual history + meta per run, then builds and checks the crossing table.
# julia --project=<repo> paper/reproduction/run.jl [case] [k form sect coup]

# ---- include path ------------------------------------------------------
const _REPO_ROOT = normpath(get(ENV, "LIRION_REPO_ROOT",
    normpath(joinpath(@__DIR__, "..", ".."))))   # @__DIR__/../.. = <repo>
const _LIB_DIR = joinpath(@__DIR__, "scripts", "lib")
isfile(joinpath(_LIB_DIR, "tuning_common.jl")) || error(
    "cannot find tuning_common.jl under\n    $(_LIB_DIR)\n" *
    "run this script from a complete paper/reproduction/ folder.")
include(joinpath(_LIB_DIR, "tuning_common.jl"))

using LinearAlgebra
using Logging
BLAS.set_num_threads(1)

# per-case parameter files: cases/<case>/gold_params.toml
const CASES_DIR = joinpath(@__DIR__, "cases")

# bound-activity detector: default primal_dual here (src default is :primal)
Lirion._BOUND_DETECT_MODE[] = Symbol(get(ENV, "LIRION_BOUND_DETECT", "primal_dual"))

# LIRION_SD_PER_ITER=1: record s_d every iteration (slow)
get(ENV, "LIRION_SD_PER_ITER", "0") == "1" && (Lirion._SD_PER_ITER[] = true)
const OUT_DIR = get(ENV, "LIRION_REPRO_OUT", joinpath(@__DIR__, "results", "runs"))
const TABLE_DIR = joinpath(@__DIR__, "results", "table")
const CENTRAL_CSV = joinpath(@__DIR__, "scripts", "repro_centralized_baseline.csv")
const SMOKE = get(ENV, "LIRION_REPRO_SMOKE", "0") == "1"

# ---- env guards ---------------------------------------------------------
# opt-ins the paper runs kept off; any value aborts
for (var, why) in (
    ("LIRION_GAUGEFIX",       "switches the coordinator-QP gauge regularisation (ALADIN trajectory)"),
    ("LIRION_ZONE_WARMSTART", "enables the opt-in ADMM zone warm-start (ADMM trajectory)"),
    ("LIRION_ZONE_TOL_TIGHT", "tightens the zone solvers off the paper's AC/DC config"),
    ("LIRION_AC_PSD_MODE",    "changes the AC PSD safeguard off the paper's :shift default"),
)
    haskey(ENV, var) && error(
        "repro guard: $(var)=$(ENV[var]) is set. It $(why) and is NOT part of the " *
        "paper configuration — unset it (the paper runs ran with it unset).")
end

# HEAD SHA goes into every meta TOML
function git_sha()
    try
        return strip(read(`git -C $(_REPO_ROOT) rev-parse HEAD`, String))
    catch
        return "unknown"
    end
end
const GIT_COMMIT = git_sha()
# official runs are STRICT — solver fallbacks become errors (tuning_common)
SMOKE || (ENV["LIRION_STRICT_SOLVERS"] = "1")
ensure_hsl!(require = !SMOKE)

_arg(i) = (length(ARGS) >= i && ARGS[i] != "-") ? ARGS[i] : ""
const F_CASE = _arg(1); const F_K = _arg(2); const F_FORM = _arg(3)
const F_SECT = _arg(4); const F_COUP = _arg(5)

# ---- parameter files ----------------------------------------------------
# case filter matches either way ("case14" ↔ "case14_ieee")
function gold_files()
    haskey(ENV, "LIRION_REPRO_GOLD") && return [ENV["LIRION_REPRO_GOLD"]]
    isdir(CASES_DIR) || error("missing $(CASES_DIR)")
    dirs = sort(filter(d -> isfile(joinpath(CASES_DIR, d, "gold_params.toml")),
                       readdir(CASES_DIR)))
    if !isempty(F_CASE)
        hit = filter(d -> occursin(d, F_CASE) || occursin(F_CASE, d), dirs)
        isempty(hit) && error("no cases/<case>/gold_params.toml matches \"$(F_CASE)\" " *
                              "(available: $(join(dirs, ", ")))")
        dirs = hit
    end
    return [joinpath(CASES_DIR, d, "gold_params.toml") for d in dirs]
end

# ---- s_d: Ipopt dual-infeasibility scaling, final iterate --------------
# s_d = max(s_max, (‖λ‖₁+‖κ‖₁)/(m+n)) / s_max
const S_MAX = 100.0

function dual_norms_and_counts(ws)
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
        # bound duals: numerator only (Ipopt's z), not counted in m
        for v in zone.all_vars
            try
                JuMP.has_lower_bound(v) && (kap1 += abs(JuMP.dual(JuMP.LowerBoundRef(v))))
                JuMP.has_upper_bound(v) && (kap1 += abs(JuMP.dual(JuMP.UpperBoundRef(v))))
            catch
            end
        end
    end
    return lam1, kap1, m_loc + nlam, n_vars
end

# ---- ‖Aᵀλ‖₁ (triangle majorant; NaN if a dummy class is populated) -----
const _PM1_CLASSES = (:theta, :v_mag, :v2, :r, :t)

function _flow_row_l1(ws, zone, le, f::AbstractString, which::Symbol)
    if f == "dc"
        return 2.0 * abs(zone.edge_b[le])
    elseif f == "socwr" || f == "socwrnl"
        # socwrnl shares the socwr W-space zone data
        coefs = which === :p ? zone.edge_pflow_coefs[le] : zone.edge_qflow_coefs[le]
        return sum(abs, coefs)
    else # ac
        lin = Lirion._ac_flow_linearisation(zone, le)
        return which === :p ?
            abs(lin.dpf_dvf) + abs(lin.dpf_dvt) + abs(lin.dpf_dtheta_f) + abs(lin.dpf_dtheta_t) :
            abs(lin.dqf_dvf) + abs(lin.dqf_dvt) + abs(lin.dqf_dtheta_f) + abs(lin.dqf_dtheta_t)
    end
end

function _atlambda_l1(ws, f::AbstractString)::Float64
    total = 0.0
    for (cls, st) in pairs(ws.vars)
        lam = st.lambda
        if cls in _PM1_CLASSES
            total += 2.0 * sum(abs, lam)
        elseif cls === :p_flow || cls === :q_flow
            which = cls === :p_flow ? :p : :q
            for (k, pair) in enumerate(ws.interzonal_line_pairs)
                row = 0.0
                for (zp, le) in ((pair.zone_i, pair.line_local_i),
                                 (pair.zone_j, pair.line_local_j))
                    row += _flow_row_l1(ws, ws.zones[zp], le, f, which)
                end
                total += abs(lam[k]) * row
            end
        elseif isempty(lam)
            # empty dummy-class slot (e.g. tie_line): zero contribution
            continue
        else
            return NaN   # populated dummy/unknown class: not covered
        end
    end
    return total
end

# ---- algorithm builders ------------------------------------------------
function build_gold_algorithm(algo, f, s, c, blk)
    if algo == "admm"
        return Lirion.admm(; coupling = Symbol(c), sectioning = Symbol(s),
            max_iter = Int(blk["max_iter"]), rho = Float64(blk["rho"]),
            eps_abs = 1e-4, eps_rel = 1e-4, adaptive_rho = Bool(blk["adaptive"]))
    else
        rho = Float64(blk["rho"]); mu = Float64(blk["mu"])
        kw = (; coupling = Symbol(c), sectioning = Symbol(s),
            max_iter = Int(blk["max_iter"]), rho = rho, mu = mu,
            rho_growth = Float64(blk["rho_growth"]), mu_growth = Float64(blk["mu_growth"]),
            adaptive = Bool(blk["adaptive"]), proximal_mode = Symbol(blk["mode"]),
            rho_upper = max(1e10, rho), mu_upper = max(1e10, mu),
            eps_abs = 1e-4, eps_rel = 1e-4,
            alpha3 = Float64(get(blk, "alpha3", 1.0)),
            globalization = Symbol(get(blk, "globalization", "armijo_l2")),
            qp_optimizer = aladin_coord_for(f))
        # per-cell fixed λ̄ (socwrnl only; 0.0 = auto)
        return f == "socwrnl" ?
            Lirion.aladin(; kw..., alg3_lam_bar = Float64(get(blk, "lam_bar", 0.0))) :
            Lirion.aladin(; kw...)
    end
end

# ---- history extraction ----------------------------------------------
_histcol(h, name) = hasproperty(h, name) ? getproperty(h, name) : Float64[]

function history_rows(hist, is_aladin)
    n = length(hist.r_primal)
    rows = NamedTuple[]
    rd = _histcol(hist, :r_dual); ed = _histcol(hist, :eps_dual)
    rs = _histcol(hist, :r_step_primal)
    # ALADIN: lambda_step = ‖Δλ‖, r_dual = μ·lambda_step
    ls = _histcol(hist, :lambda_step)
    es = _histcol(hist, :eps_step_primal)
    # per-iteration instrumentation series (hist.meta)
    meta_d = hasproperty(hist, :meta) ? hist.meta : Dict{Symbol, Any}()
    sd_v = get(meta_d, :sd, Float64[])::Vector{Float64}
    pd_v = get(meta_d, :psd_delta, Float64[])::Vector{Float64}
    for i in 1:n
        push!(rows, (
            iter = i,
            r_primal = hist.r_primal[i],
            eps_primal = hist.eps_primal[i],
            r_dual = i <= length(rd) ? rd[i] : NaN,
            eps_dual = i <= length(ed) ? ed[i] : NaN,
            r_step_primal = (is_aladin && i <= length(rs)) ? rs[i] : NaN,
            eps_step_primal = (is_aladin && i <= length(es)) ? es[i] : NaN,
            lambda_step = (is_aladin && i <= length(ls)) ? ls[i] : NaN,
            sd = (is_aladin && i <= length(sd_v)) ? sd_v[i] : NaN,
            psd_delta = (is_aladin && i <= length(pd_v)) ? pd_v[i] : NaN,
        ))
    end
    return rows
end

const HIST_COLUMNS = [:iter, :r_primal, :eps_primal, :r_dual, :eps_dual,
                      :r_step_primal, :eps_step_primal, :lambda_step, :sd, :psd_delta]

# ---- one gold run, low level (keeps ws for history + duals) -----------
function run_cell(case, k, f, s, c, algo, blk)
    path = partitioned_case_path(case, "std", k)
    alg = build_gold_algorithm(algo, f, s, c, blk)
    is_aladin = algo != "admm"

    # ALADIN toggles come from the parameter block; all restored in the finally
    save_reg       = Lirion._QP_HESS_REG[]
    save_noprox_ac = Lirion._AC_QP_NO_PROX[]
    save_noprox_sn = Lirion._SOCWRNL_QP_NO_PROX[]
    save_modgrad   = Lirion._MODIFIED_GRADIENT[]
    save_bdetect   = Lirion._BOUND_DETECT_MODE[]

    # lifted flow: on only for blocks with lifted_flow=true, otherwise forced unset
    want_lifted  = get(blk, "lifted_flow", false) == true
    had_lifted   = haskey(ENV, "LIRION_SOCWRNL_LIFTED_FLOW")
    save_lifted  = get(ENV, "LIRION_SOCWRNL_LIFTED_FLOW", "")
    # lifted flow is socwrnl-ADMM only
    want_lifted && (algo != "admm" || f != "socwrnl") && error(
        "repro guard: lifted_flow=true on a $(f)/$(algo) cell ($(case) k$(k) $(s)/$(c)); " *
        "LIRION_SOCWRNL_LIFTED_FLOW is SOCWRNL-ADMM only.")

    # coordinator-QP proximal reg is part of the parameters (ALADIN only)
    Lirion._QP_HESS_REG[] = is_aladin ? Float64(get(blk, "qp_reg", 0.0)) : 0.0

    if is_aladin
        if f == "socwrnl"
            # canonical ⇒ no_prox/modgrad true + primal_dual; prox_in_B ⇒ false/:primal
            canon = get(blk, "prox_mode", "prox_in_B") == "canonical"
            Lirion._SOCWRNL_QP_NO_PROX[] = canon
            Lirion._MODIFIED_GRADIENT[]  = canon
            Lirion._BOUND_DETECT_MODE[]  =
                Symbol(get(blk, "bound_detect", canon ? "primal_dual" : "primal"))
        elseif f == "ac"
            # paper cells: true/true/primal_dual
            Lirion._AC_QP_NO_PROX[]     = Bool(get(blk, "no_prox", true))
            Lirion._MODIFIED_GRADIENT[] = Bool(get(blk, "modgrad", true))
            Lirion._BOUND_DETECT_MODE[] = Symbol(get(blk, "bound_detect", "primal_dual"))
        elseif f == "dc"
            # inert on DC; the bound detector is what matters
            Lirion._AC_QP_NO_PROX[]     = Bool(get(blk, "no_prox", false))
            Lirion._MODIFIED_GRADIENT[] = Bool(get(blk, "modgrad", false))
            Lirion._BOUND_DETECT_MODE[] = Symbol(get(blk, "bound_detect", "primal_dual"))
        end
        # other formulations: module defaults
    end

    # scope the lifted-flow env to this cell (before the manifest)
    if want_lifted
        ENV["LIRION_SOCWRNL_LIFTED_FLOW"] = "1"
    else
        delete!(ENV, "LIRION_SOCWRNL_LIFTED_FLOW")
    end

    # manifest captured after the toggles are set, before the reset
    manifest = runtime_manifest(f, algo; alg = alg)

    status = ""; obj = NaN; iters = 0; rt = NaN; fail = ""
    lam1 = NaN; kap1 = NaN; mn = 0; sd = NaN; atl = NaN; npx = 0
    hrows = NamedTuple[]
    # what was actually applied (ALADIN only)
    runm = Dict{String, Any}()
    t_wall = time()
    try
        inputs = Lirion._build_inputs(path, alg; model = model_ctor(f)(),
            optimizer = official_dist_optimizer(f, algo), outputs = false)
        t0 = time_ns()
        parts, ws, iter_done = redirect_stdout(devnull) do
            with_logger(NullLogger()) do
                Lirion.solve_distributed(inputs)
            end
        end
        r = redirect_stdout(devnull) do
            with_logger(NullLogger()) do
                Lirion.build_lirion_result_distributed(
                    parts, ws, inputs.model_type, iter_done, t0, inputs.system)
            end
        end
        rt = time() - t_wall
        obj = Float64(r.objective); status = string(r.status)
        iters = iter_done === nothing ? 0 : iter_done
        hrows = history_rows(ws.hist, algo != "admm")
        if algo != "admm"
            try
                lam1, kap1, m, n = dual_norms_and_counts(ws)
                mn = m + n   # m constraints + n variables
                sd = max(S_MAX, (lam1 + kap1) / mn) / S_MAX
            catch e
                fail = "s_d extraction failed: " * failure_reason_of(e)
            end
            atl = try
                _atlambda_l1(ws, f)
            catch
                NaN
            end
            # n_prox into the meta
            npx = try
                Lirion._structural_nprox(ws)
            catch
                0
            end
            # inputs.model_type is the resolved ALADINConfig
            runm = try
                Lirion.aladin_run_manifest(ws, inputs.model_type)
            catch e
                Dict{String, Any}("error" => failure_reason_of(e))
            end
        end
    catch e
        rt = time() - t_wall
        status = "error"; fail = failure_reason_of(e)
    finally
        # restore EVERY togglable Ref and the lifted-flow env
        Lirion._QP_HESS_REG[]        = 0.0
        Lirion._AC_QP_NO_PROX[]      = save_noprox_ac
        Lirion._SOCWRNL_QP_NO_PROX[] = save_noprox_sn
        Lirion._MODIFIED_GRADIENT[]  = save_modgrad
        Lirion._BOUND_DETECT_MODE[]  = save_bdetect
        if had_lifted
            ENV["LIRION_SOCWRNL_LIFTED_FLOW"] = save_lifted
        else
            delete!(ENV, "LIRION_SOCWRNL_LIFTED_FLOW")
        end
    end
    return (; status, obj, iters, rt, fail, lam1, kap1, mn, sd, atl, npx, hrows,
              manifest, runm)
end

# ---- runtime manifest -------------------------------------------------
# every toggle that can alter a solve, dumped into the meta TOML
function _attrs_of(opt)
    opt isa JuMP.MOI.OptimizerWithAttributes || return ""
    parts = String[]
    for p_ in opt.params
        nm = p_.first isa JuMP.MOI.RawOptimizerAttribute ? p_.first.name : string(p_.first)
        nm in ("verbose", "print_level", "sb", "print_timing_statistics",
               "print_user_options", "file_print_level", "OutputFlag") && continue
        push!(parts, "$(nm)=$(p_.second)")
    end
    return join(parts, ";")
end

function _resolved_solver_name(opt, f::AbstractString = "")
    # record RESOLVED attributes (tolerances etc.), never a bare label
    if opt === nothing
        # nothing = Lirion default: socwrnl -> Ipopt+MA57, socwr -> Clarabel
        if f == "socwrnl"
            resolved = try
                Lirion.resolve_optimizer(Val(:socwrnl), nothing)
            catch e
                return "UNRESOLVED socwrnl default: " *
                       first(split(sprint(showerror, e), "\n"))
            end
            at = _attrs_of(resolved)
            nm = try
                m = JuMP.Model(resolved); JuMP.set_silent(m); JuMP.solver_name(m)
            catch
                "Ipopt"
            end
            return nm * "[" * at * "] (Lirion _socwrnl_tight_factory)"
        end
        base = Lirion._SILENT_FACTORY_SOCWR
        return "Clarabel[" * _attrs_of(base) * "] (Lirion _SILENT_FACTORY_SOCWR)"
    end
    name = try
        m = JuMP.Model(opt)
        JuMP.set_silent(m)
        JuMP.solver_name(m)
    catch e
        "UNRESOLVED: " * first(split(sprint(showerror, e), "\n"))
    end
    at = _attrs_of(opt)
    return isempty(at) ? name * "[defaults]" : name * "[" * at * "]"
end

function runtime_manifest(f::AbstractString, algo::AbstractString; alg = nothing)
    zones = _resolved_solver_name(official_dist_optimizer(f, algo), f)
    coord = algo == "admm" ? "n/a" : _resolved_solver_name(aladin_coord_for(f))
    m = Dict{String, Any}(
        "SOLVER_ZONES_RESOLVED" => zones,
        "SOLVER_COORD_QP_RESOLVED" => coord,
        "STRICT_SOLVERS" => get(ENV, "LIRION_STRICT_SOLVERS", "0"),
        # opt-in socwrnl ADMM lifted flow; "unset" on the default path
        "LIRION_SOCWRNL_LIFTED_FLOW" => get(ENV, "LIRION_SOCWRNL_LIFTED_FLOW", "unset"),
    )
    # module-level ALADIN toggles
    merge!(m, Lirion.aladin_runtime_manifest())
    # resolve the config the way _build_inputs does
    cfg = alg isa Lirion.ALADINAlgorithm ? Lirion._make_config(alg, model_ctor(f)()) : alg
    if cfg isa Lirion.ALADINConfig
        m["GLOBALIZATION"] = String(cfg.globalization)
        m["ALG3_LAM_BAR_CFG"] = Float64(cfg.alg3_lam_bar)
    end
    return m
end

# ---- post-step: crossing table + comparison with the expected files ------
# separate julia processes; skipped with LIRION_REPRO_NO_TABLE=1
function post_step(gold_paths)
    jl = Base.julia_cmd()
    proj = Base.active_project()
    projflag = proj === nothing ? "--project" : "--project=$(dirname(proj))"
    build = joinpath(@__DIR__, "scripts", "build_crossing_table.jl")
    verify = joinpath(@__DIR__, "scripts", "verify_crossing.jl")
    @info "building the crossing table" hist = OUT_DIR out = TABLE_DIR
    run(`$jl $projflag $build $OUT_DIR $CENTRAL_CSV $TABLE_DIR`)
    expected = [joinpath(dirname(g), "expected_crossing.csv") for g in gold_paths]
    expected = filter(isfile, expected)
    isempty(expected) && (@warn "no expected_crossing.csv next to the parameter files; comparison skipped"; return)
    ok = success(`$jl $projflag $verify $(joinpath(TABLE_DIR, "crossing_table.csv")) $expected`)
    ok || @error "the regenerated crossing table differs from the expected one — see the lines above"
end

# ---- main -------------------------------------------------------------
function main()
    gold_paths = gold_files()
    cells = Any[]
    for g in gold_paths
        isfile(g) || error("missing parameter file $g")
        append!(cells, get(TOML.parsefile(g), "cell", Any[]))
    end
    mkpath(OUT_DIR)
    n_run = 0; n_err = 0; n_ph = 0
    for cell in cells
        case = cell["case"]; k = Int(cell["k"]); f = cell["formulation"]
        s = cell["sectioning"]; c = cell["coupling"]
        (isempty(F_CASE) || occursin(F_CASE, case)) || continue
        (isempty(F_K) || string(k) == F_K) || continue
        (isempty(F_FORM) || f == F_FORM) || continue
        (isempty(F_SECT) || s == F_SECT) || continue
        (isempty(F_COUP) || c == F_COUP) || continue
        # "admm" + every ALADIN block ("aladin", and future "aladin_pflow" etc.)
        for algo_key in sort(collect(keys(cell)))
            (algo_key == "admm" || startswith(algo_key, "aladin")) || continue
            blk = cell[algo_key]
            blk isa AbstractDict || continue
            algo = algo_key == "admm" ? "admm" : "aladin"
            mode = algo == "admm" ? "" : get(blk, "mode", "primal")
            stem = join(filter(!isempty,
                [case, "k$k", f, s, c, algo_key, mode == "primal" ? "" : mode]), "_")

            # ---- non-run statuses: placeholder meta (no solve) ----
            cell_status = get(blk, "status", "pass")
            if lowercase(string(cell_status)) in ("non_convergent", "pending")
                write_toml(joinpath(OUT_DIR, stem * ".toml"), Dict(
                    "case" => case, "k" => k, "formulation" => f,
                    "sectioning" => s, "coupling" => c, "algorithm" => algo,
                    "algo_key" => algo_key, "mode" => mode,
                    "params" => Dict(String(kk) => vv for (kk, vv) in blk),
                    "status" => string(cell_status), "objective" => NaN,
                    "iterations" => 0, "runtime_s" => 0.0,
                    "failure" => get(blk, "note",
                        "declared $(cell_status) in gold_params.toml"),
                    "git_commit" => GIT_COMMIT,
                ))
                n_ph += 1
                @info @sprintf("%-16s k%-2d %-7s %-10s %-7s %-13s %s (placeholder)",
                    case, k, f, s, c, algo_key, string(cell_status))
                continue
            end

            res = run_cell(case, k, f, s, c, algo, blk)
            n_run += 1; res.status == "error" && (n_err += 1)
            write_csv(joinpath(OUT_DIR, stem * ".csv"), res.hrows; columns = HIST_COLUMNS)
            write_toml(joinpath(OUT_DIR, stem * ".toml"), Dict(
                "case" => case, "k" => k, "formulation" => f,
                "sectioning" => s, "coupling" => c, "algorithm" => algo,
                "algo_key" => algo_key, "mode" => mode,
                "params" => Dict(String(kk) => vv for (kk, vv) in blk),
                "status" => res.status, "objective" => res.obj,
                "iterations" => res.iters, "runtime_s" => res.rt,
                "failure" => res.fail,
                "lambda_l1" => res.lam1, "kappa_l1" => res.kap1,
                "m_plus_n" => res.mn, "s_d" => res.sd, "s_max" => S_MAX,
                "atlambda_l1" => res.atl,
                "n_prox" => res.npx,
                "manifest" => res.manifest,
                "run_manifest" => res.runm,
                "git_commit" => GIT_COMMIT,
            ))
            @info @sprintf("%-16s k%-2d %-7s %-10s %-7s %-13s %-14s it=%-6d obj=%.6g s_d=%s %s",
                case, k, f, s, c, algo_key, res.status, res.iters, res.obj,
                isnan(res.sd) ? "-" : @sprintf("%.3f", res.sd), res.fail)
        end
    end
    (n_run == 0 && n_ph == 0) && error("no cell matched the filters")
    @info "run.jl done" runs = n_run errors = n_err placeholders = n_ph out = OUT_DIR
    get(ENV, "LIRION_REPRO_NO_TABLE", "0") == "1" || post_step(gold_paths)
end
main()
