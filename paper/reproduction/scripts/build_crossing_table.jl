# Post-hoc crossing table: reads the histories from run.jl and the oracle
# objectives, finds the first Tier-1/Tier-2 iteration per (cell, algorithm).
# julia --project=<repo> scripts/build_crossing_table.jl [hist_dir] [centralized_csv] [out_dir]

include(joinpath(@__DIR__, "lib", "experiment_common.jl"))
using TOML

const HIST_DIR = length(ARGS) >= 1 ? ARGS[1] : normpath(joinpath(@__DIR__, "..", "results", "runs"))
const CENTRAL_CSV = length(ARGS) >= 2 ? ARGS[2] :
    joinpath(@__DIR__, "repro_centralized_baseline.csv")
const OUT_DIR = length(ARGS) >= 3 ? ARGS[3] : normpath(joinpath(@__DIR__, "..", "results", "table"))

const TIER1_MISMATCH = 0.01   # PMADA: ℓ2-norm of the mismatch < 0.01 (rad, p.u.)
const GAP_TOL_PCT = 1.0       # |gap| < 1% for either tier's PASS

# ---- oracle objectives ------------------------------------------------
# only verified optima count as gap references; anything else -> NaN gap
const ORACLE_OK_STATUS = Set(["optimal", "local_optimal", "locally_solved"])

function load_oracle(path)
    d = Dict{Tuple{String, String}, Float64}()
    isfile(path) || error("missing centralized baseline: $path")
    lines = readlines(path)
    hdr = split(lines[1], ",")
    ic = findfirst(==("case"), hdr); iff = findfirst(==("formulation"), hdr)
    io = findfirst(==("objective"), hdr); is = findfirst(==("status"), hdr)
    (ic === nothing || iff === nothing || io === nothing || is === nothing) &&
        error("unexpected header in $path")
    for ln in lines[2:end]
        fl = split(ln, ",")
        length(fl) >= max(io, is) || continue
        obj = tryparse(Float64, fl[io])
        obj === nothing && continue
        st = lowercase(String(fl[is]))
        if !(st in ORACLE_OK_STATUS)
            @warn "oracle row excluded (unverified status)" case = fl[ic] formulation = fl[iff] status = fl[is]
            continue
        end
        d[(String(fl[ic]), String(fl[iff]))] = obj
    end
    return d
end

# ---- power (relaxation) oracle ---------------------------------------
# optional power_oracle_baseline.csv next to the centralized baseline
function load_power_oracle(central_csv)
    path = joinpath(dirname(central_csv), "power_oracle_baseline.csv")
    d = Dict{Tuple{String, String, String}, Float64}()
    isfile(path) || return d
    lines = readlines(path)
    hdr = split(lines[1], ",")
    col = Dict(h => i for (i, h) in enumerate(hdr))
    for ln in lines[2:end]
        fl = split(ln, ",")
        length(fl) >= 4 || continue
        obj = tryparse(Float64, fl[col["objective"]])
        obj === nothing && continue
        st = lowercase(String(fl[col["status"]]))
        st in ORACLE_OK_STATUS || continue
        d[(String(fl[col["case"]]), String(fl[col["formulation"]]),
           String(fl[col["sectioning"]]))] = obj
    end
    return d
end

# ---- history loading --------------------------------------------------
function load_history(path)
    lines = readlines(path)
    hdr = Symbol.(split(lines[1], ","))
    cols = Dict(h => Float64[] for h in hdr)
    for ln in lines[2:end]
        fl = split(ln, ",")
        for (j, h) in enumerate(hdr)
            v = j <= length(fl) ? something(tryparse(Float64, fl[j]), NaN) : NaN
            push!(cols[h], v)
        end
    end
    return cols
end

_first_iter(pred, n) = something(findfirst(pred, 1:n), 0)

# Start of the final contiguous block over which `pred` holds through the last
# iteration; 0 if the last iteration fails (no sustained crossing).
function _sustained_iter(pred, n)
    n == 0 && return 0
    pred(n) || return 0
    s = n
    while s > 1 && pred(s - 1)
        s -= 1
    end
    return s
end

# ---- main -------------------------------------------------------------
function main()
    isdir(HIST_DIR) || error("missing history dir: $HIST_DIR — run run.jl first")
    oracle = load_oracle(CENTRAL_CSV)
    poracle = load_power_oracle(CENTRAL_CSV)
    rows = NamedTuple[]
    metas = sort(filter(f -> endswith(f, ".toml"), readdir(HIST_DIR)))
    for mf in metas
        meta = TOML.parsefile(joinpath(HIST_DIR, mf))
        case = meta["case"]; k = Int(meta["k"]); f = meta["formulation"]
        s = meta["sectioning"]; c = meta["coupling"]; algo = meta["algorithm"]
        cen = get(oracle, (case, f), NaN)
        isnan(cen) && @warn "no verified oracle for (case, formulation) — gap=NaN, tiers fail-safe" case f meta = mf
        obj = Float64(get(meta, "objective", NaN))
        gap = (isnan(obj) || isnan(cen) || cen == 0) ? NaN : 100 * (obj - cen) / abs(cen)
        gap_ok = !isnan(gap) && abs(gap) < GAP_TOL_PCT

        hist_path = joinpath(HIST_DIR, replace(mf, ".toml" => ".csv"))
        h = isfile(hist_path) ? load_history(hist_path) : nothing
        n = h === nothing ? 0 : length(h[:iter])

        t1 = 0; t2 = 0; t1_touch = 0; t2_touch = 0
        t1_ok = false; t2_ok = false
        # t2_rel NaN = no Aᵀλ in the meta -> blank column
        t2_raw = 0; t2_raw_ok = false; t2_rel = NaN
        t2_raw_touch = 0; sust_frac = NaN
        estac_final = NaN; estac_min = NaN
        sd = Float64(get(meta, "s_d", NaN))
        if n > 0
            converged = meta["status"] == "converged"
            # iteration_limit runs must hold over the last min(50, 1% of n) iterations
            min_block = converged ? 1 : min(50, max(1, round(Int, 0.01 * n)))
            rp = h[:r_primal]; ep = h[:eps_primal]
            p1 = i -> rp[i] <= TIER1_MISMATCH
            t1_touch = _first_iter(p1, n)
            t1 = _sustained_iter(p1, n)
            t1_ok = t1 > 0 && (n - t1 + 1) >= min_block
            p2 = nothing
            # stationarity bar: eps_step_primal history, else n_prox from the meta
            es_hist = haskey(h, :eps_step_primal) ? h[:eps_step_primal] : Float64[]
            np_meta = Float64(get(meta, "n_prox", NaN))
            # no silent fallback: a run without a bar is refused
            bar_at = i -> begin
                if i <= length(es_hist) && !isnan(es_hist[i]) && es_hist[i] > 0
                    es_hist[i]
                elseif !isnan(np_meta)
                    np_meta * 1e-4
                else
                    error("derived stationarity bar unavailable for $mf " *
                          "(no eps_step_primal in history, no n_prox in meta) — " *
                          "refusing silent legacy 1e-2; re-run or annotate the meta")
                end
            end
            # rel-term column: bar + ε_rel·‖Aᵀλ‖₁ (blank without atlambda_l1)
            atl_meta = Float64(get(meta, "atlambda_l1", NaN))
            relterm = isnan(atl_meta) ? NaN : 1e-4 * atl_meta
            bar_rel_at = i -> bar_at(i) + relterm
            if algo == "admm"
                rd = h[:r_dual]; ed = h[:eps_dual]
                p2 = i -> rp[i] <= ep[i] && rd[i] <= ed[i]
            else
                rs = h[:r_step_primal]
                estac_final = rs[n]
                # init=Inf (NaN would poison the min); iteration 1 records estac = 0
                estac_min = minimum(x for x in rs if !isnan(x) && x > 0; init = Inf)
                isfinite(estac_min) || (estac_min = NaN)
                if !isnan(sd)
                    # skip iteration 1 (estac = 0 artifact)
                    p2 = i -> i > 1 && rp[i] <= ep[i] && rs[i] / sd <= bar_at(i)
                end
            end
            if p2 !== nothing
                t2_touch = _first_iter(p2, n)
                t2 = _sustained_iter(p2, n)
                t2_ok = t2 > 0 && (n - t2 + 1) >= min_block
            end
            # raw column: ADMM raw ≡ scaled (no s_d); ALADIN raw uses s_d=1
            p2r = nothing
            if algo == "admm"
                t2_raw = t2; t2_raw_ok = t2_ok
                t2_rel = t2
                p2r = p2
            else
                rs = h[:r_step_primal]
                p2r = i -> i > 1 && rp[i] <= ep[i] && rs[i] <= bar_at(i)
                t2_raw = _sustained_iter(p2r, n)
                t2_raw_ok = t2_raw > 0 && (n - t2_raw + 1) >= min_block
                # symmetric rel-term column (post-hoc only; blank without Aᵀλ)
                if !isnan(relterm)
                    p2rel = i -> i > 1 && rp[i] <= ep[i] && rs[i] <= bar_rel_at(i)
                    t2_rel = _sustained_iter(p2rel, n)
                end
            end
            # first touch on the raw column = the published Tier-2 iteration
            if p2r !== nothing
                t2_raw_touch = _first_iter(p2r, n)
                if t2_raw_touch > 0
                    sust_frac = count(p2r, t2_raw_touch:n) / (n - t2_raw_touch + 1)
                end
            end
        end

        rp_final = n > 0 ? h[:r_primal][n] : NaN
        ep_final = n > 0 ? h[:eps_primal][n] : NaN
        # power coupling = relaxation: its gap measures slack, not quality
        relaxation = (c == "power")
        # gap vs the relaxed oracle when available
        pcen = get(poracle, (case, f, s), NaN)
        gap_relax = (relaxation && !isnan(obj) && !isnan(pcen) && pcen != 0) ?
            100 * (obj - pcen) / abs(pcen) : NaN
        optimality_gated = !relaxation || !isnan(gap_relax)
        # convex formulation beating its oracle by >0.1%: flag it
        red_flag = meta["status"] == "converged" && f in ("dc", "socwr", "socwrnl") &&
                   ((!relaxation && !isnan(gap) && gap < -0.1) ||
                    (relaxation && !isnan(gap_relax) && gap_relax < -0.1))

        # power rows gate vs the relaxed oracle only; non-power vs the full one
        opt_ok = relaxation ?
            (isnan(gap_relax) ? true : abs(gap_relax) < GAP_TOL_PCT) : gap_ok

        # zone-solver KKT budget per solver vs the bar; ADMM is not gated
        npx_meta = Float64(get(meta, "n_prox", NaN))
        barra_meta = npx_meta * 1e-4
        kap1 = Float64(get(meta, "kappa_l1", NaN))
        is_aladin = startswith(meta["algo_key"], "aladin")
        # socwrnl zones run Ipopt at dual_inf_tol=1e-9
        kkt_budget = !is_aladin ? NaN :
            f == "socwr" ? 1e-4 * kap1 :
            f == "socwrnl" ? 1e-9 * kap1 :
            f == "dc" ? 1e-6 * npx_meta :
            f == "ac" ? 1e-7 * sd : NaN
        # fail closed: an unmeasurable budget never reaches "pass"
        is_conic_aladin = is_aladin && (f == "socwr" || f == "socwrnl")
        budget_unmeasurable = is_conic_aladin && (isnan(kkt_budget) || isnan(barra_meta))
        gated_nc = budget_unmeasurable ||
                   (is_aladin && !isnan(kkt_budget) && !isnan(barra_meta) &&
                    kkt_budget > barra_meta)
        verdict, verdict_reason = if meta["status"] == "non_convergent"
            ("non_convergent", get(meta, "failure", ""))
        elseif gated_nc
            ("nao_certificavel",
             budget_unmeasurable ?
             @sprintf("B4: kkt_budget nao-mensuravel (‖κ‖₁ ou barra = NaN) para ALADIN %s — certificado de estacionariedade indisponivel (fail-closed)",
                      f) :
             @sprintf("B4: kkt_budget medido %.3g > barra %.3g — certificado de estacionariedade invalido (convergencia frouxa = artefato)",
                      kkt_budget, barra_meta))
        else
            # first-touch verdict; the sustained columns are diagnostics only
            pass_full = t1_touch > 0 && t2_raw_touch > 0 && opt_ok
            (pass_full ? "pass" : "fail", "")
        end
        push!(rows, (
            case = case, k = k, formulation = f, sectioning = s, coupling = c,
            algorithm = meta["algo_key"], mode = get(meta, "mode", ""),
            verdict = verdict, verdict_reason = verdict_reason,
            kkt_budget = kkt_budget,
            status = meta["status"], objective = obj, centralized = cen,
            gap_pct = gap,
            r_primal_final = rp_final, eps_primal_final = ep_final,
            # *_pass are first-touch; tier*_iter are the sustained diagnostics
            tier1_iter = t1, tier1_pass = (t1_touch > 0 && opt_ok),
            tier2_iter = t2,
            tier2_pass = (t2_raw_touch > 0 && opt_ok),
            tier2_iter_raw = t2_raw,
            tier2_pass_sustained = (t1_ok && t2_ok && t2_raw_ok && opt_ok),
            tier2_iter_reltol = t2_rel,
            tier1_first_touch = t1_touch, tier2_first_touch = t2_touch,
            tier2_first_touch_raw = t2_raw_touch,
            sustained_frac = sust_frac,
            iters_total = Int(get(meta, "iterations", 0)),
            estac_final = estac_final, estac_min = estac_min, s_d = sd,
            relaxation = relaxation, gap_relax_pct = gap_relax,
            optimality_gated = optimality_gated, red_flag = red_flag,
            runtime_s = Float64(get(meta, "runtime_s", NaN)),
            failure = get(meta, "failure", ""),
        ))
    end
    isempty(rows) && error("no run metas found under $HIST_DIR")
    sort!(rows, by = r -> (r.case, r.k, r.formulation, r.sectioning, r.coupling, r.algorithm))
    cols = [:case, :k, :formulation, :sectioning, :coupling, :algorithm, :mode,
            :verdict, :verdict_reason, :kkt_budget,
            :status, :objective, :centralized, :gap_pct,
            :r_primal_final, :eps_primal_final,
            :tier2_first_touch_raw, :tier2_pass, :tier1_first_touch, :tier1_pass,
            :tier2_first_touch, :sustained_frac,
            :tier1_iter, :tier2_iter, :tier2_iter_raw, :tier2_pass_sustained,
            :tier2_iter_reltol, :iters_total,
            :estac_final, :estac_min, :s_d, :relaxation, :gap_relax_pct,
            :optimality_gated, :red_flag,
            :runtime_s, :failure]
    mkpath(OUT_DIR)
    write_csv(joinpath(OUT_DIR, "crossing_table.csv"), rows; columns = cols)
    write_markdown(joinpath(OUT_DIR, "crossing_table.md"), rows; columns = cols,
        title = "Tier-1 / Tier-2 crossing table (post-hoc, frozen criterion) — " *
                "power-coupling gaps measure the interzonal-KVL relaxation slack, " *
                "not algorithm quality (relaxation = true rows)")
    n1 = count(r -> r.tier1_pass, rows); n2 = count(r -> r.tier2_pass, rows)
    nflag = count(r -> r.red_flag, rows)
    npass = count(r -> r.verdict == "pass", rows)
    nnc = count(r -> r.verdict == "nao_certificavel", rows)
    nnconv = count(r -> r.verdict == "non_convergent", rows)
    nfail = count(r -> r.verdict == "fail", rows)
    @info @sprintf("crossing table: %d runs | verdict: pass %d, nao_certificavel %d, non_convergent %d, fail %d | Tier-1 %d | Tier-2 %d | red flags %d | out=%s",
        length(rows), npass, nnc, nnconv, nfail, n1, n2, nflag, OUT_DIR)
end
main()
