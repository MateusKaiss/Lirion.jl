# ============================================================
# paper/scripts/experiments/run_case14_gate.jl
# ============================================================
# Stage D1a — package model review + case14 gate.
# Runs ONLY case14_ieee (std). Computes the case14 centralized oracle
# (official solvers), then tests every supported distributed cell
# (formulation × sectioning × coupling) with BOTH ADMM and ALADIN, with
# robust failure capture, and classifies each cell (both_converged / ...)
# plus a paper_status. Also tests rho=0 as a diagnostic and pflow.
#
# Official runtime: JULIA_NUM_THREADS=8 (env) + BLAS=1.
# Official distributed solvers: DC=OSQP, SOCWR=Clarabel(hardened),
# AC=Ipopt+MA57. Centralized oracle: DC=Gurobi, SOCWR=Clarabel(hardened),
# AC=Ipopt+MA57. tol = 1e-4.
#
# Outputs:
#   results/tuning/case14_gate_raw.csv            (scratch)
#   artifacts/tables/D1a_case14_gate.csv / .md    (curated)
#   artifacts/model_support_matrix.toml           (curated)
# ============================================================

include(joinpath(@__DIR__, "lib", "tuning_common.jl"))

using LinearAlgebra
using Ipopt, Clarabel
using JuMP: optimizer_with_attributes
BLAS.set_num_threads(1)

const RES_TUNING = joinpath(RESULTS_DIR, "tuning")
const CASE = "case14_ieee"
const KS = [2, 3, 4]
const FORMS = ["dc", "socwr", "ac"]
const SECTS = ["tie_line", "dummy_bus", "split_bus"]
const COUPS = ["state", "hybrid", "power"]
const GAP_TARGET = 1.0
const ADMM_MAXIT = 500
const ALADIN_MAXIT = 200

# ---- Centralized oracle for case14 (official solvers) ----------------
const _GRB = try; @eval import Gurobi; Gurobi; catch; nothing; end
_grbenv = nothing
function _gurobi_opt()
    _GRB === nothing && return optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
    try
        global _grbenv; _grbenv === nothing && (_grbenv = _GRB.Env())
        return optimizer_with_attributes(() -> _GRB.Optimizer(_grbenv), "OutputFlag"=>0)
    catch; return optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0); end
end
_clarabel_hard() = optimizer_with_attributes(Clarabel.Optimizer, "verbose"=>false,
    "tol_feas"=>1e-4, "tol_gap_abs"=>1e-4, "tol_gap_rel"=>1e-4, "max_iter"=>1200,
    "direct_solve_method"=>:qdldl, "equilibrate_enable"=>true)
_ipopt_ac() = have_hsl() ?
    optimizer_with_attributes(Ipopt.Optimizer,"print_level"=>0,"linear_solver"=>"ma57","hsllib"=>hsllib_path()) :
    optimizer_with_attributes(Ipopt.Optimizer,"print_level"=>0)
_central_opt(f) = f=="dc" ? _gurobi_opt() : (f=="socwr" ? _clarabel_hard() : _ipopt_ac())

function compute_oracle()
    oracle = Dict{String,Float64}()
    for f in FORMS
        try
            r = Lirion.solve(raw_case_path(CASE,"std"); model=model_ctor(f)(),
                algorithm=Lirion.centralized(), optimizer=_central_opt(f), outputs=false, plots=false)
            oracle[f] = Float64(r.objective)
            @info @sprintf("  central %-6s -> %-14s obj=%.4f", f, string(r.status), r.objective)
        catch e
            @warn "central failed" f exc=first(split(sprint(showerror,e),"\n"))
        end
    end
    return oracle
end

# ---- Algorithm build (supports proximal_mode for pflow) --------------
function build_alg(algo, sect, coup, p; proximal_mode=nothing)
    if algo == "admm"
        return Lirion.admm(; coupling=Symbol(coup), sectioning=Symbol(sect),
            max_iter=ADMM_MAXIT, rho=p.rho, eps_abs=1e-4, eps_rel=1e-4, adaptive_rho=p.adapt)
    else
        kw = (; coupling=Symbol(coup), sectioning=Symbol(sect), max_iter=ALADIN_MAXIT,
              rho=p.rho, mu=p.mu, rho_growth=p.gr)
        proximal_mode !== nothing && (kw = merge(kw, (; proximal_mode=Symbol(proximal_mode))))
        return Lirion.aladin(; kw...)
    end
end

# ---- One robust attempt ----------------------------------------------
function attempt(k, f, algo, sect, coup, p, aid, oracle; proximal_mode=nothing)
    path = partitioned_case_path(CASE, "std", k); cen = get(oracle, f, NaN)
    status=""; obj=NaN; iters=0; pr=NaN; dr=NaN; rt=NaN; conv=false; support=""; fail=""
    coup_lbl = proximal_mode === nothing ? coup : "$(coup)+$(proximal_mode)"
    plabel = algo=="admm" ? @sprintf("rho=%.0e,adapt=%s",p.rho,p.adapt) :
             @sprintf("rho=%.0e,mu=%.0e,gr=%.2f",p.rho,p.mu,p.gr)
    t0=time()
    try
        alg = build_alg(algo, sect, coup, p; proximal_mode=proximal_mode)
        r = Lirion.solve(path; model=model_ctor(f)(), algorithm=alg,
            optimizer=official_dist_optimizer(f), outputs=false, plots=false)
        rt=time()-t0; obj=Float64(r.objective); status=string(r.status)
        iters = r.iter_done===nothing ? 0 : r.iter_done
        pr,dr = extract_residuals(r); conv=(r.status==:converged)
        support = isfinite(obj) ? "supported" : "unstable_or_fails"
    catch e
        rt=time()-t0; support=classify_attempt(e); fail=failure_reason_of(e); status="error"
        p.rho==0 && support=="unstable_or_fails" && (support="invalid_parameter")
    end
    gap = (isnan(obj)||isnan(cen)||cen==0) ? NaN : 100*(obj-cen)/abs(cen)
    return Dict{Symbol,Any}(:case=>CASE,:k=>k,:formulation=>f,:algorithm=>algo,
        :sectioning=>sect,:coupling=>coup_lbl,:parameter_set=>plabel,:attempt_id=>aid,
        :rho=>p.rho,:status=>status,:objective=>obj,:centralized_objective=>cen,
        :gap_pct=>gap,:iterations=>iters,:primal_residual=>pr,:dual_residual=>dr,
        :runtime_s=>rt,:solver=>official_solver_label(f),:convergence_flag=>conv,
        :support_status=>support,:failure_reason=>fail,:timestamp=>string(Dates.now()))
end

# grids
admm_params()   = [(; rho=rho, adapt=ad) for rho in [0.0,1e2,1e3,1e4,1e5,1e6] for ad in (false,true)]
aladin_params() = [(; rho=rho, mu=mu, gr=1.05) for rho in [0.0,1e3,1e4] for mu in [1e3,1e4]]

# run a cell's grid with early-stop; return (rows, status)
function run_grid(k,f,algo,sect,coup,params,oracle; proximal_mode=nothing)
    rows=Dict{Symbol,Any}[]
    for (i,p) in enumerate(params)
        r = attempt(k,f,algo,sect,coup,p,i,oracle; proximal_mode=proximal_mode)
        push!(rows,r)
        r[:support_status] in ("not_implemented","skipped_by_design") && return rows, r[:support_status]
        if r[:convergence_flag] && !isnan(r[:gap_pct]) && abs(r[:gap_pct])<GAP_TARGET
            return rows, "converged"
        end
    end
    any(x->x[:support_status]=="supported", rows) && return rows, "failed"
    return rows, "not_implemented"
end

cell_class(a,l) =
    (a=="not_implemented" && l=="not_implemented") ? "not_implemented" :
    (a in ("skipped_by_design",) && l in ("skipped_by_design",)) ? "skipped_by_design" :
    (a=="converged" && l=="converged") ? "both_converged" :
    (a=="converged" && l!="converged") ? "admm_converged_aladin_failed" :
    (a!="converged" && l=="converged") ? "aladin_converged_admm_failed" :
    "both_failed"

paper_class(cls) =
    cls=="both_converged" ? "paper_ready" :
    cls=="admm_converged_aladin_failed" ? "needs_aladin_tuning" :
    cls=="aladin_converged_admm_failed" ? "needs_admm_tuning" :
    cls=="both_failed" ? "exclude_from_main" :
    "support_only"

function run_case14_gate()
    mkpath(RES_TUNING); mkpath(ARTIFACTS_DIR)
    @info "D1a case14 gate config" JULIA_NUM_THREADS=Threads.nthreads() BLAS_threads=BLAS.get_num_threads() ma57=have_hsl() gurobi=(_GRB!==nothing)
    oracle = compute_oracle()

    all_rows = Dict{Symbol,Any}[]; matrix = Dict{String,Any}[]
    for k in KS, f in FORMS, s in SECTS, c in COUPS
        arows, ast = run_grid(k,f,"admm",s,c,admm_params(),oracle)
        lrows, lst = run_grid(k,f,"aladin",s,c,aladin_params(),oracle)
        append!(all_rows, arows); append!(all_rows, lrows)
        cls = cell_class(ast,lst)
        push!(matrix, Dict("k"=>k,"formulation"=>f,"sectioning"=>s,"coupling"=>c,
            "ADMM_status"=>ast,"ALADIN_status"=>lst,"classification"=>cls,
            "paper_status"=>paper_class(cls)))
        @info @sprintf("  k%d %-6s %-10s %-6s -> ADMM=%-16s ALADIN=%-16s => %s",
            k,f,s,c,ast,lst,cls)
    end
    # pflow probes (ALADIN only; state+pflow is skipped_by_design)
    for k in KS, f in FORMS, s in SECTS, c in ["hybrid","power"]
        lrows, lst = run_grid(k,f,"aladin",s,c,aladin_params(),oracle; proximal_mode="pflow")
        append!(all_rows, lrows)
        push!(matrix, Dict("k"=>k,"formulation"=>f,"sectioning"=>s,"coupling"=>"$(c)+pflow",
            "ADMM_status"=>"n/a","ALADIN_status"=>lst,"classification"=>(lst=="converged" ? "aladin_only_pflow" : "pflow_"*lst),
            "paper_status"=>(lst=="converged" ? "paper_ready" : "support_only")))
    end

    cols=[:case,:k,:formulation,:algorithm,:sectioning,:coupling,:parameter_set,:attempt_id,
          :rho,:status,:objective,:centralized_objective,:gap_pct,:iterations,:primal_residual,
          :dual_residual,:runtime_s,:solver,:convergence_flag,:support_status,:failure_reason,:timestamp]
    nts = [(; (c=>r[c] for c in cols)...) for r in all_rows]
    write_csv(joinpath(RES_TUNING,"case14_gate_raw.csv"), nts; columns=cols)
    write_csv(joinpath(TABLES_DIR,"D1a_case14_gate.csv"), nts; columns=cols)
    write_markdown(joinpath(TABLES_DIR,"D1a_case14_gate.md"), nts; columns=cols, title="D1a - case14 gate (all attempts)")
    write_toml(joinpath(ARTIFACTS_DIR,"model_support_matrix.toml"), Dict(
        "_meta"=>Dict("generated_at"=>string(Dates.now()),"stage"=>"D1a","case"=>CASE,
            "centralized_models"=>["DC","AC","ACR","IVR","LPAC","SOCWR","QC","SDPWRM","SOCBF","BFA","NFA"],
            "distributed_models"=>["DC","AC","SOCWR"],
            "centralized_only_models"=>["ACR","IVR","LPAC","QC","SDPWRM","SOCBF","BFA","NFA"],
            "algorithms"=>["centralized","admm","aladin"],
            "sectionings"=>["tie_line","dummy_bus","split_bus"],
            "couplings"=>["state","hybrid","power","pflow(aladin proximal_mode)"]),
        "cell"=>matrix))
    @info "Wrote D1a artifacts" n_attempts=length(all_rows) n_cells=length(matrix)
    return all_rows, matrix
end

run_case14_gate()
