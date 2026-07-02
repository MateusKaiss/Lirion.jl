# ============================================================
# paper/scripts/experiments/run_distributed_smoke.jl
# ============================================================
# Stage B — distributed smoke test. Validates the distributed pipeline on
# a tiny controlled set: reads partitioned cases, applies sectioning at
# runtime, runs ADMM & ALADIN in one schema, computes gap vs the T1
# centralized oracle, and classifies each cell into a support matrix.
#
# NOT tuning. Small declared limits. Solver per formulation is Lirion's
# own tuned distributed default (optimizer = nothing):
#   DC -> OSQP,  SOCWR -> Clarabel (hardened),  AC -> Ipopt/MUMPS
# (reported explicitly, not silently changed).
#
# Outputs:
#   results/distributed/distributed_smoke_raw.csv   (scratch, git-ignored)
#   artifacts/tables/B_distributed_smoke.csv         (curated)
#   artifacts/tables/B_distributed_smoke.md          (curated)
#   artifacts/support_matrix.toml                    (curated)
#
# Reproducible via: include(".../experiments/run_distributed_smoke.jl")
# ============================================================

include(joinpath(@__DIR__, "lib", "distributed_common.jl"))

using LinearAlgebra
BLAS.set_num_threads(1)   # deterministic (JULIA_NUM_THREADS=1 in env)

const SMOKE_ADMM_MAXITER   = 300
const SMOKE_ALADIN_MAXITER = 40

# Small controlled grid (std only).
const B_CASES       = [("case57_ieee", 2), ("case118_ieee", 2)]
const B_SECTIONINGS = ["tie_line", "dummy_bus", "split_bus"]
const B_COUPLINGS   = ["state", "hybrid"]
const B_FORMS       = ["dc", "socwr", "ac"]
const B_ALGOS       = ["admm", "aladin"]

_run(case, k, f, a, s, c; proximal_mode = nothing, t1) =
    run_distributed_one(case, "std", k, f, a, s, c;
        proximal_mode = proximal_mode, t1 = t1,
        admm_max_iter = SMOKE_ADMM_MAXITER, aladin_max_iter = SMOKE_ALADIN_MAXITER)

function warmup(t1)
    for a in B_ALGOS, f in B_FORMS
        try
            _run("case57_ieee", 2, f, a, "tie_line", "state"; t1 = t1)
        catch
        end
    end
    @info "Warm-up done (discarded)"
end

function build_support_matrix(rows)
    # precedence: usable (supported) > unstable_or_fails > not_implemented > skipped
    rank = Dict("supported" => 4, "unstable_or_fails" => 3,
                "not_implemented" => 2, "skipped_by_design" => 1)
    best = Dict{NTuple{4,String},String}()
    for r in rows
        key = (r.formulation, r.algorithm, r.sectioning, r.coupling)
        s = r.support_status
        if !haskey(best, key) || rank[s] > rank[best[key]]
            best[key] = s
        end
    end
    cells = [Dict("formulation" => k[1], "algorithm" => k[2],
                  "sectioning" => k[3], "coupling" => k[4], "status" => v)
             for (k, v) in best]
    counts = Dict{String,Int}()
    for (_, v) in best; counts[v] = get(counts, v, 0) + 1; end
    return cells, counts
end

function run_distributed_smoke()
    mkpath(RES_DIST)
    t1 = load_t1_objectives()
    @info "Distributed smoke config" cases = B_CASES sectionings = B_SECTIONINGS couplings = B_COUPLINGS forms = B_FORMS algos = B_ALGOS JULIA_NUM_THREADS = Threads.nthreads() BLAS_threads = BLAS.get_num_threads() n_t1 = length(t1)

    warmup(t1)

    rows = NamedTuple[]
    # main grid
    for (case, k) in B_CASES, s in B_SECTIONINGS, c in B_COUPLINGS, f in B_FORMS, a in B_ALGOS
        row = _run(case, k, f, a, s, c; t1 = t1)
        push!(rows, row)
        @info @sprintf("  %-14s k%d %-6s %-6s %-10s %-6s -> %-16s obj=%s gap=%s it=%d",
            case, k, f, a, s, c, row.support_status,
            isnan(row.objective) ? "NaN" : @sprintf("%.2f", row.objective),
            isnan(row.gap_pct) ? "NaN" : @sprintf("%.3f%%", row.gap_pct), row.iterations)
    end

    # probes: register power + pflow support on case57 k2, tie_line only
    for f in B_FORMS, a in B_ALGOS
        push!(rows, _run("case57_ieee", 2, f, a, "tie_line", "power"; t1 = t1))
    end
    for f in B_FORMS, c in ["state", "power", "hybrid"]   # pflow is ALADIN-only
        push!(rows, _run("case57_ieee", 2, f, "aladin", "tie_line", c;
                         proximal_mode = "pflow", t1 = t1))
    end

    # write results
    raw = write_csv(joinpath(RES_DIST, "distributed_smoke_raw.csv"), rows; columns = DIST_COLUMNS)
    csv = write_csv(joinpath(TABLES_DIR, "B_distributed_smoke.csv"), rows; columns = DIST_COLUMNS)
    md  = write_markdown(joinpath(TABLES_DIR, "B_distributed_smoke.md"), rows;
                         columns = DIST_COLUMNS, title = "B — distributed smoke")

    cells, counts = build_support_matrix(rows)
    write_toml(joinpath(ARTIFACTS_DIR, "support_matrix.toml"), Dict(
        "_meta" => Dict("generated_at" => string(Dates.now()),
                        "script" => "run_distributed_smoke.jl",
                        "note" => "smoke classification; cells over formulation×algorithm×sectioning×coupling",
                        "counts" => counts),
        "cell" => cells,
    ))
    @info "Wrote" raw = relpath(raw, PAPER_DIR) csv = relpath(csv, PAPER_DIR) md = relpath(md, PAPER_DIR) counts
    return rows
end

run_distributed_smoke()
