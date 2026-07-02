# ============================================================
# paper/scripts/experiments/run_tuning_smoke.jl
# ============================================================
# Stage D0 — tuning infrastructure smoke. Validates the declared-grid
# tuning mechanism on a tiny set of cells: reads grids from TOML, runs
# every attempt with robust failure capture, selects the best per cell
# by the declared lexicographic rule, and writes all attempts.
# NOT the full tuning, NOT the main matrix.
#
# Official runtime (C0): JULIA_NUM_THREADS=1 (env) + BLAS threads = 1.
# Official distributed solvers: DC=OSQP, SOCWR=Clarabel(hardened),
# AC=Ipopt+MA57 via LIRION_HSLLIB.
#
# Outputs:
#   results/tuning/tuning_smoke_raw.csv            (scratch, git-ignored)
#   artifacts/tables/D0_tuning_smoke.csv / .md     (curated)
#   artifacts/tuning/D0_tuning_smoke_summary.toml  (curated)
# ============================================================

include(joinpath(@__DIR__, "lib", "tuning_common.jl"))

using LinearAlgebra
BLAS.set_num_threads(1)

const RES_TUNING = joinpath(RESULTS_DIR, "tuning")

# Cells: (case, k, formulation, algorithm, sectioning, coupling)
const D0_CELLS = [
    ("case118_ieee", 2, "dc", "admm",   "tie_line", "hybrid"),  # MUST exercise failure capture (rho=1e3)
    ("case57_ieee",  2, "dc", "admm",   "tie_line", "hybrid"),  # should work
    ("case57_ieee",  2, "ac", "admm",   "tie_line", "hybrid"),  # AC + MA57 (official)
    ("case57_ieee",  2, "dc", "aladin", "tie_line", "hybrid"),  # ALADIN small
    ("case57_ieee",  2, "ac", "aladin", "tie_line", "hybrid"),  # ALADIN AC
]

function run_tuning_smoke()
    mkpath(RES_TUNING); mkpath(joinpath(ARTIFACTS_DIR, "tuning"))
    t1 = load_t1_objectives()
    admm_grid = load_grid("tuning_grid_admm_smoke.toml")
    aladin_grid = load_grid("tuning_grid_aladin_smoke.toml")
    crit = load_criteria("selection_criteria_smoke.toml")
    gap_target = crit["accuracy"]["gap_target_pct"]

    @info "D0 tuning smoke config" JULIA_NUM_THREADS = Threads.nthreads() BLAS_threads = BLAS.get_num_threads() n_cells = length(D0_CELLS) ma57 = have_hsl() gap_target n_t1 = length(t1)

    all_rows = Dict{Symbol,Any}[]
    per_cell_selected = Dict{String,Any}[]

    for cell in D0_CELLS
        (case, k, f, algo, sect, coup) = cell
        grid = algo == "admm" ? admm_grid : aladin_grid
        psets = param_sets(algo, grid)
        cell_rows = Dict{Symbol,Any}[]
        for (i, p) in enumerate(psets)
            row = run_attempt(cell, p, i; t1 = t1)
            push!(cell_rows, row)
            @info @sprintf("  %-14s k%d %-6s %-6s %-22s -> %-16s conv=%s gap=%s it=%d %s",
                case, k, f, algo, row[:parameter_set], row[:support_status],
                row[:convergence_flag], isnan(row[:gap_pct]) ? "NaN" : @sprintf("%.3f%%", row[:gap_pct]),
                row[:iterations], isempty(row[:failure_reason]) ? "" : "FAIL=" * row[:failure_reason][1:min(end,60)])
        end
        select_best!(cell_rows, gap_target)
        sel = cell_rows[findfirst(r -> r[:is_selected], cell_rows)]
        push!(per_cell_selected, Dict(
            "cell" => "$(case)/k$(k)/$(f)/$(algo)/$(sect)/$(coup)",
            "selected_parameter_set" => sel[:parameter_set],
            "convergence_flag" => sel[:convergence_flag],
            "gap_pct" => isnan(sel[:gap_pct]) ? "NaN" : sel[:gap_pct],
            "support_status" => sel[:support_status],
            "failure_reason" => sel[:failure_reason],
            "n_attempts" => length(cell_rows),
            "n_failed" => count(r -> r[:support_status] != "supported", cell_rows),
        ))
        append!(all_rows, cell_rows)
    end

    nts = as_nt.(all_rows)
    raw = write_csv(joinpath(RES_TUNING, "tuning_smoke_raw.csv"), nts; columns = TUNING_COLUMNS)
    csv = write_csv(joinpath(TABLES_DIR, "D0_tuning_smoke.csv"), nts; columns = TUNING_COLUMNS)
    md  = write_markdown(joinpath(TABLES_DIR, "D0_tuning_smoke.md"), nts;
                         columns = TUNING_COLUMNS, title = "D0 — tuning infrastructure smoke")
    write_toml(joinpath(ARTIFACTS_DIR, "tuning", "D0_tuning_smoke_summary.toml"), Dict(
        "_meta" => Dict("generated_at" => string(Dates.now()), "stage" => "D0",
            "note" => "tuning-infrastructure smoke; validates declared grid + robust " *
                      "failure capture + lexicographic selection. Not final results.",
            "n_attempts" => length(all_rows),
            "runtime" => Dict("JULIA_NUM_THREADS" => 1, "BLAS_threads" => 1)),
        "selected" => per_cell_selected,
    ))
    @info "Wrote" raw = relpath(raw, PAPER_DIR) csv = relpath(csv, PAPER_DIR) md = relpath(md, PAPER_DIR) n_attempts = length(all_rows)
    return all_rows
end

run_tuning_smoke()
