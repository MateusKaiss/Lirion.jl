# ============================================================
# paper/scripts/prepare/02_build_partitions.jl
# ============================================================
# Prepare step 2: build the FINAL partitioned cases used by the paper.
#
# The Lirion solver applies sectioning (TieLine/DummyBus/SplitBus) at
# RUNTIME from full.bus.zone (prepare_distributed_partitions). Therefore
# a partition is ONE case.m carrying only the `zone` column — NO
# per-sectioning files are produced. The sectionings to be tested later
# are recorded as metadata (sectionings_planned) in meta.toml.
#
# Per (case, k) we search seeds x refine (on the reference variant std)
# and keep the cleanest partition by the declared criterion:
#   1. all zones connected   2. all zones have generation
#   3. lowest imbalance      4. fewest tie-lines
#   5. refine=true / smallest seed
# The winning bus->zone map is shared across std/api/sad (identical
# topology) and written to:
#   cases/partitioned/<case>/<variant>/k<k>/case.m
#   cases/partitioned/<case>/<variant>/k<k>/meta.toml
# Roll-ups: cases/partitioned/{_summary,_diagnostic}.toml
#
# Requires Metis/Graphs/SimpleWeightedGraphs. Runs NO solver.
# Reproducible via: include(".../prepare/02_build_partitions.jl")
# ============================================================

include(joinpath(@__DIR__, "..", "lib", "common.jl"))

using Metis, Graphs, SimpleWeightedGraphs   # activate Lirion MetisExt

const SECTIONINGS_PLANNED = ["tie_line", "dummy_bus", "split_bus"]

function eval_partition(c::Lirion.CaseRaw)
    members = zone_members(c)
    ties = tie_lines(c)
    bal = balance_metrics(members)
    conn = zone_connectivity(c)
    gen = zone_generation(c)
    return (
        k_actual = length(members),
        n_tie_lines = length(ties),
        n_boundary = length(boundary_buses(ties)),
        cv = bal["cv"],
        imbalance = bal["imbalance_max_over_mean"],
        all_connected = conn.all_connected,
        all_gen = gen.all_zones_have_gen,
    )
end

# Lower is better — encodes the declared priority order.
selection_key(a) = (
    a.all_connected ? 0 : 1,
    a.all_gen ? 0 : 1,
    a.cv,
    a.n_tie_lines,
    a.refine ? 0 : 1,
    a.seed,
)

function apply_zone_map!(c::Lirion.CaseRaw, zmap::Dict{Int,Int})
    for i in eachindex(c.bus.bus_i)
        c.bus.zone[i] = zmap[c.bus.bus_i[i]]
    end
    return c
end

status_of(a) = !a.all_connected ? "not_clean" : (a.all_gen ? "clean" : "connected_no_gen")

function write_variant!(name, variant, k, chosen, src, c::Lirion.CaseRaw)
    m = eval_partition(c)
    members = zone_members(c)
    ties = tie_lines(c)
    conn = zone_connectivity(c)
    gen = zone_generation(c)
    warnings = String[]
    m.all_connected || push!(warnings, "not all zones connected")
    m.all_gen || push!(warnings, "some zone has no generator")
    m.k_actual == k || push!(warnings, "k_actual=$(m.k_actual) != k_requested=$k")

    outdir = joinpath(PARTITIONED_DIR, name, variant, "k$(k)")
    mkpath(outdir)
    mpath = joinpath(outdir, "case.m")
    write_matpower(c, mpath)

    meta = Dict{String,Any}(
        "case" => name, "variant" => variant,
        "k_requested" => k, "k_actual" => m.k_actual,
        "seed" => chosen.seed, "refine" => chosen.refine, "method" => "metis",
        "status" => status_of(m),
        "clean_for_main" => (m.all_connected && m.all_gen),
        # connected but a zone lacks a dispatchable generator (imports via
        # tie-lines): valid, but kept out of the headline matrix.
        "optional_or_sensitivity" => (status_of(m) == "connected_no_gen"),
        "shared_partition" => true, "shared_with_reference" => "std",
        "path_original" => src,
        "path_copied" => relpath(joinpath(RAW_DIR, variant, basename(src)), PAPER_DIR),
        "path_partitioned" => relpath(mpath, PAPER_DIR),
        "n_bus" => n_bus(c), "n_branch" => n_branch(c), "n_gen" => n_gen(c),
        "n_tie_lines" => m.n_tie_lines, "n_boundary_buses" => m.n_boundary,
        "bus_to_zone" => Dict(string(b) => z for (b, z) in bus_zone_map(c)),
        "tie_lines" => ties, "boundary_buses" => boundary_buses(ties),
        "zones" => Dict(string(z) => v for (z, v) in members),
        "balance" => balance_metrics(members),
        "connectivity" => Dict("all_connected" => conn.all_connected,
                               "per_zone" => conn.per_zone),
        "generation" => Dict("all_zones_have_gen" => gen.all_zones_have_gen,
                             "per_zone" => gen.per_zone),
        # Planning only — sectioning is applied at RUNTIME, not materialized.
        "sectionings_planned" => SECTIONINGS_PLANNED,
        "sectioning_coupling_counts" => Dict(
            "tie_line"  => Dict("coupling_edges" => m.n_tie_lines),
            "dummy_bus" => Dict("coupling_edges" => m.n_tie_lines),
            "split_bus" => Dict("split_buses" => m.n_boundary),
        ),
        "warnings" => warnings,
    )
    write_toml(joinpath(outdir, "meta.toml"), meta)
    return m, status_of(m)
end

function run_build_partitions()
    mkpath(LOG_DIR)
    io = open(joinpath(LOG_DIR, "02_build_partitions.log"), "w")
    summary = Vector{Dict{String,Any}}()
    diagnostic = Vector{Dict{String,Any}}()
    try
        root = pglib_root()
        pcfg = load_partitioning_config()
        settings = pcfg["settings"]
        kmap = pcfg["k"]
        seeds = Int.(settings["seeds"])
        refine_opts = Bool.(settings["refine"])
        max_passes = Int(settings["max_refine_passes"])
        cases = load_selected_cases()["case"]
        log_line(io, "seeds=$seeds refine=$refine_opts")

        for case in cases
            name = case["name"]
            haskey(kmap, name) || (log_line(io, "[skip] $name : no k"); continue)
            ks = Int.(kmap[name])

            avail = [v for v in VARIANTS if isfile(pglib_variant_path(root, name, v))]
            ref = "std" in avail ? "std" : avail[1]
            refcase0 = load_inputs(pglib_variant_path(root, name, ref))

            for k in ks
                attempts = NamedTuple[]
                for refine in refine_opts, seed in seeds
                    c = deepcopy(refcase0)
                    try
                        create_partition(c; k = k, seed = seed, refine = refine,
                                          max_refine_passes = max_passes)
                    catch e
                        log_line(io, "    [warn] $name k=$k seed=$seed refine=$refine failed: " *
                                 sprint(showerror, e)[1:min(end,120)])
                        continue
                    end
                    m = eval_partition(c)
                    push!(attempts, (seed = seed, refine = refine, m...))
                    push!(diagnostic, Dict(
                        "case" => name, "k" => k, "seed" => seed, "refine" => refine,
                        "k_actual" => m.k_actual, "all_connected" => m.all_connected,
                        "all_zones_have_gen" => m.all_gen, "n_tie_lines" => m.n_tie_lines,
                        "cv" => m.cv, "imbalance" => m.imbalance,
                    ))
                end
                isempty(attempts) && (log_line(io, "    [FAIL] $name k=$k : no attempt"); continue)

                chosen = attempts[argmin(map(selection_key, attempts))]
                n_clean = count(a -> a.all_connected && a.all_gen, attempts)
                n_conn  = count(a -> a.all_connected, attempts)

                cref = deepcopy(refcase0)
                create_partition(cref; k = k, seed = chosen.seed,
                                 refine = chosen.refine, max_refine_passes = max_passes)
                zmap = bus_zone_map(cref)

                st = "?"
                for v in avail
                    src = pglib_variant_path(root, name, v)
                    cv = load_inputs(src)
                    apply_zone_map!(cv, zmap)
                    m, st = write_variant!(name, v, k, chosen, src, cv)
                    push!(summary, Dict(
                        "case" => name, "variant" => v, "k" => k, "k_actual" => m.k_actual,
                        "seed" => chosen.seed, "refine" => chosen.refine, "status" => st,
                        "clean_for_main" => (m.all_connected && m.all_gen),
                        "optional_or_sensitivity" => (st == "connected_no_gen"),
                        "all_connected" => m.all_connected, "all_zones_have_gen" => m.all_gen,
                        "n_tie_lines" => m.n_tie_lines, "cv" => m.cv, "imbalance" => m.imbalance,
                    ))
                end
                log_line(io, @sprintf("  %-16s k=%-2d -> seed=%-3d refine=%-5s %-16s ties=%-4d cv=%.3f  [%d/%d conn, %d clean]",
                    name, k, chosen.seed, string(chosen.refine), st,
                    chosen.n_tie_lines, chosen.cv, n_conn, length(attempts), n_clean))
            end
        end

        # Status counts over distinct (case, k) combos (variants share status).
        combo = Dict{Tuple{String,Int},String}()
        for p in summary
            combo[(p["case"], p["k"])] = p["status"]
        end
        counts = Dict("clean" => 0, "connected_no_gen" => 0, "not_clean" => 0)
        for (_, st) in combo
            counts[st] += 1
        end
        write_toml(joinpath(PARTITIONED_DIR, "_summary.toml"), Dict(
            "_meta" => Dict("generated_at" => string(now()),
                            "script" => "prepare/02_build_partitions.jl",
                            "n_partitions" => length(summary),
                            "n_case_k" => length(combo),
                            "status_counts_by_case_k" => counts,
                            "main_matrix" => "clean only (clean_for_main = true)",
                            "note" => "connected_no_gen => optional/sensitivity; not_clean should be 0"),
            "partition" => summary,
        ))
        write_toml(joinpath(PARTITIONED_DIR, "_diagnostic.toml"), Dict(
            "_meta" => Dict("generated_at" => string(now()),
                            "note" => "every seed x refine attempt on the reference (std) variant"),
            "attempt" => diagnostic,
        ))
        log_line(io, "Chosen $(length(summary)) partitions; summary + diagnostic in cases/partitioned/")
    finally
        close(io)
    end
    return nothing
end

run_build_partitions()
