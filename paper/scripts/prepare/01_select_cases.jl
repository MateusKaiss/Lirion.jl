# ============================================================
# paper/scripts/prepare/01_select_cases.jl
# ============================================================
# Prepare step 1: resolve the selected PGLib cases, COPY every available
# variant (std/api/sad) into paper/cases/raw/pglib/, copy the PGLib
# LICENSE/README (CC-BY 4.0 attribution), load each case to record real
# bus/branch/gen counts + a topology signature, and emit
# results/cases/case_manifest.toml.
#
# Never touches the original PGLib files (copy only). No solver.
# Reproducible via: include(".../prepare/01_select_cases.jl")
# ============================================================

include(joinpath(@__DIR__, "..", "lib", "common.jl"))

function copy_pglib_license(root)
    for (src, dst) in (("LICENSE", "PGLIB_LICENSE"), ("README.md", "PGLIB_README.md"))
        s = joinpath(root, src)
        isfile(s) && cp(s, joinpath(RAW_DIR, dst); force = true)
    end
end

function run_select_cases()
    mkpath(LOG_DIR); mkpath(RAW_DIR)
    io = open(joinpath(LOG_DIR, "01_select_cases.log"), "w")
    try
        root = pglib_root()
        log_line(io, "PGLib root: $root")
        copy_pglib_license(root)
        log_line(io, "Copied PGLib LICENSE + README into cases/raw/pglib/")

        cases = load_selected_cases()["case"]
        entries = Vector{Dict{String,Any}}()

        for case in cases
            name = case["name"]
            for variant in VARIANTS
                src = pglib_variant_path(root, name, variant)
                if !isfile(src)
                    log_line(io, "  [skip] $name/$variant : not found")
                    continue
                end
                dstdir = joinpath(RAW_DIR, variant); mkpath(dstdir)
                dst = joinpath(dstdir, basename(src))
                cp(src, dst; force = true)

                c = load_inputs(src)
                sig = topology_signature(c)
                push!(entries, Dict{String,Any}(
                    "name" => name, "variant" => variant,
                    "family" => get(case, "family", ""),
                    "role" => get(case, "role", ""),
                    "group" => get(case, "group", ""),
                    "path_original" => src,
                    "path_copied" => relpath(dst, PAPER_DIR),
                    "n_bus" => n_bus(c), "n_branch" => n_branch(c), "n_gen" => n_gen(c),
                    "topology_hash" => string(sig.hash),
                ))
                log_line(io, @sprintf("  [ok]   %-16s %-3s  nbus=%-5d nbranch=%-5d ngen=%-4d",
                    name, variant, n_bus(c), n_branch(c), n_gen(c)))
            end
        end

        manifest = Dict{String,Any}(
            "_meta" => Dict("generated_at" => string(now()),
                            "script" => "prepare/01_select_cases.jl",
                            "pglib_root" => root, "n_entries" => length(entries)),
            "case" => entries,
        )
        out = joinpath(RESULTS_DIR, "cases", "case_manifest.toml")
        write_toml(out, manifest)
        log_line(io, "Wrote manifest: $(relpath(out, PAPER_DIR))  ($(length(entries)) entries)")
    finally
        close(io)
    end
    return nothing
end

run_select_cases()
