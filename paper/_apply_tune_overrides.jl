# ============================================================
# _apply_tune_overrides.jl
# ============================================================
# Reads `_tune_case14_admm.toml` + `_tune_case14_aladin.toml`,
# selects the best (ρ, adaptive) per (variant, base, sectioning,
# coupling), and rewrites the `ADMM_OVERRIDES` / `ALADIN_OVERRIDES`
# blocks in the corresponding paper_results_*.jl files for case14.
#
# Selection rule: among runs that (a) converged and (b) satisfy
# gap < GAP_TARGET, pick the one with the fewest iterations and
# break ties by smallest |gap|. Prefer fixed ρ if both fixed and
# adaptive satisfy the target (simpler, more reproducible).
#
# Usage:
#     julia> include(joinpath(pkgdir(Lirion), "paper", "_apply_tune_overrides.jl"))
# ============================================================

using Printf
using TOML

const PAPER_DIR  = @__DIR__
const CASES_DIR  = joinpath(PAPER_DIR, "cases")
const CASE       = get(ENV, "LIRION_TUNE_CASE", "case14_ieee")
const GAP_TARGET = 1.0

_variants() = ("std", "api", "sad")
_bases()    = ("dc", "ac", "socwr")
_combos()   = (
    ("tie_line",  "state"),
    ("tie_line",  "power"),
    ("tie_line",  "hybrid"),
    ("dummy_bus", "power"),
    ("dummy_bus", "hybrid"),
    ("split_bus", "hybrid"),
)

# Pick a best sweep entry. Returns nothing if none converged w/ gap<target.
function _pick_best(sweep::Dict)
    best = nothing
    for (rkey, rec) in sweep
        rec isa Dict || continue
        get(rec, "status", "") == "converged" || continue
        gap = Float64(get(rec, "gap_pct", NaN))
        isfinite(gap) && abs(gap) < GAP_TARGET || continue
        iters = Int(get(rec, "iters", typemax(Int)))
        if best === nothing ||
           iters < best.iters ||
           (iters == best.iters && abs(gap) < abs(best.gap))
            best = (rho = Float64(get(rec, "rho", NaN)),
                    gap = gap, iters = iters,
                    time_s = Float64(get(rec, "time_s", 0.0)))
        end
    end
    return best
end

# For a combo: compare fixed and adaptive, prefer fixed when both hit target.
function _pick_combo(entry::Dict)
    fixed_sweep    = get(get(entry, "fixed",    Dict()), "sweep", Dict())
    adaptive_sweep = get(get(entry, "adaptive", Dict()), "sweep", Dict())
    bf = _pick_best(fixed_sweep)
    ba = _pick_best(adaptive_sweep)
    if bf !== nothing
        return (adaptive = false, bf...)
    elseif ba !== nothing
        return (adaptive = true,  ba...)
    else
        return nothing
    end
end

function _collect_admm(toml::Dict)
    rows = Tuple[]
    for v in _variants(), b in _bases()
        haskey(toml, v) || continue
        haskey(toml[v], b) || continue
        for (sect, coup) in _combos()
            key = "$(sect)_$(coup)"
            haskey(toml[v][b], key) || continue
            pick = _pick_combo(toml[v][b][key])
            pick === nothing && (push!(rows, (v, b, sect, coup, nothing)); continue)
            push!(rows, (v, b, sect, coup, pick))
        end
    end
    return rows
end

function _format_admm_lines(rows::Vector)
    io = IOBuffer()
    println(io, "    # --- $(CASE) (auto-generated from _tune_$(CASE)_admm.toml) -----")
    for (v, b, sect, coup, pick) in rows
        pick === nothing && continue
        @printf(io, "    (\"%s\", \"%s\", :%s, :%s, :%s) => Dict(:rho => %.0e",
                CASE, v, b, sect, coup, pick.rho)
        if pick.adaptive
            print(io, ", :adaptive_rho => true")
        end
        @printf(io, "),  # gap=%+6.3f%%, iters=%d\n", pick.gap, pick.iters)
    end
    return String(take!(io))
end

function _format_aladin_lines(rows::Vector)
    io = IOBuffer()
    println(io, "    # --- $(CASE) (auto-generated from _tune_$(CASE)_aladin.toml) ----")
    for (v, b, sect, coup, pick) in rows
        pick === nothing && continue
        mu = pick.rho / 10
        @printf(io, "    (\"%s\", \"%s\", :%s, :%s, :%s) => Dict(:rho => %.0e, :mu => %.0e",
                CASE, v, b, sect, coup, pick.rho, mu)
        if pick.adaptive
            print(io, ", :adaptive => true")
        else
            print(io, ", :adaptive => false")
        end
        @printf(io, "),  # gap=%+6.3f%%, iters=%d\n", pick.gap, pick.iters)
    end
    return String(take!(io))
end

const _ADMM_MARKER_START = "const ADMM_OVERRIDES = Dict"
const _ALAD_MARKER_START = "const ALADIN_OVERRIDES = Dict"
const _END_MARKER        = "\n)\n"

# Merge new case-specific lines into an existing Dict block, preserving
# every entry for OTHER cases. The caller passes a line-generator that
# returns just the lines for `CASE`; this function rebuilds the block
# as: marker header + existing lines for other cases + new lines for
# CASE + closing paren.
function _merge_block!(path::String, marker::String, new_case_lines::String)
    src = read(path, String)
    i = findfirst(marker, src)
    i === nothing && error("marker $(repr(marker)) not found in $path")
    start = first(i)
    rest = SubString(src, start)
    j = findfirst(_END_MARKER, rest)
    j === nothing && error("end marker not found after $(repr(marker)) in $path")
    stop = start + last(j) - 1

    # Collect existing block body (between first "(" after marker and final ")")
    block = String(SubString(src, start, stop))
    body_lines = split(block, '\n')
    header = body_lines[1]  # `const ... Dict{...}(`

    # Keep every existing line whose tuple key does not mention CASE and
    # whose comment banner does not mention CASE. Drop the block's own
    # closing paren (we re-add it after splicing case lines).
    case_prefix = "(\"$(CASE)\","
    case_banner = "--- $(CASE) "
    kept = String[]
    for line in body_lines[2:end]
        s = String(line)
        stripped = strip(s)
        stripped == ")" && continue
        occursin(case_prefix, s) && continue
        occursin(case_banner, s) && continue
        push!(kept, s)
    end

    io = IOBuffer()
    println(io, header)
    for ln in kept
        isempty(strip(ln)) && continue
        println(io, ln)
    end
    print(io, new_case_lines)  # already ends with newline per _format_*
    println(io, ")")
    new_block = String(take!(io))

    new_src = src[1:start-1] * new_block * src[stop+2:end]
    write(path, new_src)
    return nothing
end

function main()
    out_admm   = joinpath(CASES_DIR, "_tune_$(CASE)_admm.toml")
    out_aladin = joinpath(CASES_DIR, "_tune_$(CASE)_aladin.toml")
    paper_admm = joinpath(PAPER_DIR, "paper_results_admm.jl")
    paper_aladin = joinpath(PAPER_DIR, "paper_results_aladin.jl")

    if isfile(out_admm)
        toml = TOML.parsefile(out_admm)
        rows = _collect_admm(toml)
        lines = _format_admm_lines(rows)
        _merge_block!(paper_admm, _ADMM_MARKER_START, lines)
        n_ok = count(r -> r[5] !== nothing, rows)
        println("ADMM:   wrote ", n_ok, "/", length(rows), " $(CASE) overrides into ", paper_admm)
    else
        @warn "No ADMM tune file at $out_admm — skipping."
    end

    if isfile(out_aladin)
        toml = TOML.parsefile(out_aladin)
        rows = _collect_admm(toml)  # same shape
        lines = _format_aladin_lines(rows)
        _merge_block!(paper_aladin, _ALAD_MARKER_START, lines)
        n_ok = count(r -> r[5] !== nothing, rows)
        println("ALADIN: wrote ", n_ok, "/", length(rows), " $(CASE) overrides into ", paper_aladin)
    else
        @warn "No ALADIN tune file at $out_aladin — skipping."
    end
end

main()
