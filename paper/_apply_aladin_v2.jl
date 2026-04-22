# ============================================================
# _apply_aladin_v2.jl
# ============================================================
# Reads `_tune_<CASE>_aladin_v2.toml` and merges the winner entries
# into the `ALADIN_OVERRIDES` block in `paper_results_aladin.jl`.
# Preserves entries for any other case already in the block.
# Set `LIRION_TUNE_CASE` to target a case other than case14_ieee.
# ============================================================

using Printf
using TOML

const PAPER_DIR  = @__DIR__
const CASE       = get(ENV, "LIRION_TUNE_CASE", "case14_ieee")
const TOML_IN    = joinpath(PAPER_DIR, "cases", "_tune_$(CASE)_aladin_v2.toml")
const PAPER_FILE = joinpath(PAPER_DIR, "paper_results_aladin.jl")
const MARKER     = "const ALADIN_OVERRIDES = Dict"
const END_MARKER = "\n)\n"

_variants() = ("std", "api", "sad")
_bases()    = ("dc", "ac", "socwr")
_combos()   = (
    ("tie_line",  "state"),  ("tie_line",  "power"),  ("tie_line",  "hybrid"),
    ("dummy_bus", "power"),  ("dummy_bus", "hybrid"), ("split_bus", "hybrid"),
)

function _format_case_lines(toml::Dict)
    io = IOBuffer()
    println(io, "    # --- $(CASE) (from _tune_$(CASE)_aladin_v2.toml) --------------")
    n_hit = 0
    for v in _variants(), b in _bases()
        haskey(toml, v) && haskey(toml[v], b) || continue
        for (sect, coup) in _combos()
            k = "$(sect)_$(coup)"
            haskey(toml[v][b], k) || continue
            cell = toml[v][b][k]
            get(cell, "hit_target", false) || continue
            w = cell["winner"]
            n_hit += 1
            parts = String[]
            push!(parts, @sprintf(":rho => %.3g", Float64(w["rho"])))
            push!(parts, @sprintf(":mu => %.3g", Float64(w["mu"])))
            if get(w, "adaptive", false)
                push!(parts, ":adaptive => true")
                isfinite(get(w, "rho_growth", NaN)) &&
                    push!(parts, @sprintf(":rho_growth => %.3g", Float64(w["rho_growth"])))
                isfinite(get(w, "mu_growth", NaN)) &&
                    push!(parts, @sprintf(":mu_growth => %.3g", Float64(w["mu_growth"])))
            else
                push!(parts, ":adaptive => false")
            end
            @printf(io, "    (\"%s\", \"%s\", :%s, :%s, :%s) => Dict(%s),  # %s gap=%+.3f%%, iters=%d\n",
                    CASE, v, b, sect, coup, join(parts, ", "),
                    w["strategy"], Float64(w["gap_pct"]), Int(w["iters"]))
        end
    end
    return String(take!(io)), n_hit
end

function _merge!(file::String, new_case_lines::String)
    src = read(file, String)
    i = findfirst(MARKER, src)
    i === nothing && error("marker not found in $file")
    start = first(i)
    rest = SubString(src, start)
    j = findfirst(END_MARKER, rest)
    j === nothing && error("end marker not found")
    stop = start + last(j) - 1

    block = String(SubString(src, start, stop))
    body_lines = split(block, '\n')
    header = String(body_lines[1])

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
    print(io, new_case_lines)
    println(io, ")")
    new_block = String(take!(io))

    new_src = src[1:start-1] * new_block * src[stop+2:end]
    write(file, new_src)
end

function main()
    isfile(TOML_IN) || error("missing $TOML_IN")
    toml = TOML.parsefile(TOML_IN)
    lines, n_hit = _format_case_lines(toml)
    _merge!(PAPER_FILE, lines)
    println("ALADIN v2: wrote $n_hit $(CASE) overrides into $PAPER_FILE")
end

main()
