# ============================================================
# _apply_aladin_v2.jl
# ============================================================
# Reads `_tune_case14_aladin_v2.toml` and rewrites the
# `ALADIN_OVERRIDES` block in `paper_results_aladin.jl`. Each winner
# entry records the full strategy kwargs (ρ, μ, rho_growth, mu_growth,
# adaptive) so paper reproduction uses exactly the tuned config.
# ============================================================

using Printf
using TOML

const PAPER_DIR = @__DIR__
const CASE = "case14_ieee"
const TOML_IN = joinpath(PAPER_DIR, "cases", "_tune_case14_aladin_v2.toml")
const PAPER_FILE = joinpath(PAPER_DIR, "paper_results_aladin.jl")
const MARKER = "const ALADIN_OVERRIDES = Dict"

_variants() = ("std", "api", "sad")
_bases()    = ("dc", "ac", "socwr")
_combos()   = (
    ("tie_line",  "state"),  ("tie_line",  "power"),  ("tie_line",  "hybrid"),
    ("dummy_bus", "power"),  ("dummy_bus", "hybrid"), ("split_bus", "hybrid"),
)

function _format_kwarg(k::Symbol, v)
    if v isa Bool
        return @sprintf(":%s => %s", k, v)
    elseif v isa AbstractFloat
        return @sprintf(":%s => %.3g", k, v)
    else
        return @sprintf(":%s => %s", k, v)
    end
end

function _format_block(toml::Dict)
    io = IOBuffer()
    println(io, "const ALADIN_OVERRIDES = Dict{Tuple{String, String, Symbol, Symbol, Symbol}, Dict{Symbol, Any}}(")
    println(io, "    # --- case14_ieee (from _tune_case14_aladin_v2.toml) --------------")
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
    n_hit == 0 && println(io, "    # (no combos hit target)")
    println(io, ")")
    return String(take!(io)), n_hit
end

function _splice!(file::String, block::String)
    src = read(file, String)
    i = findfirst(MARKER, src)
    i === nothing && error("marker not found in $file")
    start = first(i)
    rest = SubString(src, start)
    j = findfirst("\n)\n", rest)
    j === nothing && error("end marker not found")
    stop = start + last(j) - 1
    new_src = src[1:start-1] * block * src[stop+1:end]
    write(file, new_src)
end

function main()
    isfile(TOML_IN) || error("missing $TOML_IN")
    toml = TOML.parsefile(TOML_IN)
    block, n_hit = _format_block(toml)
    _splice!(PAPER_FILE, block)
    println("ALADIN v2: wrote $n_hit overrides into $PAPER_FILE")
end

main()
