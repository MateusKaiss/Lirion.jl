# Compares a regenerated crossing table with cases/<case>/expected_crossing.csv.
# julia scripts/verify_crossing.jl <crossing_table.csv> <expected.csv>...
# Exit 0 when every expected row reproduces, 1 otherwise.

const KEY_COLS = ["case", "k", "formulation", "sectioning", "coupling", "algorithm"]
const CHECK_COLS = ["tier1_first_touch", "tier2_iter", "tier2_first_touch_raw", "verdict"]

# minimal RFC-4180 reader: quoted fields may contain commas and doubled quotes
function split_csv_line(ln::AbstractString)
    out = String[]; buf = IOBuffer(); inq = false
    i = firstindex(ln); n = lastindex(ln)
    while i <= n
        ch = ln[i]
        if inq
            if ch == '"'
                if i < n && ln[nextind(ln, i)] == '"'
                    write(buf, '"'); i = nextind(ln, i)
                else
                    inq = false
                end
            else
                write(buf, ch)
            end
        elseif ch == '"'
            inq = true
        elseif ch == ','
            push!(out, String(take!(buf)))
        else
            write(buf, ch)
        end
        i = nextind(ln, i)
    end
    push!(out, String(take!(buf)))
    return out
end

function read_table(path)
    lines = filter(!isempty, readlines(path))
    isempty(lines) && error("empty csv: $path")
    hdr = split_csv_line(lines[1])
    col = Dict(h => i for (i, h) in enumerate(hdr))
    for c in vcat(KEY_COLS, CHECK_COLS)
        haskey(col, c) || error("$path has no column \"$c\"")
    end
    rows = Dict{NTuple{6, String}, Dict{String, String}}()
    for ln in lines[2:end]
        fl = split_csv_line(ln)
        length(fl) == length(hdr) || error("ragged row in $path: $(first(ln, 80))")
        d = Dict(h => fl[i] for (h, i) in col)
        rows[Tuple(String(d[c]) for c in KEY_COLS)] = d
    end
    return rows
end

_gap(d) = something(tryparse(Float64, get(d, "gap_pct", "")), NaN)

function gap_note(ge, gg)
    (isnan(ge) || isnan(gg)) && return ""
    d = abs(ge - gg)
    return d < 1e-9 ? "  gap identical" : "  gap expected $(ge) got $(gg) (diff $(d))"
end

function main()
    length(ARGS) >= 2 || error("usage: verify_crossing.jl <crossing_table.csv> <expected.csv>...")
    got = read_table(ARGS[1])
    n_ok = 0; n_bad = 0; n_missing = 0
    for exp_path in ARGS[2:end]
        expected = read_table(exp_path)
        for key in sort(collect(keys(expected)))
            e = expected[key]
            label = join(key, "/")
            if !haskey(got, key)
                n_missing += 1
                println("missing   $label")
                continue
            end
            g = got[key]
            diffs = [c for c in CHECK_COLS if g[c] != e[c]]
            if isempty(diffs)
                n_ok += 1
                println("ok        $label  a=$(e["tier1_first_touch"]) b=$(e["tier2_iter"]) " *
                        "$(e["verdict"])" * gap_note(_gap(e), _gap(g)))
            else
                n_bad += 1
                println("MISMATCH  $label")
                for c in diffs
                    println("            $c: expected $(e[c]) got $(g[c])")
                end
            end
        end
    end
    println()
    println("verify_crossing: $n_ok reproduce, $n_bad mismatch, $n_missing missing")
    exit(n_bad == 0 ? 0 : 1)
end

main()
