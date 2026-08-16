# ============================================================
# ALADIN Algorithm 3 — globalization branch (c): dual-only step
# ============================================================

function _vrho_resolve_sum!(ws::AbstractALADINWorkspace)
    parts = ws.partitions.parts
    vals = zeros(Float64, length(parts))
    Threads.@threads for i in eachindex(parts)
        part = parts[i]
        update_objective!(part, ws)
        JuMP.optimize!(part.local_model)
        _require_optimal!(part.local_model; iter = ws.iter, zone_id = part.id)
        vals[i] = JuMP.objective_value(part.local_model)
    end
    return sum(vals)
end

function _golden_max(f, a::Float64, b::Float64; iters::Int = 6)
    invphi = 0.6180339887498949            # 1/φ
    c = b - invphi * (b - a)
    d = a + invphi * (b - a)
    fc = f(c)
    fd = f(d)
    @inbounds for _ in 1:iters
        if fc < fd
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd = f(d)
        else
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc = f(c)
        end
    end
    return fc > fd ? c : d
end

function _branch_c_dual_linesearch!(ws::AbstractALADINWorkspace)
    states = collect(values(ws.vars))
    saved = [copy(s.lambda) for s in states]
    dir = [s.nu .- s.lambda for s in states]

    function setlambda!(alpha)
        @inbounds for (s, s0, d) in zip(states, saved, dir)
            @. s.lambda = s0 + alpha * d
        end
        return nothing
    end

    vrho(alpha) = (setlambda!(alpha); _vrho_resolve_sum!(ws))

    alpha3 = _golden_max(vrho, 0.0, 1.0; iters = 6)
    @inbounds for (s, s0) in zip(states, saved)
        copyto!(s.lambda, s0)
    end

    ws.alpha1 = 0.0
    ws.alpha2 = 0.0
    ws.alpha3 = alpha3
    return nothing
end
