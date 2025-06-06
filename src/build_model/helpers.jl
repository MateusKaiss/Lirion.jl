"""
Selects a heuristic scaling factor for conic constraints (e.g., SOCP)
based on the voltage magnitude bounds in the system.

This function was empirically developed. It observes that large discrepancies
between `vmin` and `vmax` (delta) tend to correlate with numerical instability
when using conic relaxations. A larger delta generally benefits from higher
scaling values to prevent solver issues (e.g., infeasibility or NaNs in cones).

- If the system is large (>1000 buses) and has a wide voltage range (delta > 1.17), scale = 100
- If the voltage range is tight (delta < 1.13), scale = 20 (safe and precise)
- Moderate cases (delta < 1.16) use scale = 70
- Default fallback is scale = 100

These values were chosen based on testing and iterative tuning.
"""
function select_cone_scale(params)
    vmin_vals = collect(values(params.vmin))
    vmax_vals = collect(values(params.vmax))

    vmin = clamp(minimum(vmin_vals), 0.80, 1.05)
    vmax = clamp(maximum(vmax_vals), 0.90, 1.15)
    delta = vmax / vmin

    nb = length(params.vmin)
    if nb > 1000 && delta > 1.17
        return 100
    elseif delta < 1.13
        return 20
    elseif delta < 1.16
        return 70
    else
        return 100
    end
end