"""
Returns the validated IPOPT linear solver name.

Validates and logs the linear solver choice. Returns the solver name as a string.

For installation, after downloading the HSL library, run the following commands on your shell:

# on your Downloads folder (or any other folder), extract the tar.gz file:
tar -xzf HSL_jll.jl.v2025.7.21.tar.gz

# then, run the following command in your Julia's project folder:
julia --project=. -e 'using Pkg; Pkg.develop(path="/home/username/Downloads/HSL_jll.jl.v2025.7.21")'

Replace the path as appropriate for your system.

Arguments:
- linear_solver: String - one of "ma27", "ma57", "ma77", "ma86", "ma97", "mumps" (default: "ma57")

Returns:
- String with the validated linear solver name
"""
function get_ipopt_linear_solver_option(linear_solver::String = "ma57")
    hsl_solvers = ["ma27", "ma57", "ma77", "ma86", "ma97"]
    
    if linear_solver in hsl_solvers
        @info "[IPOPT] Configuring HSL linear solver: $linear_solver"
        return linear_solver
    elseif linear_solver == "mumps"
        @info "[IPOPT] Configuring MUMPS linear solver"
        return "mumps"
    else
        @warn "[IPOPT] Unknown linear solver '$linear_solver'. Using default: ma57"
        @info "[IPOPT] Configuring HSL linear solver: ma57"
        return "ma57"
    end
end

"""
Configures IPOPT optimizer with HSL linear solver.

This is the shared function used by all IPOPT configurations (centralized and ADMM).
It ensures consistent linear solver configuration across the codebase.

Arguments:
- optimizer: The optimizer (should be Ipopt.Optimizer or nothing)
- ipopt_linear_solver: String - one of "ma27", "ma57", "ma77", "ma86", "ma97", "mumps" (default: "ma57")

Returns:
- Configured optimizer with linear solver option, or unchanged optimizer if not IPOPT
"""
function configure_ipopt_with_hsl(optimizer, ipopt_linear_solver::String = "ma57")
    if optimizer !== nothing && optimizer == Ipopt.Optimizer
        linear_solver_name = get_ipopt_linear_solver_option(ipopt_linear_solver)
        return optimizer_with_attributes(Ipopt.Optimizer, "linear_solver" => linear_solver_name)
    end
    return optimizer
end

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
