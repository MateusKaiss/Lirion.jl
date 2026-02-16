# ============================================================
# GLOBAL FAILURE REGISTRY (genérico)
# ============================================================

const GLOBAL_FAILED_CASES = Dict{String, Dict{String, Vector{String}}}()

function track_failure(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] = Dict{String, Vector{String}}()
    end
    if !haskey(GLOBAL_FAILED_CASES[label], origin)
        GLOBAL_FAILED_CASES[label][origin] = String[]
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

function print_failures()
    println("\n===== GLOBAL SUMMARY OF FAILURES BY MODEL =====")
    for (label, failures) in GLOBAL_FAILED_CASES
        println("\n[$label]")
        for (origin, cases) in failures
            println("  $origin error in $(length(cases)) case(s):")
            for c in cases
                println("    - $c")
            end
        end
    end
    return println("==============================================\n")
end

# ============================================================
# FILE RESOLUTION
# ============================================================

function resolve_case_files(test_dir::String, case_list::Vector{String})
    return [
        joinpath(test_dir, fname)
        for fname in case_list
        if isfile(joinpath(test_dir, fname))
    ]
end

# ============================================================
# SOLVER HELPERS
# ============================================================

function get_available_optimizer(preferred, fallback)
    try
        preferred()
        return preferred
    catch
        @warn "Preferred optimizer not available. Falling back to $(fallback)"
        return fallback
    end
end

function _solve_maybe_optimizer(file::String; optimizer = nothing, kwargs...)
    if optimizer === nothing
        return Lirion.solve(file; kwargs...)
    else
        return Lirion.solve(file; optimizer = optimizer, kwargs...)
    end
end

# ============================================================
# SOCP solver picker
# ============================================================

const _SOCP_CHAIN_MAP = Dict(
    :auto => (Gurobi.Optimizer, SCS.Optimizer),
    :gurobi => (Gurobi.Optimizer, SCS.Optimizer),
    :scs => (SCS.Optimizer, Gurobi.Optimizer),
)

function pick_socp_optimizer(;
    prefer::Symbol = Symbol(lowercase(get(ENV, "LIRION_SOCP_PREFER", "auto"))),
)
    chain = get(_SOCP_CHAIN_MAP, prefer, _SOCP_CHAIN_MAP[:auto])
    for factory in chain
        try
            _ = factory()
            return factory
        catch
        end
    end
    return SCS.Optimizer
end

# ============================================================
# GENERIC DISTRIBUTED TEST RUNNER
# ============================================================

function run_case_suite(;
    label::String,
    test_dir::String,
    case_list::Vector{String}, model, central_optimizer,
    distributed_optimizer, objective_central::Function,
    objective_dist::Function, make_algorithm::Function,
    defaults::Function = _ -> NamedTuple(),
    overrides::Function = _ -> NamedTuple(), rtol::Real = 0.01,
    atol::Real = 1e-6,
)
    files = resolve_case_files(test_dir, case_list)
    n = length(files)

    if n == 0
        @info "[$label] No valid files found."
        return true
    end

    for (i, file) in enumerate(files)
        case_name = basename(file)

        base = defaults(case_name)
        over = overrides(case_name)
        params = merge(base, over)

        central_opt_case =
            hasproperty(params, :central_optimizer) ?
            getproperty(params, :central_optimizer) :
            central_optimizer

        dist_opt_case =
            hasproperty(params, :distributed_optimizer) ?
            getproperty(params, :distributed_optimizer) :
            distributed_optimizer

        @info "[$label] Running case $i/$n: $case_name | params=$(params) | " *
              "central_opt=$(central_opt_case === nothing ? "default" : string(central_opt_case)) | " *
              "dist_opt=$(dist_opt_case === nothing ? "default" : string(dist_opt_case))"

        @testset "$(label) - $(case_name)" begin
            centralized_obj = try
                out = _solve_maybe_optimizer(
                    file;
                    model = model,
                    algorithm = Lirion.centralized(),
                    optimizer = central_opt_case,
                    outputs = false,
                    plots = false,
                )
                objective_central(out)
            catch e
                @warn "Centralized failed on $case_name: $e"
                track_failure(label, "Centralized", case_name)
                return
            end

            distributed_obj = try
                alg = make_algorithm(params)
                out = _solve_maybe_optimizer(
                    file;
                    model = model,
                    algorithm = alg,
                    optimizer = dist_opt_case,
                    outputs = false,
                    plots = false,
                )
                objective_dist(out)
            catch e
                @warn "Distributed failed on $case_name: $e"
                track_failure(label, "Distributed", case_name)
                return
            end

            @info "OK $case_name | centralized=$(round(centralized_obj,digits=6)) distributed=$(round(distributed_obj,digits=6))"
            @test isapprox(distributed_obj, centralized_obj; rtol = rtol, atol = atol)
        end
    end

    return true
end
