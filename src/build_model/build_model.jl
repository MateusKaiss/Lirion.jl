# =========================
# BUILD ENTRYPOINT
# =========================
# Import the shared configuration function
# (configure_ipopt_with_hsl is defined in helpers.jl)

function _select_default_optimizer(model_type::Type{<:OPFModel}, opt, ipopt_linear_solver::String = "ma57")
    if opt !== nothing
        # Configure IPOPT if it's the selected optimizer
        return configure_ipopt_with_hsl(opt, ipopt_linear_solver)
    end
    if model_type == DC
        return OSQP.Optimizer
    elseif model_type == SOCP
        return SCS.Optimizer
    else
        # AC model defaults to IPOPT - configure with HSL
        return configure_ipopt_with_hsl(Ipopt.Optimizer, ipopt_linear_solver)
    end
end

function build_model(data, model_type; optimizer = nothing, build_only::Bool = false, ipopt_linear_solver::String = "ma57")
    return build_and_optimize_model(
        data,
        model_type;
        optimizer,
        build_only = build_only,
        ipopt_linear_solver = ipopt_linear_solver,
    )
end

# =========================
# CENTRALIZED (DC/AC/SOCP)
# =========================
function build_and_optimize_model(
    matpower_data::Dict{String, Any},
    model_type::Type{<:OPFModel};
    optimizer = nothing,
    build_only::Bool = false,
    ipopt_linear_solver::String = "ma57",
)
    # define optimizer default por tipo, se não veio de Args
    optimizer = _select_default_optimizer(model_type, optimizer, ipopt_linear_solver)

    # timeseries? single x multi
    times = get_model_times(matpower_data)
    if length(unique(times)) > 1
        return build_and_optimize_multiperiod_model(
            matpower_data, model_type; optimizer, times = times, build_only = build_only, ipopt_linear_solver = ipopt_linear_solver,
        )
    else
        return build_and_optimize_singleperiod_model(
            matpower_data, model_type; optimizer, build_only = build_only, ipopt_linear_solver = ipopt_linear_solver,
        )
    end
end

# ============================================
# SINGLE-PERIOD
# ============================================
function build_and_optimize_singleperiod_model(
    matpower_data::Dict{String, Any},
    model_type::Type{<:OPFModel};
    optimizer = nothing,
    build_only::Bool = false,
    ipopt_linear_solver::String = "ma57",
)
    optimizer === nothing && (optimizer = _select_default_optimizer(model_type, optimizer, ipopt_linear_solver))
    model = JuMP.Model(optimizer)
    model.ext[:time_periods] = 1
    model.ext[:times] = [1]

    process_subproblem(model, matpower_data, model_type, SubproblemBuild)

    if !build_only
        JuMP.optimize!(model)
        status = JuMP.termination_status(model)
        if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
            @info "Objective value = $(JuMP.objective_value(model)) | status=$status"
        else
            @warn "Modelo não resolveu com sucesso. Status: $status"
            # compute_conflict!(model)  # habilite se quiser depurar
        end
    end

    return model
end

# ============================================
# MULTI-PERIOD
# ============================================
function build_and_optimize_multiperiod_model(
    matpower_data::Dict{String, Any},
    model_type::Type{<:OPFModel};
    optimizer = nothing,
    times::AbstractVector,
    build_only::Bool = false,
    ipopt_linear_solver::String = "ma57",
)
    optimizer === nothing && (optimizer = _select_default_optimizer(model_type, optimizer, ipopt_linear_solver))
    model = JuMP.Model(optimizer)
    model.ext[:time_periods] = length(times)
    model.ext[:times] = times

    process_subproblem(model, matpower_data, model_type, SubproblemBuild)

    if !build_only
        JuMP.optimize!(model)
        status = JuMP.termination_status(model)
        if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
            @info "Objective value (multi) = $(JuMP.objective_value(model)) | status=$status"
        else
            @warn "Modelo multi-período não resolveu com sucesso. Status: $status"
        end
    end

    return model
end

# ============================================
# SUBPROBLEM PIPELINE
# ============================================
function process_subproblem(args...)
    buses!(args...)
    branches!(args...)
    generators!(args...)
    batteries!(args...)
    prosumer(args...)
    active_power_balance!(args...)
    reactive_power_balance!(args...)
    objective!(args...)
    return nothing
end
