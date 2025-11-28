# =========================
# BUILD ENTRYPOINT
# =========================
function build_model(data, model_type; optimizer, build_only::Bool = false)
    if model_type == SOCP
        return try_with_fallback_scale(data, model_type; optimizer, build_only = build_only)
    else
        return build_and_optimize_model(
            data,
            model_type;
            optimizer,
            build_only = build_only,
        )
    end
end

# =========================
# CENTRALIZED (DC/AC/SOCP)
# =========================
function build_and_optimize_model(
    matpower_data,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
    build_only::Bool = false,
)
    optimizer === nothing && (optimizer = Ipopt.Optimizer)
    times = get_model_times(matpower_data)
    if length(unique(times)) > 1
        return build_and_optimize_multiperiod_model(
            matpower_data, model_type; optimizer, times = times, build_only = build_only,
        )
    else
        return build_and_optimize_singleperiod_model(
            matpower_data, model_type; optimizer, build_only = build_only,
        )
    end
end

# =========================================================
# SOCP: tenta escalas diferentes (a menos que build_only)
# =========================================================
"""
Attempts to solve the model using a sequence of scaling factors
to improve numerical stability for conic models (e.g., SOCP).

Tries up to 5 scale values: selected, 150, 1, 300, 500.
When `build_only=true`, it only builds once (no scaling loop, no optimize).
"""
function try_with_fallback_scale(
    matpower_data::Dict,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
    build_only::Bool = false,
)
    if build_only
        matpower_data_once = deepcopy(matpower_data)
        matpower_data_once["selected_scale"] = get(matpower_data, "selected_scale", nothing)
        return build_and_optimize_model(
            matpower_data_once,
            model_type;
            optimizer,
            build_only = true,
        )
    end

    # 1) escala selecionada a partir dos buses
    scale = select_cone_scale(matpower_data["bus_params"])
    matpower_data_attempt = deepcopy(matpower_data)
    matpower_data_attempt["selected_scale"] = scale
    model = build_and_optimize_model(matpower_data_attempt, model_type; optimizer)
    status = JuMP.termination_status(model)

    # 2) fallback 150
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Model with scale=$scale failed (status: $status). Trying scale=150..."
        matpower_data_fallback = deepcopy(matpower_data)
        matpower_data_fallback["selected_scale"] = 150
        model = build_and_optimize_model(matpower_data_fallback, model_type; optimizer)
        status = JuMP.termination_status(model)
    end

    # 3) minimal 1
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Model still failed with scale=150 (status: $status). Trying scale=1..."
        matpower_data_final_attempt = deepcopy(matpower_data)
        matpower_data_final_attempt["selected_scale"] = 1
        model = build_and_optimize_model(matpower_data_final_attempt, model_type; optimizer)
        status = JuMP.termination_status(model)
    end

    # 4) 300
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Model still failed with scale=1 (status: $status). Trying scale=300..."
        matpower_data_last_attempt = deepcopy(matpower_data)
        matpower_data_last_attempt["selected_scale"] = 300
        model = build_and_optimize_model(matpower_data_last_attempt, model_type; optimizer)
        status = JuMP.termination_status(model)
    end

    # 5) 500
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Model still failed with scale=300 (status: $status). Trying scale=500 (final attempt)..."
        matpower_data_extreme_attempt = deepcopy(matpower_data)
        matpower_data_extreme_attempt["selected_scale"] = 500
        model =
            build_and_optimize_model(matpower_data_extreme_attempt, model_type; optimizer)
    end

    return model
end

# ============================================
# SINGLE-PERIOD
# ============================================
function build_and_optimize_singleperiod_model(
    matpower_data,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
    build_only::Bool = false,
)
    model = JuMP.Model(optimizer)
    model.ext[:time_periods] = 1
    model.ext[:times] = [1]
    model.ext[:selected_scale] = get(matpower_data, "selected_scale", nothing)

    process_subproblem(model, matpower_data, model_type, SubproblemBuild)

    if !build_only
        JuMP.optimize!(model)
        print(model)
        termination_status = JuMP.termination_status(model)
        if termination_status == MOI.OPTIMAL || termination_status == MOI.LOCALLY_SOLVED
            @show JuMP.objective_value(model)
        else
            @warn "Modelo não resolveu com sucesso. Status: $termination_status"
            print(model)
            # compute_conflict!(model)  # habilite se quiser depurar
            # error("Falha na otimização do modelo single-period.")
        end
    end

    return model
end

# ============================================
# MULTI-PERIOD
# ============================================
function build_and_optimize_multiperiod_model(
    matpower_data,
    model_type::Type{<:OPFModel};
    optimizer = Ipopt.Optimizer,
    times::AbstractVector,
    build_only::Bool = false,
)
    model = JuMP.Model(optimizer)
    model.ext[:time_periods] = length(times)
    model.ext[:times] = times
    model.ext[:selected_scale] = get(matpower_data, "selected_scale", nothing)

    process_subproblem(model, matpower_data, model_type, SubproblemBuild)

    if !build_only
        JuMP.optimize!(model)
        @show JuMP.objective_value(model)
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
