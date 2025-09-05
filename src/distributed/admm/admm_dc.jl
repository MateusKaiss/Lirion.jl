# =============== DC Utils  ===============
# ---------- Solver ----------

function resolve_optimizer_dc(config)
    opt = try
        getfield(config, :optimizer)
    catch
        nothing
    end
    if opt !== nothing
        return opt
    end

    eps_abs = try
        getfield(config, Symbol("ε_primal"))
    catch
        1e-6
    end
    eps_rel = try
        getfield(config, Symbol("ε_dual"))
    catch
        1e-6
    end

    return optimizer_with_attributes(OSQP.Optimizer,
        "eps_abs" => eps_abs,
        "eps_rel" => eps_rel,
        "max_iter" => 50_000,
        "verbose" => false,
        "warm_start" => true,
        "adaptive_rho" => true,
        "linsys_solver" => "qdldl",
    )
end

@inline function _dc_flow_val_model(
    m::Model,
    sub::Dict{String, Any},
    l_local::Int;
    t::Int = 1,
)
    l_global = sub["branch_local_indices"][l_local]
    f = sub["branch_full"]["f_bus"][l_global]
    tb = sub["branch_full"]["t_bus"][l_global]
    x = sub["branch_full"]["br_x"][l_global]
    @assert abs(x) > 1e-8 "br_x≈0 no ramo $l_global (DC)."
    return (value(m[:theta][f, t]) - value(m[:theta][tb, t])) / x
end

# ---------- Penalidades ----------
function add_penalty_bus_dc!(m, subdata, shared_buses::Vector{Int},
    θ_avg::Dict{Int, Float64}, duals_zone::Dict{Int, Float64}, ρ::Float64)
    local_shared = get(subdata, "shared_buses_local", shared_buses)
    pen = @expression(m, 0.0)
    @inbounds for b in local_shared
        θ = m[:theta][b, 1]
        θ̄ = get(θ_avg, b, 0.0)
        λ = get(duals_zone, b, 0.0)
        pen += λ * (θ - θ̄) + (ρ / 2) * (θ - θ̄)^2
    end
    return pen
end

function add_penalty_branch_dc!(m, subdata, zone_edges::Vector{Int},
    z_flow::Dict{Int, Float64}, duals_zone::Dict{Int, Float64}, ρ::Float64)
    pen = @expression(m, 0.0)
    @inbounds for e in zone_edges
        l_local = get(subdata["edge_to_local"], e, nothing)
        @assert l_local !== nothing "aresta global $e não pertence a esta partição"
        p_loc = _dc_flow_expr(m, subdata, l_local)   # assumido definido noutro arquivo
        λ = get(duals_zone, e, 0.0)
        z = z_flow[e]
        pen += λ * (p_loc - z) + (ρ / 2) * (p_loc - z)^2
    end
    return pen
end

# ---------- Consensos (2 + 3) ----------
function compute_theta_consensus!(θ_avg, partitions, shared_buses)
    sums = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()
    for (_, sub) in partitions
        m = sub["model"]
        @inbounds for b in sub["shared_buses_local"]
            θb = value(m[:theta][b, 1])
            sums[b] = get(sums, b, 0.0) + θb
            cnts[b] = get(cnts, b, 0) + 1
        end
    end
    @inbounds for b in shared_buses
        n = get(cnts, b, 0)
        if n > 0
            θ_avg[b] = sums[b] / n
        end
    end
end

function compute_z_flow_consensus!(
    z_flow::Dict{Int, Float64},
    partitions,
    cut_edges::Vector{Int},
)
    sums = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()
    for (_, sub) in partitions
        m = sub["model"]
        @inbounds for e in cut_edges
            l_local = get(sub["edge_to_local"], e, nothing)
            if l_local !== nothing
                p = _dc_flow_val_model(m, sub, l_local)
                sums[e] = get(sums, e, 0.0) + p
                cnts[e] = get(cnts, e, 0) + 1
            end
        end
    end
    @inbounds for e in cut_edges
        n = get(cnts, e, 0)
        if n > 0
            z_flow[e] = sums[e] / n
        end
    end
end

# ---------- Duais com over-relaxation opcional ----------
function update_duals_theta!(
    duals::Dict{Int, Dict{Int, Float64}}, partitions,
    θ_avg::Dict{Int, Float64}, θ_avg_prev::Dict{Int, Float64},
    _shared_buses::Vector{Int}, ρ::Float64, α::Union{Nothing, Float64},
)
    for (z, sub) in partitions
        m = sub["model"]
        @inbounds for b in sub["shared_buses_local"]
            θb = value(m[:theta][b, 1])
            if isnothing(α)
                duals[z][b] = get(duals[z], b, 0.0) + ρ * (θb - θ_avg[b])
            else
                θ_prev̄ = get(θ_avg_prev, b, get(θ_avg, b, 0.0))
                θ_rel = α * θb + (1 - α) * θ_prev̄
                duals[z][b] = get(duals[z], b, 0.0) + ρ * (θ_rel - θ_avg[b])
            end
        end
    end
end

function update_duals_flow!(
    duals::Dict{Int, Dict{Int, Float64}}, partitions,
    z_flow::Dict{Int, Float64}, z_flow_prev::Dict{Int, Float64},
    cut_edges::Vector{Int}, ρ::Float64, α::Union{Nothing, Float64},
)
    for (z, sub) in partitions
        m = sub["model"]
        @inbounds for e in cut_edges
            l_local = get(sub["edge_to_local"], e, nothing)
            if l_local !== nothing
                p_val = _dc_flow_val_model(m, sub, l_local)
                if isnothing(α)
                    duals[z][e] = get(duals[z], e, 0.0) + ρ * (p_val - z_flow[e])
                else
                    z_prev = get(z_flow_prev, e, z_flow[e])
                    p_rel = α * p_val + (1 - α) * z_prev
                    duals[z][e] = get(duals[z], e, 0.0) + ρ * (p_rel - z_flow[e])
                end
            end
        end
    end
end

# ---------- Resíduos (2 + 3) ----------
function primal_residual_bus_dc(
    partitions,
    θ_avg::Dict{Int, Float64},
    _shared_buses::Vector{Int},
)
    acc = 0.0
    for (_, sub) in partitions
        m = sub["model"]
        @inbounds for b in sub["shared_buses_local"]
            acc += (value(m[:theta][b, 1]) - get(θ_avg, b, 0.0))^2
        end
    end
    return sqrt(acc)
end

dual_residual_bus_dc(θ_avg::Dict{Int, Float64}, θ_avg_prev::Dict{Int, Float64}) =
    sqrt(sum((θ_avg[b] - get(θ_avg_prev, b, 0.0))^2 for b in keys(θ_avg)))

function eps_primal_bus_dc(partitions, θ_avg::Dict{Int, Float64}, shared_buses::Vector{Int},
    ϵ_abs::Float64, ϵ_rel::Float64)
    N = max(1, length(shared_buses))
    norm_x_sq = 0.0
    for (_, sub) in partitions
        m = sub["model"]
        @inbounds for b in sub["shared_buses_local"]
            norm_x_sq += (value(m[:theta][b, 1]))^2
        end
    end
    norm_z_sq = sum((get(θ_avg, b, 0.0))^2 for b in shared_buses)
    return sqrt(N) * ϵ_abs + ϵ_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function eps_dual_bus_dc(θ_avg::Dict{Int, Float64}, θ_avg_prev::Dict{Int, Float64},
    ϵ_abs::Float64, ϵ_rel::Float64)
    N = max(1, length(θ_avg))
    return sqrt(N) * ϵ_abs +
           ϵ_rel * sqrt(sum((θ_avg[b] - get(θ_avg_prev, b, 0.0))^2 for b in keys(θ_avg)))
end

function primal_residual_branch_dc(
    partitions,
    z_flow::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    acc = 0.0
    for (_, sub) in partitions
        m = sub["model"]
        @inbounds for e in cut_edges
            l_local = get(sub["edge_to_local"], e, nothing)
            if l_local !== nothing
                acc += (_dc_flow_val_model(m, sub, l_local) - z_flow[e])^2
            end
        end
    end
    return sqrt(acc)
end

dual_residual_branch_dc(z_flow::Dict{Int, Float64}, z_flow_prev::Dict{Int, Float64}) =
    sqrt(sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in keys(z_flow)))

function eps_primal_branch_dc(partitions, z_flow::Dict{Int, Float64},
    cut_edges::Vector{Int}, ϵ_abs::Float64, ϵ_rel::Float64)
    dim_pri = sum(length(_zone_cut_edges(sub, cut_edges)) for (_, sub) in partitions)
    return sqrt(max(1, dim_pri)) * ϵ_abs +
           ϵ_rel * max(sqrt(sum(z_flow[e]^2 for e in cut_edges)), 1e-12)
end

function eps_dual_branch_dc(z_flow::Dict{Int, Float64}, z_flow_prev::Dict{Int, Float64},
    ϵ_abs::Float64, ϵ_rel::Float64)
    dim_du = max(1, length(z_flow))
    return sqrt(dim_du) * ϵ_abs +
           ϵ_rel * sqrt(sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in keys(z_flow)))
end

# ---------- ρ adaptativo ----------
function adaptive_rho_dc!(ρ::Base.RefValue{Float64}, r::Float64, s::Float64;
    μ_incr = 2.0, μ_decr = 2.0, τ_incr = 1.05, τ_decr = 0.98)
    if r > μ_incr * s
        ρ[] *= τ_incr
    elseif s > μ_decr * r
        ρ[] *= τ_decr
    end
end

# ---------- Plots (wrappers) ----------
function save_dc_bus_plots(
    res_pri,
    res_du,
    eps_pri,
    eps_du,
    ρ_hist,
    θ_avg_history;
    stem::String = "admm",
)
    plot_admm_convergence(
        res_pri,
        res_du,
        eps_pri,
        eps_du,
        "plot_convergencia_$(stem)_dc.html",
    )
    plot_rho(ρ_hist, "plot_rho_$(stem)_dc.html")
    return plot_theta_avg_heatmap(θ_avg_history, "plot_theta_avg_heatmap_$(stem)_dc.html")
end

function save_dc_branch_plots(
    res_pri,
    res_du,
    eps_pri,
    eps_du,
    ρ_hist;
    stem::String = "admm_branch",
)
    plot_admm_convergence(
        res_pri,
        res_du,
        eps_pri,
        eps_du,
        "plot_convergencia_$(stem).html",
    )
    return plot_rho(ρ_hist, "plot_rho_$(stem).html")
end

# ---------- ADMM: BUS ----------
function run_admm(config::DistributedOPFModel, data)
    partitions, shared_buses = prepare_admm_partitions(data)
    @assert config.base_model == DC "run_admm (bus/DC) requer DC."

    ρ = Ref(config.rho_theta)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    duals =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    θ_avg = Dict{Int, Float64}()
    θ_avg_prev = Dict{Int, Float64}()

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rho_vec = Float64[]
    θ_avg_history = Dict{Int, Vector{Float64}}()
    θ_local_history = Dict{Int, Dict{Int, Vector{Float64}}}()

    optimizer_factory = resolve_optimizer_dc(config)

    for (zone_id, sub) in partitions
        m = build_model(sub, DC; optimizer = optimizer_factory, build_only=true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["theta_var"] = m[:theta]
        sub["has_solution"] = false
        _set_model_silent!(m)
        @inbounds for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
    end

    for k in 1:config.max_iter
        println("\nIteração $k (bus/DC)")

        for (zone_id, sub) in partitions
            m = sub["model"]

            pen = add_penalty_bus_dc!(m, sub, shared_buses, θ_avg, duals[zone_id], ρ[])
            @objective(m, Min, sub["obj_base"] + pen)
            optimize!(m)
            sub["has_solution"] = true

            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] =
                compute_real_generation_cost(m, sub["objective_params"])

            @inbounds for b in sub["shared_buses_local"]
                push!(
                    get!(
                        get!(θ_local_history, zone_id, Dict{Int, Vector{Float64}}()),
                        b,
                        Float64[],
                    ),
                    value(sub["theta_var"][b, 1]))
            end
        end

        θ_avg_prev = isempty(θ_avg) ? θ_avg : copy(θ_avg)
        compute_theta_consensus!(θ_avg, partitions, shared_buses)
        for b in shared_buses
            push!(get!(θ_avg_history, b, Float64[]), get(θ_avg, b, 0.0))
        end

        update_duals_theta!(
            duals,
            partitions,
            θ_avg,
            θ_avg_prev,
            shared_buses,
            ρ[],
            config.alpha_overrelax,
        )

        res_primal = primal_residual_bus_dc(partitions, θ_avg, shared_buses)
        res_dual = dual_residual_bus_dc(θ_avg, θ_avg_prev)
        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        ϵ_pri = eps_primal_bus_dc(partitions, θ_avg, shared_buses, ϵ_abs, ϵ_rel)
        ϵ_du = eps_dual_bus_dc(θ_avg, θ_avg_prev, ϵ_abs, ϵ_rel)
        push!(eps_primal_vec, ϵ_pri)
        push!(eps_dual_vec, ϵ_du)

        println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        println("‖s‖ = $res_dual   | ε_dual = $ϵ_du")

        push!(rho_vec, ρ[])

        if res_primal < ϵ_pri && res_dual < ϵ_du
            println("Stopping criterion met.")
            break
        end

        adaptive_rho_dc!(ρ, res_primal, res_dual)
    end

    save_dc_bus_plots(
        res_primal_vec,
        res_dual_vec,
        eps_primal_vec,
        eps_dual_vec,
        rho_vec,
        θ_avg_history;
        stem = "admm",
    )

    total_real_cost = sum(sub["real_generation_cost"] for sub in values(partitions))
    total_expr_cost = sum(sub["generation_cost_model_expr"] for sub in values(partitions))

    println("\nResumo (bus/DC):")
    println("→ Total (value(:generation_cost)) = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ inicial: $(config.rho_theta) | ρθ final: $(ρ[])")

    return total_real_cost, partitions
end

# ---------- ADMM: BRANCH ----------
function run_admm_branch(
    config::DistributedOPFModel,
    data::Dict{String, Any};
    warm_start_theta::Bool = false,
    use_legacy_rho::Bool = false,
    μ_incr::Float64 = 2.0,
    μ_decr::Float64 = 2.0,
    τ_incr::Float64 = 1.05,
    τ_decr::Float64 = 0.98,
)
    @assert config.base_model == DC "ADMM-branch/DC requer base_model == DC."

    partitions, _ = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    ρ = Ref(config.rho_flow)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    z_flow = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    z_flow_prev = Dict{Int, Float64}()
    duals =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rho_vec = Float64[]

    optimizer_factory = resolve_optimizer_dc(config)
    for (zone_id, sub) in partitions
        m = build_model(sub, DC; optimizer = optimizer_factory, build_only=true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["theta_var"] = m[:theta]
        sub["zone_edges"] = _zone_cut_edges(sub, cut_edges)
        sub["has_solution"] = false
        _set_model_silent!(m)
        @inbounds for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
    end

    for k in 1:config.max_iter
        println("\n[ADMM-branch/DC] Iteração $k")

        for (zone_id, sub) in partitions
            m = sub["model"]

            # if warm_start_theta && sub["has_solution"] === true
            #     @inbounds for I in keys(sub["theta_var"])
            #         set_start_value(sub["theta_var"][I], value(sub["theta_var"][I]))
            #     end
            # end

            pen = add_penalty_branch_dc!(
                m,
                sub,
                sub["zone_edges"],
                z_flow,
                duals[zone_id],
                ρ[],
            )
            @objective(m, Min, sub["obj_base"] + pen)
            optimize!(m)
            sub["has_solution"] = true

            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] =
                compute_real_generation_cost(m, sub["objective_params"])
        end

        z_flow_prev = copy(z_flow)
        compute_z_flow_consensus!(z_flow, partitions, cut_edges)
        update_duals_flow!(
            duals,
            partitions,
            z_flow,
            z_flow_prev,
            cut_edges,
            ρ[],
            config.alpha_overrelax,
        )

        res_primal = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        res_dual = dual_residual_branch_dc(z_flow, z_flow_prev)
        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        ϵ_pri = eps_primal_branch_dc(partitions, z_flow, cut_edges, ϵ_abs, ϵ_rel)
        ϵ_du = eps_dual_branch_dc(z_flow, z_flow_prev, ϵ_abs, ϵ_rel)
        push!(eps_primal_vec, ϵ_pri)
        push!(eps_dual_vec, ϵ_du)

        println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        println("‖s‖ = $res_dual   | ε_dual = $ϵ_du")

        push!(rho_vec, ρ[])

        if res_primal < ϵ_pri && res_dual < ϵ_du
            println("Stopping criterion met.")
            break
        end

        if use_legacy_rho
            if res_primal > μ_incr * res_dual
                ρ[] *= τ_incr
            elseif res_dual > μ_decr * res_primal
                ρ[] *= τ_decr
            end
        else
            adaptive_rho_dc!(
                ρ,
                res_primal,
                res_dual;
                μ_incr = μ_incr,
                μ_decr = μ_decr,
                τ_incr = τ_incr,
                τ_decr = τ_decr,
            )
        end
    end

    save_dc_branch_plots(
        res_primal_vec,
        res_dual_vec,
        eps_primal_vec,
        eps_dual_vec,
        rho_vec;
        stem = "admm_branch",
    )

    total_real_cost = sum(sub["real_generation_cost"] for sub in values(partitions))
    total_expr_cost = sum(sub["generation_cost_model_expr"] for sub in values(partitions))

    println("\n[ADMM-branch/DC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρe_final = $(ρ[])")

    return total_real_cost, partitions
end

# ---------- ADMM: HÍBRIDO ----------
function run_admm_hybrid(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == DC "Hybrid ADMM/DC requer DC."

    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    ρθ = Ref(config.rho_theta)
    ρe = Ref(config.rho_flow)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    θ_avg = Dict{Int, Float64}()
    θ_avg_prev = Dict{Int, Float64}()
    z_flow = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    z_flow_prev = Dict{Int, Float64}()

    duals_theta =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    duals_flow =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rhoθ_vec = Float64[]
    rhoe_vec = Float64[]
    θ_avg_history = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_dc(config)

    for (zone_id, sub) in partitions
        m = build_model(sub, DC; optimizer = optimizer_factory, build_only=true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["theta_var"] = m[:theta]
        sub["zone_edges"] = _zone_cut_edges(sub, cut_edges)
        sub["has_solution"] = false
        _set_model_silent!(m)
        @inbounds for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
    end

    for k in 1:config.max_iter
        println("\n[ADMM-hybrid/DC] Iteração $k")

        for (zone_id, sub) in partitions
            m = sub["model"]

            pen_bus =
                add_penalty_bus_dc!(m, sub, shared_buses, θ_avg, duals_theta[zone_id], ρθ[])
            pen_flow = add_penalty_branch_dc!(
                m,
                sub,
                sub["zone_edges"],
                z_flow,
                duals_flow[zone_id],
                ρe[],
            )
            @objective(m, Min, sub["obj_base"] + pen_bus + pen_flow)
            optimize!(m)
            sub["has_solution"] = true

            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] =
                compute_real_generation_cost(m, sub["objective_params"])
        end

        θ_avg_prev = isempty(θ_avg) ? θ_avg : copy(θ_avg)
        z_flow_prev = copy(z_flow)

        compute_theta_consensus!(θ_avg, partitions, shared_buses)
        for b in shared_buses
            push!(get!(θ_avg_history, b, Float64[]), get(θ_avg, b, 0.0))
        end
        compute_z_flow_consensus!(z_flow, partitions, cut_edges)

        update_duals_theta!(
            duals_theta,
            partitions,
            θ_avg,
            θ_avg_prev,
            shared_buses,
            ρθ[],
            config.alpha_overrelax,
        )
        update_duals_flow!(
            duals_flow,
            partitions,
            z_flow,
            z_flow_prev,
            cut_edges,
            ρe[],
            config.alpha_overrelax,
        )

        res_pri_bus = primal_residual_bus_dc(partitions, θ_avg, shared_buses)
        res_pri_br = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        res_primal = sqrt(res_pri_bus^2 + res_pri_br^2)

        res_du_bus = dual_residual_bus_dc(θ_avg, θ_avg_prev)
        res_du_br = dual_residual_branch_dc(z_flow, z_flow_prev)
        res_dual = sqrt(res_du_bus^2 + res_du_br^2)

        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        ϵ_pri = begin
            ϵb = eps_primal_bus_dc(partitions, θ_avg, shared_buses, ϵ_abs, ϵ_rel)
            ϵe = eps_primal_branch_dc(partitions, z_flow, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(ϵb^2 + ϵe^2)
        end
        ϵ_du = begin
            ϵb = eps_dual_bus_dc(θ_avg, θ_avg_prev, ϵ_abs, ϵ_rel)
            ϵe = eps_dual_branch_dc(z_flow, z_flow_prev, ϵ_abs, ϵ_rel)
            sqrt(ϵb^2 + ϵe^2)
        end
        push!(eps_primal_vec, ϵ_pri)
        push!(eps_dual_vec, ϵ_du)

        println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        println("‖s‖ = $res_dual   | ε_dual = $ϵ_du")

        push!(rhoθ_vec, ρθ[])
        push!(rhoe_vec, ρe[])

        if res_primal < ϵ_pri && res_dual < ϵ_du
            println("Stopping criterion met.")
            break
        end

        adaptive_rho_dc!(ρθ, res_primal, res_dual)
        adaptive_rho_dc!(ρe, res_primal, res_dual)
    end

    plot_admm_convergence(
        res_primal_vec,
        res_dual_vec,
        eps_primal_vec,
        eps_dual_vec,
        "plot_convergencia_admm_hybrid.html",
    )

    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoθ_vec),
                    y = rhoθ_vec,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρθ",
                );
                p.layout.title = "Evolução de ρθ";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρθ";
                p.layout.yaxis.type = "log";
                p
            ))(),
        "plot_rho_theta_hybrid.html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoe_vec),
                    y = rhoe_vec,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρe",
                );
                p.layout.title = "Evolução de ρe";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρe";
                p.layout.yaxis.type = "log";
                p
            ))(),
        "plot_rho_flow_hybrid.html",
    )

    plot_theta_avg_heatmap(θ_avg_history, "plot_theta_avg_heatmap_hybrid.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))

    println("\n[ADMM-hybrid/DC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ_final = $(ρθ[]) | ρe_final = $(ρe[])")

    return total_real_cost, partitions
end
