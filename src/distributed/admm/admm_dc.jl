# =============== DC Utils (FAST DROP-IN)  ===============
# Mudanças:
# (1.1) Snapshot vetorial de θ + índices precomputados
# (1.2) Penalidade construída 1x (parte QUADRÁTICA via JuMP) + UPDATE só LINEAR por iteração
# (1.3) Consenso e duais com vetores densos
# (1.4) Reinstala Q SOMENTE quando ρ muda (bus, branch, hybrid)

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

# =======================================================
# HOTPATHS & SNAPSHOTS (1.1)
# =======================================================
function prepare_hotpaths!(sub::Dict{String, Any})
    bus_ids = sub["bus"]["bus_i"]
    bus_local_index = Dict{Int, Int}()
    @inbounds for (i, b) in enumerate(bus_ids)
        bus_local_index[b] = i
    end
    sub["bus_local_index"] = bus_local_index

    local_shared = get(sub, "shared_buses_local", Int[])
    shared_pos_local = Vector{Int}(undef, length(local_shared))
    @inbounds for (k, b) in enumerate(local_shared)
        shared_pos_local[k] = bus_local_index[b]
    end
    sub["shared_pos_local"] = shared_pos_local
    sub["theta_buf"] = zeros(Float64, length(bus_ids))
    return nothing
end

@inline function snapshot_theta!(sub::Dict{String, Any})
    θ = sub["model"][:theta]
    buf = sub["theta_buf"]
    bus_ids = sub["bus"]["bus_i"]
    @inbounds for i in eachindex(bus_ids)
        buf[i] = value(θ[bus_ids[i], 1])
    end
    return buf
end

# =======================================================
# FLUXO DC (para branch/híbrido)
# =======================================================
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

# =======================================================
# OBJETIVO PERSISTENTE — JuMP puro (1.2 + 1.4)
#  - Cache da parte QUADRÁTICA; reconstruída apenas quando ρ muda
#  - Parte LINEAR atualizada a cada iteração
# =======================================================

# BUS: constrói e cacheia PARTE QUADRÁTICA fixa: (ρ/2)*θ_b^2
function _build_quad_bus_expr!(sub::Dict{String, Any}, ρ::Float64)
    m = sub["model"]
    θ = m[:theta]
    q = zero(JuMP.QuadExpr)
    for b in sub["shared_buses_local"]
        q += 0.5 * ρ * θ[b, 1]^2
    end
    sub["quad_bus_expr"] = q
    sub["ρ_bus_cached"] = ρ
    return q
end

# BRANCH: constrói e cacheia PARTE QUADRÁTICA fixa: (ρ/2)*((θf-θt)/x)^2
function _build_quad_branch_expr!(
    sub::Dict{String, Any},
    ρ::Float64,
    zone_edges::Vector{Int},
)
    m = sub["model"]
    θ = m[:theta]
    q = zero(JuMP.QuadExpr)
    for e in zone_edges
        l_local = get(sub["edge_to_local"], e, nothing)
        l_local === nothing && continue
        l_global = sub["branch_local_indices"][l_local]
        f = sub["branch_full"]["f_bus"][l_global]
        t = sub["branch_full"]["t_bus"][l_global]
        x = sub["branch_full"]["br_x"][l_global]
        @assert abs(x) > 1e-10 "br_x≈0 no ramo $l_global (DC)."
        p = (θ[f, 1] - θ[t, 1]) / x
        q += 0.5 * ρ * p^2
    end
    sub["quad_branch_expr"] = q
    sub["ρ_branch_cached"] = ρ
    return q
end

# BUS: parte LINEAR por iteração: ∑ (λ - ρ*θ̄) * θ_b
function _build_linear_bus_expr(sub::Dict{String, Any},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    θ̄_vec::Vector{Float64},
    λ_vec_zone::Vector{Float64},
    ρ::Float64)
    m = sub["model"]
    θ = m[:theta]
    lin = zero(JuMP.AffExpr)
    for b in sub["shared_buses_local"]
        j = b2pos[b]
        lin += (λ_vec_zone[j] - ρ * θ̄_vec[j]) * θ[b, 1]
    end
    return lin
end

# BRANCH: parte LINEAR por iteração: ∑ (λe - ρ*z_e) * ((θf-θt)/x)
function _build_linear_branch_expr(sub::Dict{String, Any},
    zone_edges::Vector{Int},
    z_flow::Dict{Int, Float64},
    λe_zone::Dict{Int, Float64},
    ρ::Float64)
    m = sub["model"]
    θ = m[:theta]
    lin = zero(JuMP.AffExpr)
    for e in zone_edges
        l_local = get(sub["edge_to_local"], e, nothing)
        l_local === nothing && continue
        l_global = sub["branch_local_indices"][l_local]
        f = sub["branch_full"]["f_bus"][l_global]
        t = sub["branch_full"]["t_bus"][l_global]
        x = sub["branch_full"]["br_x"][l_global]
        @assert abs(x) > 1e-10 "br_x≈0 no ramo $l_global (DC)."
        coef = get(λe_zone, e, 0.0) - ρ * z_flow[e]
        lin += coef * (θ[f, 1] - θ[t, 1]) / x
    end
    return lin
end

# =======================================================
# CONSENSO / DUAIS DENSOS (1.3)
# =======================================================
function prepare_dense_maps(shared_buses::Vector{Int})
    shared_buses_sorted = sort(shared_buses)
    b2pos = Dict{Int, Int}(b => i for (i, b) in enumerate(shared_buses_sorted))
    return shared_buses_sorted, b2pos
end

function compute_theta_consensus_vec!(θ̄_vec::Vector{Float64},
    partitions::Dict{Int, Dict{String, Any}},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int})
    fill!(θ̄_vec, 0.0)
    counts = zeros(Int, length(shared_buses_sorted))
    for (_, sub) in partitions
        θbuf = snapshot_theta!(sub)
        pos = sub["bus_local_index"]
        @inbounds for b in sub["shared_buses_local"]
            j = b2pos[b]
            θ̄_vec[j] += θbuf[pos[b]]
            counts[j] += 1
        end
    end
    @inbounds for j in eachindex(θ̄_vec)
        counts[j] != 0 && (θ̄_vec[j] /= counts[j])
    end
end

function compute_z_flow_consensus!(
    z_flow::Dict{Int, Float64},
    partitions::Dict{Int, Dict{String, Any}},
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
        n > 0 && (z_flow[e] = sums[e] / n)
    end
end

@inline function update_duals_theta_vec!(dual_vec_zone::Vector{Float64},
    sub::Dict{String, Any},
    θ̄_vec::Vector{Float64},
    θ̄_prev_vec::Vector{Float64},
    ρ::Float64,
    α::Union{Nothing, Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int})
    θbuf = snapshot_theta!(sub)
    pos = sub["bus_local_index"]
    if isnothing(α)
        @inbounds for b in sub["shared_buses_local"]
            j = b2pos[b]
            dual_vec_zone[j] += ρ * (θbuf[pos[b]] - θ̄_vec[j])
        end
    else
        αv = something(α, 1.0)
        @inbounds for b in sub["shared_buses_local"]
            j = b2pos[b]
            θrel = αv * θbuf[pos[b]] + (1 - αv) * θ̄_prev_vec[j]
            dual_vec_zone[j] += ρ * (θrel - θ̄_vec[j])
        end
    end
end

@inline function gather_zone_lambda!(λ_vec_zone::Vector{Float64},
    duals_dense_zone::Vector{Float64})
    λ_vec_zone .= duals_dense_zone
    return λ_vec_zone
end

# =======================================================
# RESÍDUOS / TOLERÂNCIAS
# =======================================================
function primal_residual_bus_dc(
    partitions::Dict{Int, Dict{String, Any}},
    θ̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    acc = 0.0
    for (_, sub) in partitions
        θbuf = snapshot_theta!(sub)
        pos = sub["bus_local_index"]
        @inbounds for b in sub["shared_buses_local"]
            j = b2pos[b]
            acc += (θbuf[pos[b]] - θ̄_vec[j])^2
        end
    end
    return sqrt(acc)
end

dual_residual_bus_dc(θ̄_vec::Vector{Float64}, θ̄_prev_vec::Vector{Float64}) =
    sqrt(sum((θ̄_vec[j] - θ̄_prev_vec[j])^2 for j in eachindex(θ̄_vec)))

function eps_primal_bus_dc(partitions::Dict{Int, Dict{String, Any}},
    θ̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    ϵ_abs::Float64, ϵ_rel::Float64)
    N = max(1, length(shared_buses_sorted))
    norm_x_sq = 0.0
    for (_, sub) in partitions
        θbuf = snapshot_theta!(sub)
        pos = sub["bus_local_index"]
        @inbounds for b in sub["shared_buses_local"]
            norm_x_sq += (θbuf[pos[b]])^2
        end
    end
    norm_z_sq = sum((θ̄_vec[j])^2 for j in eachindex(θ̄_vec))
    return sqrt(N) * ϵ_abs + ϵ_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function eps_dual_bus_dc(θ̄_vec::Vector{Float64}, θ̄_prev_vec::Vector{Float64},
    ϵ_abs::Float64, ϵ_rel::Float64)
    N = max(1, length(θ̄_vec))
    return sqrt(N) * ϵ_abs +
           ϵ_rel * sqrt(sum((θ̄_vec[j] - θ̄_prev_vec[j])^2 for j in eachindex(θ̄_vec)))
end

function primal_residual_branch_dc(
    partitions::Dict{Int, Dict{String, Any}},
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

function eps_primal_branch_dc(partitions::Dict{Int, Dict{String, Any}},
    z_flow::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ϵ_abs::Float64, ϵ_rel::Float64)
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
    μ_incr = 3.0, μ_decr = 3.0, τ_incr = 1.05, τ_decr = 0.98)
    if r > μ_incr * s
        ρ[] *= τ_incr
    elseif s > μ_decr * r
        ρ[] *= τ_decr
    end
end

# ---------- Duais (flow) ----------
function update_duals_flow!(
    duals::Dict{Int, Dict{Int, Float64}},
    partitions::Dict{Int, Dict{String, Any}},
    z_flow::Dict{Int, Float64},
    z_flow_prev::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ρ::Float64,
    α::Union{Nothing, Float64} = nothing,
)
    for (z, sub) in partitions
        m = sub["model"]
        @inbounds for e in cut_edges
            l_local = get(sub["edge_to_local"], e, nothing)
            if l_local !== nothing
                p_val = _dc_flow_val_model(m, sub, l_local)
                if α === nothing
                    duals[z][e] = get(duals[z], e, 0.0) + ρ * (p_val - z_flow[e])
                else
                    z_prev = get(z_flow_prev, e, z_flow[e])
                    p_rel = α * p_val + (1 - α) * z_prev
                    duals[z][e] = get(duals[z], e, 0.0) + ρ * (p_rel - z_flow[e])
                end
            end
        end
    end
    return nothing
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

# =======================================================
# ADMM: BUS (com objetivo persistente)
# =======================================================
function run_admm(config::DistributedOPFModel, data)
    partitions, shared_buses = prepare_admm_partitions(data)
    @assert config.base_model == DC "run_admm (bus/DC) requer DC."

    for (_, sub) in partitions
        prepare_hotpaths!(sub)
    end
    shared_buses_sorted, b2pos = prepare_dense_maps(shared_buses)

    ρ = Ref(config.rho_theta)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    duals_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted))
        for z in keys(partitions)
    )
    θ̄_vec = zeros(Float64, length(shared_buses_sorted))
    θ̄_prev_vec = similar(θ̄_vec)

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rho_vec = Float64[]
    θ_avg_history = Dict{Int, Vector{Float64}}()
    θ_local_history = Dict{Int, Dict{Int, Vector{Float64}}}()

    optimizer_factory = resolve_optimizer_dc(config)

    # build + cache quadrático (bus)
    for (zone_id, sub) in partitions
        m = build_model(sub, DC; optimizer = optimizer_factory, build_only = true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["theta_var"] = m[:theta]
        sub["has_solution"] = false
        _set_model_silent!(m)
        for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
        _build_quad_bus_expr!(sub, ρ[])
        # objetivo inicial = base + quad_bus (linear = 0)
        @objective(m, Min, sub["obj_base"] + sub["quad_bus_expr"])
    end

    zone_ids = collect(keys(partitions))
    λ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted))
        for z in zone_ids
    )

    prev_rho = ρ[]  # para detectar mudança e reinstalar Q

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\nIteração $k (bus/DC)")

        if ρ[] != prev_rho
            for (_, sub) in partitions
                _build_quad_bus_expr!(sub, ρ[])
            end
            prev_rho = ρ[]
        end

        # ---------- resolve zonas em paralelo ----------
        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        @threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            λz = gather_zone_lambda!(λ_zone_buf[z], duals_dense[z])
            lin = _build_linear_bus_expr(sub, shared_buses_sorted, b2pos, θ̄_vec, λz, ρ[])
            # objetivo: base + quad_bus_cache + linear novo
            @objective(sub["model"], Min, sub["obj_base"] + sub["quad_bus_expr"] + lin)
            optimize!(sub["model"])

            m = sub["model"]
            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, sub["objective_params"])

            θbuf = snapshot_theta!(sub)
            θsnap = Dict{Int, Float64}()
            @inbounds for b in sub["shared_buses_local"]
                θsnap[b] = θbuf[sub["bus_local_index"][b]]
            end
            results[i] =
                (gen_cost_expr = gen_cost_expr, real_cost = real_cost, θsnap = θsnap)
        end

        for (i, z) in enumerate(zone_ids)
            sub = partitions[z]
            out = results[i]
            sub["has_solution"] = true
            sub["generation_cost_model_expr"] = out.gen_cost_expr
            sub["real_generation_cost"] = out.real_cost
            @inbounds for (b, θb) in out.θsnap
                push!(
                    get!(
                        get!(θ_local_history, z, Dict{Int, Vector{Float64}}()),
                        b,
                        Float64[],
                    ),
                    θb,
                )
            end
        end

        θ̄_prev_vec .= θ̄_vec
        compute_theta_consensus_vec!(θ̄_vec, partitions, shared_buses_sorted, b2pos)
        @inbounds for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_history, b, Float64[]), θ̄_vec[j])
        end

        for z in zone_ids
            update_duals_theta_vec!(duals_dense[z], partitions[z], θ̄_vec, θ̄_prev_vec, ρ[],
                config.alpha_overrelax, shared_buses_sorted, b2pos)
        end

        res_primal = primal_residual_bus_dc(partitions, θ̄_vec, shared_buses_sorted, b2pos)
        res_dual = dual_residual_bus_dc(θ̄_vec, θ̄_prev_vec)
        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        ϵ_pri =
            eps_primal_bus_dc(partitions, θ̄_vec, shared_buses_sorted, b2pos, ϵ_abs, ϵ_rel)
        ϵ_du = eps_dual_bus_dc(θ̄_vec, θ̄_prev_vec, ϵ_abs, ϵ_rel)
        push!(eps_primal_vec, ϵ_pri)
        push!(eps_dual_vec, ϵ_du)

        _FAST_USE_LOGS && println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        _FAST_USE_LOGS && println("‖s‖ = $res_dual   | ε_dual = $ϵ_du")

        push!(rho_vec, ρ[])

        if res_primal < ϵ_pri && res_dual < ϵ_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
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

# =======================================================
# ADMM: BRANCH (objetivo persistente)
# =======================================================
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

    for (_, sub) in partitions
        prepare_hotpaths!(sub)
        sub["zone_edges"] = _zone_cut_edges(sub, cut_edges)
    end

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
        m = build_model(sub, DC; optimizer = optimizer_factory, build_only = true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["theta_var"] = m[:theta]
        sub["has_solution"] = false
        _set_model_silent!(m)
        for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
        _build_quad_branch_expr!(sub, ρ[], sub["zone_edges"])
        @objective(m, Min, sub["obj_base"] + sub["quad_branch_expr"])
    end

    zone_ids = collect(keys(partitions))
    prev_rho = ρ[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-branch/DC] Iteração $k")

        if ρ[] != prev_rho
            for (_, sub) in partitions
                _build_quad_branch_expr!(sub, ρ[], sub["zone_edges"])
            end
            prev_rho = ρ[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        @threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            lin = _build_linear_branch_expr(sub, sub["zone_edges"], z_flow, duals[z], ρ[])
            @objective(sub["model"], Min, sub["obj_base"] + sub["quad_branch_expr"] + lin)
            optimize!(sub["model"])
            m = sub["model"]
            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, sub["objective_params"])
            results[i] = (gen_cost_expr = gen_cost_expr, real_cost = real_cost)
        end

        for (i, z) in enumerate(zone_ids)
            sub = partitions[z]
            out = results[i]
            sub["has_solution"] = true
            sub["generation_cost_model_expr"] = out.gen_cost_expr
            sub["real_generation_cost"] = out.real_cost
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

        _FAST_USE_LOGS && println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        _FAST_USE_LOGS && println("‖s‖ = $res_dual   | ε_dual = $ϵ_du")

        push!(rho_vec, ρ[])

        if res_primal < ϵ_pri && res_dual < ϵ_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
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

# =======================================================
# ADMM: HÍBRIDO (bus + branch) — objetivo persistente
# =======================================================
function run_admm_hybrid(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == DC "Hybrid ADMM/DC requer DC."

    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    for (_, sub) in partitions
        prepare_hotpaths!(sub)
        sub["zone_edges"] = _zone_cut_edges(sub, cut_edges)
    end

    shared_buses_sorted, b2pos = prepare_dense_maps(shared_buses)

    ρθ = Ref(config.rho_theta)
    ρe = Ref(config.rho_flow)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    θ̄_vec = zeros(Float64, length(shared_buses_sorted))
    θ̄_prev_vec = similar(θ̄_vec)

    z_flow = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    z_flow_prev = Dict{Int, Float64}()

    duals_theta = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
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
        m = build_model(sub, DC; optimizer = optimizer_factory, build_only = true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["theta_var"] = m[:theta]
        sub["has_solution"] = false
        _set_model_silent!(m)
        for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
        _build_quad_bus_expr!(sub, ρθ[])
        _build_quad_branch_expr!(sub, ρe[], sub["zone_edges"])
        @objective(m, Min, sub["obj_base"] + sub["quad_bus_expr"] + sub["quad_branch_expr"])
    end

    zone_ids = collect(keys(partitions))
    λθ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    prev_ρθ = ρθ[]
    prev_ρe = ρe[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-hybrid/DC] Iteração $k")

        # Reinstala quadráticas se rhos mudaram
        if ρθ[] != prev_ρθ
            for (_, sub) in partitions
                _build_quad_bus_expr!(sub, ρθ[])
            end
            prev_ρθ = ρθ[]
        end
        if ρe[] != prev_ρe
            for (_, sub) in partitions
                _build_quad_branch_expr!(sub, ρe[], sub["zone_edges"])
            end
            prev_ρe = ρe[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        @threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            λθ = gather_zone_lambda!(λθ_zone_buf[z], duals_theta[z])

            # Linear = soma BUS + BRANCH
            lin_bus =
                _build_linear_bus_expr(sub, shared_buses_sorted, b2pos, θ̄_vec, λθ, ρθ[])
            lin_branch = _build_linear_branch_expr(
                sub,
                sub["zone_edges"],
                z_flow,
                duals_flow[z],
                ρe[],
            )
            lin_total = lin_bus + lin_branch

            @objective(
                sub["model"],
                Min,
                sub["obj_base"] + sub["quad_bus_expr"] + sub["quad_branch_expr"] +
                lin_total
            )
            optimize!(sub["model"])

            m = sub["model"]
            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, sub["objective_params"])

            results[i] = (gen_cost_expr = gen_cost_expr, real_cost = real_cost)
        end

        for (i, z) in enumerate(zone_ids)
            sub = partitions[z]
            out = results[i]
            sub["has_solution"] = true
            sub["generation_cost_model_expr"] = out.gen_cost_expr
            sub["real_generation_cost"] = out.real_cost
        end

        θ̄_prev_vec .= θ̄_vec
        z_flow_prev = copy(z_flow)

        compute_theta_consensus_vec!(θ̄_vec, partitions, shared_buses_sorted, b2pos)
        for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_history, b, Float64[]), θ̄_vec[j])
        end
        compute_z_flow_consensus!(z_flow, partitions, cut_edges)

        for z in zone_ids
            update_duals_theta_vec!(
                duals_theta[z],
                partitions[z],
                θ̄_vec,
                θ̄_prev_vec,
                ρθ[],
                config.alpha_overrelax,
                shared_buses_sorted,
                b2pos,
            )
        end
        update_duals_flow!(
            duals_flow,
            partitions,
            z_flow,
            z_flow_prev,
            cut_edges,
            ρe[],
            config.alpha_overrelax,
        )

        res_pri_bus = primal_residual_bus_dc(partitions, θ̄_vec, shared_buses_sorted, b2pos)
        res_pri_br = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        res_primal = sqrt(res_pri_bus^2 + res_pri_br^2)

        res_du_bus = dual_residual_bus_dc(θ̄_vec, θ̄_prev_vec)
        res_du_br = dual_residual_branch_dc(z_flow, z_flow_prev)
        res_dual = sqrt(res_du_bus^2 + res_du_br^2)

        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        ϵ_pri = begin
            ϵb = eps_primal_bus_dc(partitions, θ̄_vec, shared_buses_sorted, b2pos, ϵ_abs, ϵ_rel)
            ϵe = eps_primal_branch_dc(partitions, z_flow, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(ϵb^2 + ϵe^2)
        end
        ϵ_du = begin
            ϵb = eps_dual_bus_dc(θ̄_vec, θ̄_prev_vec, ϵ_abs, ϵ_rel)
            ϵe = eps_dual_branch_dc(z_flow, z_flow_prev, ϵ_abs, ϵ_rel)
            sqrt(ϵb^2 + ϵe^2)
        end
        push!(eps_primal_vec, ϵ_pri)
        push!(eps_dual_vec, ϵ_du)

        _FAST_USE_LOGS && println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        _FAST_USE_LOGS && println("‖s‖ = $res_dual   | ε_dual = $ϵ_du")

        push!(rhoθ_vec, ρθ[])
        push!(rhoe_vec, ρe[])

        if res_primal < ϵ_pri && res_dual < ϵ_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
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
            )
        )(),
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
            )
        )(),
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
