# =============== DC Utils (FAST DROP-IN)  ===============
# Patches aplicados:
# (1.1) Snapshot vetorial de θ + índices precomputados
# (1.2) Penalidade construída de forma compacta e vetorial (sem lookups caros)
# (1.3) Consenso e duais com vetores densos

using LinearAlgebra
const _FAST_USE_LOGS = false

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
"""
    prepare_hotpaths!(sub)

Prepara índices posicionais e buffers para acesso contíguo a θ
e para trabalhar com vetores densos (sem Dict/set lookups).
"""
function prepare_hotpaths!(sub::Dict{String, Any})
    bus_ids = sub["bus"]["bus_i"]
    # mapa bus_id -> posição local na primeira dimensão de θ[:, t]
    bus_local_index = Dict{Int, Int}()
    @inbounds for (i, b) in enumerate(bus_ids)
        bus_local_index[b] = i
    end
    sub["bus_local_index"] = bus_local_index

    # posições locais (índices) apenas das compartilhadas desta zona
    local_shared = get(sub, "shared_buses_local", Int[])
    shared_pos_local = Vector{Int}(undef, length(local_shared))
    @inbounds for (k, b) in enumerate(local_shared)
        shared_pos_local[k] = bus_local_index[b]
    end
    sub["shared_pos_local"] = shared_pos_local

    # buffer contíguo para θ[:,1]
    sub["theta_buf"] = zeros(Float64, length(bus_ids))

    return nothing
end

"""
    snapshot_theta!(sub) -> Vector{Float64}

Preenche `theta_buf` com o valor de θ[:,1] de uma vez (batch value).
"""
# substitua esta função no seu arquivo FAST DROP-IN
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
# FLUXO DC
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
# PENALIDADES (1.2) – construção compacta, vetorial, sem lookups caros
# =======================================================
"""
    penalty_bus_dc_expr(m, sub, shared_buses_sorted, θ̄_vec, λ_vec, ρ)

Gera (uma única vez por iteração) a expressão de penalidade para as barras compartilhadas,
usando vetores densos (θ̄, λ) e os índices locais precomputados.
Retorna um `QuadExpr` compacto: ∑ [ λ*(θ-θ̄) + (ρ/2)*(θ-θ̄)^2 ].
"""
function penalty_bus_dc_expr(m::Model,
    sub::Dict{String, Any},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    θ̄_vec::Vector{Float64},
    λ_vec::Vector{Float64},
    ρ::Float64)
    # variável θ desta zona (DenseAxisArray indexado por bus_id)
    θvar = m[:theta]

    pen = JuMP.QuadExpr()
    @inbounds for b in sub["shared_buses_local"]
        j = b2pos[b]         # posição densa global
        θb = θvar[b, 1]       # <<< usar o ID da barra, NÃO posição local
        θ̄ = θ̄_vec[j]
        λ = λ_vec[j]

        # (ρ/2)*(θ - θ̄)^2 + λ*(θ-θ̄) = (ρ/2)θ^2 + (λ - ρ*θ̄)θ + const
        JuMP.add_to_expression!(pen, 0.5 * ρ, θb, θb)       # termo quadrático
        JuMP.add_to_expression!(pen, (λ - ρ * θ̄), θb)      # termo linear
    end
    return pen
end

function penalty_branch_dc_expr(
    m::Model,
    sub::Dict{String,Any},
    cut_edges::Vector{Int},
    z_flow::Dict{Int,Float64},
    λe_zone::Dict{Int,Float64},
    ρ::Float64,
)
    pen = JuMP.QuadExpr()
    @inbounds for e in cut_edges
        l_local = get(sub["edge_to_local"], e, nothing)
        if l_local === nothing
            continue
        end
        # p_loc é um AffExpr: (θf - θt)/x
        p_loc = _dc_flow_expr(m, sub, l_local)

        z  = z_flow[e]
        λe = get(λe_zone, e, 0.0)

        # termo afim: λe*(p_loc - z)
        # (há método add_to_expression!(QuadExpr, Number, AffExpr))
        JuMP.add_to_expression!(pen, λe, p_loc - z)

        # termo quadrático: (ρ/2)*(p_loc - z)^2  (AffExpr^2 -> QuadExpr)
        pen += 0.5 * ρ * (p_loc - z)^2
    end
    return pen
end


# =======================================================
# SOLVER POR ZONA (com objetivos compactos)
# =======================================================
function _solve_zone_bus!(sub::Dict{String, Any},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    θ̄_vec::Vector{Float64},
    λ_vec_zone::Vector{Float64},
    ρ::Float64)
    m = sub["model"]
    pen = penalty_bus_dc_expr(m, sub, shared_buses_sorted, b2pos, θ̄_vec, λ_vec_zone, ρ)
    @objective(m, Min, sub["obj_base"] + pen)
    optimize!(m)

    gen_cost_expr = value(m[:generation_cost])
    real_cost = compute_real_generation_cost(m, sub["objective_params"])

    # Snapshot local das compartilhadas (para histórico)
    θbuf = snapshot_theta!(sub)
    θsnap = Dict{Int, Float64}()
    @inbounds for b in sub["shared_buses_local"]
        θsnap[b] = θbuf[sub["bus_local_index"][b]]
    end
    return (gen_cost_expr = gen_cost_expr, real_cost = real_cost, θsnap = θsnap)
end

function _solve_zone_branch!(sub::Dict{String, Any},
    cut_edges::Vector{Int},
    z_flow::Dict{Int, Float64},
    λe_zone::Dict{Int, Float64},
    ρ::Float64)
    m = sub["model"]
    pen = penalty_branch_dc_expr(m, sub, cut_edges, z_flow, λe_zone, ρ)
    @objective(m, Min, sub["obj_base"] + pen)
    optimize!(m)
    gen_cost_expr = value(m[:generation_cost])
    real_cost = compute_real_generation_cost(m, sub["objective_params"])
    return (gen_cost_expr = gen_cost_expr, real_cost = real_cost)
end

function _solve_zone_hybrid!(sub::Dict{String, Any},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    θ̄_vec::Vector{Float64}, λθ_zone::Vector{Float64}, ρθ::Float64,
    cut_edges::Vector{Int},
    z_flow::Dict{Int, Float64}, λe_zone::Dict{Int, Float64}, ρe::Float64)
    m = sub["model"]
    pen_bus = penalty_bus_dc_expr(m, sub, shared_buses_sorted, b2pos, θ̄_vec, λθ_zone, ρθ)
    pen_flow = penalty_branch_dc_expr(m, sub, cut_edges, z_flow, λe_zone, ρe)
    @objective(m, Min, sub["obj_base"] + pen_bus + pen_flow)
    optimize!(m)
    gen_cost_expr = value(m[:generation_cost])
    real_cost = compute_real_generation_cost(m, sub["objective_params"])
    return (gen_cost_expr = gen_cost_expr, real_cost = real_cost)
end

# =======================================================
# CONSENSO / DUAIS DENSOS (1.3)
# =======================================================
"""
    prepare_dense_maps(shared_buses) -> shared_buses_sorted, b2pos

Ordena as compartilhadas e retorna mapa bus_id -> posição no vetor denso.
"""
function prepare_dense_maps(shared_buses::Vector{Int})
    shared_buses_sorted = sort(shared_buses)
    b2pos = Dict{Int, Int}(b => i for (i, b) in enumerate(shared_buses_sorted))
    return shared_buses_sorted, b2pos
end

"""
    compute_theta_consensus_vec!(θ̄_vec, partitions, shared_buses_sorted, b2pos)

Consenso de θ usando snapshots vetoriais e vetores densos.
"""
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
        if counts[j] != 0
            θ̄_vec[j] /= counts[j]
        end
    end
end

"""
    compute_z_flow_consensus!(z_flow, partitions, cut_edges)

(mantido com mesmo contrato; já era relativamente barato)
"""
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
        if n > 0
            z_flow[e] = sums[e] / n
        end
    end
end

"""
    update_duals_theta_vec!(dual_vec_zone, sub, θ̄_vec, θ̄_prev_vec, ρ, α, shared_buses_sorted, b2pos)

Atualiza duais λ (vetor) da zona para barras compartilhadas.
"""
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
            θb = θbuf[pos[b]]
            dual_vec_zone[j] += ρ * (θb - θ̄_vec[j])
        end
    else
        αv = something(α, 1.0)
        @inbounds for b in sub["shared_buses_local"]
            j = b2pos[b]
            θb = θbuf[pos[b]]
            θrel = αv * θb + (1 - αv) * θ̄_prev_vec[j]
            dual_vec_zone[j] += ρ * (θrel - θ̄_vec[j])
        end
    end
end

"""
    gather_zone_lambda!(λ_vec_zone, duals_theta_dense[z])

Copia (view) do vetor de duais da zona (já denso).
"""
@inline function gather_zone_lambda!(λ_vec_zone::Vector{Float64},
    duals_dense_zone::Vector{Float64})
    λ_vec_zone .= duals_dense_zone
    return λ_vec_zone
end

# =======================================================
# RESÍDUOS / TOLERÂNCIAS (usam snapshots quando possível)
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

# ---------- Duais (flow) com over-relaxation opcional (drop-in) ----------
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
                    p_rel  = α * p_val + (1 - α) * z_prev
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
# ADMM: BUS
# =======================================================
function run_admm(config::DistributedOPFModel, data)
    partitions, shared_buses = prepare_admm_partitions(data)
    @assert config.base_model == DC "run_admm (bus/DC) requer DC."

    # preparar hotpaths por zona (1.1)
    for (_, sub) in partitions
        prepare_hotpaths!(sub)
    end

    # mapas densos globais (1.3)
    shared_buses_sorted, b2pos = prepare_dense_maps(shared_buses)

    ρ = Ref(config.rho_theta)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    # duais densos por zona (um vetor por zona, alinhado com shared_buses_sorted) (1.3)
    duals_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted))
        for z in keys(partitions)
    )
    # médias globais densas (θ̄) e cópia anterior
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

    for (zone_id, sub) in partitions
        m = build_model(sub, DC; optimizer = optimizer_factory, build_only = true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["theta_var"] = m[:theta]
        sub["has_solution"] = false
        _set_model_silent!(m)
        @inbounds for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
    end

    zone_ids = collect(keys(partitions))
    λ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted))
        for z in zone_ids
    )

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\nIteração $k (bus/DC)")

        # ---------- resolve zonas em paralelo ----------
        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        @threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            # pega λ denso da zona
            λz = gather_zone_lambda!(λ_zone_buf[z], duals_dense[z])
            results[i] = _solve_zone_bus!(sub,
                shared_buses_sorted, b2pos,
                θ̄_vec, λz, ρ[])
        end

        # ---------- aplica resultados (sequencial) ----------
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

        # atualiza duais (vetorial, por zona)
        for z in zone_ids
            update_duals_theta_vec!(duals_dense[z], partitions[z],
                θ̄_vec, θ̄_prev_vec, ρ[],
                config.alpha_overrelax,
                shared_buses_sorted, b2pos)
        end

        # resíduos / tolerâncias
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

        # println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        # println("‖s‖ = $res_dual   | ε_dual = $ϵ_du")

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
# ADMM: BRANCH
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

    # preparar hotpaths
    for (_, sub) in partitions
        prepare_hotpaths!(sub)
        # também deixe marcado quais são as arestas do corte nesta zona
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
        @inbounds for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
    end

    zone_ids = collect(keys(partitions))

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-branch/DC] Iteração $k")

        # ---------- resolve zonas em paralelo ----------
        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        @threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            results[i] = _solve_zone_branch!(sub, sub["zone_edges"], z_flow, duals[z], ρ[])
        end

        # ---------- aplica resultados (sequencial) ----------
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
# ADMM: HÍBRIDO
# =======================================================
function run_admm_hybrid(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == DC "Hybrid ADMM/DC requer DC."

    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    # hotpaths + zone_edges
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
        z => zeros(Float64, length(shared_buses_sorted))
        for z in keys(partitions)
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
        @inbounds for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
    end

    zone_ids = collect(keys(partitions))
    λθ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted))
        for z in zone_ids
    )

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-hybrid/DC] Iteração $k")

        # ---------- resolve zonas em paralelo ----------
        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        @threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            λθ = gather_zone_lambda!(λθ_zone_buf[z], duals_theta[z])
            results[i] = _solve_zone_hybrid!(sub,
                shared_buses_sorted, b2pos,
                θ̄_vec, λθ, ρθ[],
                sub["zone_edges"], z_flow, duals_flow[z], ρe[])
        end

        # ---------- aplica resultados (sequencial) ----------
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
        @inbounds for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_history, b, Float64[]), θ̄_vec[j])
        end
        compute_z_flow_consensus!(z_flow, partitions, cut_edges)

        # duais
        for z in zone_ids
            update_duals_theta_vec!(duals_theta[z], partitions[z],
                θ̄_vec, θ̄_prev_vec, ρθ[],
                config.alpha_overrelax,
                shared_buses_sorted, b2pos)
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

        # resíduos
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
                p = Plot(); p(;
                    x = 1:length(rhoθ_vec),
                    y = rhoθ_vec,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρθ",
                ); p.layout.title = "Evolução de ρθ"; p.layout.xaxis.title = "Iterações"; p.layout.yaxis.title = "ρθ"; p.layout.yaxis.type = "log"; p
            )
        )(),
        "plot_rho_theta_hybrid.html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot(); p(;
                    x = 1:length(rhoe_vec),
                    y = rhoe_vec,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρe",
                ); p.layout.title = "Evolução de ρe"; p.layout.xaxis.title = "Iterações"; p.layout.yaxis.title = "ρe"; p.layout.yaxis.type = "log"; p
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
