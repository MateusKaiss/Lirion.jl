# ================= SOCP Utils =================
# ---- Otimizador (default SCS) ----
function resolve_socp_optimizer(config)
    if hasproperty(config, :optimizer) && getproperty(config, :optimizer) !== nothing
        return getproperty(config, :optimizer)
    else
        return SCS.Optimizer
    end
end

safe_set_optimizer_attr!(m, key, val) = (
    try
        set_optimizer_attribute(m, key, val)
    catch
        nothing
    end
)

"Configura tolerâncias/accel do SCS."
function set_socp_defaults!(m; eps_abs = 1e-6, eps_rel = 1e-6, max_iters = 200_000,
    acceleration_lookback = 10, acceleration_interval = 10,
    alpha = 1.8, normalize = true, scale = 1.0)
    safe_set_optimizer_attr!(m, "eps_abs", eps_abs)
    safe_set_optimizer_attr!(m, "eps_rel", eps_rel)
    safe_set_optimizer_attr!(m, "max_iters", max_iters)
    safe_set_optimizer_attr!(m, "acceleration_lookback", acceleration_lookback)
    safe_set_optimizer_attr!(m, "acceleration_interval", acceleration_interval)
    safe_set_optimizer_attr!(m, "alpha", alpha)
    safe_set_optimizer_attr!(m, "normalize", normalize)
    safe_set_optimizer_attr!(m, "scale", scale)
    return nothing
end

# ---- Penalidades ADMM ----
# BUS (V^2)
function add_penalty_bus_socp!(m, sub, shared_buses::Vector{Int},
    v2_avg::Dict{Int, Float64}, duals_v2_zone::Dict{Int, Float64},
    ρV2::Float64)
    pen = @expression(m, 0.0)
    local_shared = get(sub, "shared_buses_local", shared_buses)
    @inbounds for b in local_shared
        V2 = m[:v_magnitude_sq][b, 1]
        V2̄ = Float64(get(v2_avg, b, 1.0))
        λ = Float64(get(duals_v2_zone, b, 0.0))
        pen += λ * (V2 - V2̄) + (float(ρV2) / 2) * (V2 - V2̄)^2
    end
    return pen
end

# BRANCH (R,T)
function add_penalty_branch_rt_socp!(m, sub, zone_edges::Vector{Int},
    zR::Dict{Int, Float64}, zT::Dict{Int, Float64},
    dualsR_zone::Dict{Int, Float64}, dualsT_zone::Dict{Int, Float64},
    ρR::Float64, ρT::Float64)
    pen = @expression(m, 0.0)
    @inbounds for e in zone_edges
        l_local = get(sub["edge_to_local"], e, nothing)
        @assert l_local !== nothing "edge $e não está na partição"
        Rloc, Tloc = _rt_from_exprs_socp(m, l_local)
        λR = Float64(get(dualsR_zone, e, 0.0))
        λT = Float64(get(dualsT_zone, e, 0.0))
        zRe = Float64(get(zR, e, 0.0))
        zTe = Float64(get(zT, e, 0.0))
        pen += λR * (Rloc - zRe) + (float(ρR) / 2) * (Rloc - zRe)^2
        pen += λT * (Tloc - zTe) + (float(ρT) / 2) * (Tloc - zTe)^2
    end
    return pen
end

# BRANCH (P,Q)
function add_penalty_branch_pq_socp!(m, sub, zone_edges::Vector{Int},
    zP::Dict{Int, Float64}, zQ::Dict{Int, Float64},
    dualsP_zone::Dict{Int, Float64}, dualsQ_zone::Dict{Int, Float64},
    ρP::Float64, ρQ::Float64)
    pen = @expression(m, 0.0)
    @inbounds for e in zone_edges
        l_local = get(sub["edge_to_local"], e, nothing)
        @assert l_local !== nothing "edge $e não está na partição"
        P_loc, Q_loc = _socp_flow_from_exprs_existing(m, l_local)
        λP = Float64(get(dualsP_zone, e, 0.0))
        λQ = Float64(get(dualsQ_zone, e, 0.0))
        zPe = Float64(get(zP, e, 0.0))
        zQe = Float64(get(zQ, e, 0.0))
        pen += λP * (P_loc - zPe) + (float(ρP) / 2) * (P_loc - zPe)^2
        pen += λQ * (Q_loc - zQe) + (float(ρQ) / 2) * (Q_loc - zQe)^2
    end
    return pen
end

# ---- Consensos (★ evitar value.(matriz inteira)) ----
function compute_v2_consensus!(
    v2_avg::Dict{Int, Float64},
    partitions,
    shared_buses::Vector{Int},
)
    sums = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()
    for (_, sub) in partitions
        # ler só as barras compartilhadas locais
        local_shared = get(sub, "shared_buses_local", shared_buses)
        @inbounds for b in local_shared
            v2b = value(sub["model"][:v_magnitude_sq][b, 1])
            sums[b] = get(sums, b, 0.0) + v2b
            cnts[b] = get(cnts, b, 0) + 1
        end
    end
    @inbounds for b in shared_buses
        n = get(cnts, b, 0)
        if n > 0
            v2_avg[b] = sums[b] / n
        end
    end
end

function compute_rt_consensus!(zR::Dict{Int, Float64}, zT::Dict{Int, Float64},
    partitions, cut_edges::Vector{Int})
    sumsR = Dict{Int, Float64}()
    sumsT = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()
    for (_, sub) in partitions
        @inbounds for e in cut_edges
            l_local = get(sub["edge_to_local"], e, nothing)
            if l_local !== nothing
                R, T = _rt_from_exprs_socp(sub["model"], l_local)
                sumsR[e] = get(sumsR, e, 0.0) + value(R)
                sumsT[e] = get(sumsT, e, 0.0) + value(T)
                cnts[e] = get(cnts, e, 0) + 1
            end
        end
    end
    @inbounds for e in cut_edges
        n = get(cnts, e, 0)
        if n > 0
            zR[e] = sumsR[e] / n
            zT[e] = sumsT[e] / n
        end
    end
end

function compute_pq_consensus!(zP::Dict{Int, Float64}, zQ::Dict{Int, Float64},
    partitions, cut_edges::Vector{Int})
    sumsP = Dict{Int, Float64}()
    sumsQ = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()
    for (_, sub) in partitions
        @inbounds for e in cut_edges
            l_local = get(sub["edge_to_local"], e, nothing)
            if l_local !== nothing
                P, Q = _socp_flow_from_exprs_existing(sub["model"], l_local)
                sumsP[e] = get(sumsP, e, 0.0) + value(P)
                sumsQ[e] = get(sumsQ, e, 0.0) + value(Q)
                cnts[e] = get(cnts, e, 0) + 1
            end
        end
    end
    @inbounds for e in cut_edges
        n = get(cnts, e, 0)
        if n > 0
            zP[e] = sumsP[e] / n
            zQ[e] = sumsQ[e] / n
        end
    end
end

# ---- Duais ----
function update_duals_bus_socp!(duals_v2, partitions, v2_avg, shared_buses, ρV2)
    for (z, sub) in partitions
        local_shared = get(sub, "shared_buses_local", shared_buses)
        @inbounds for b in local_shared
            v2b = value(sub["model"][:v_magnitude_sq][b, 1])
            duals_v2[z][b] = get(duals_v2[z], b, 0.0) + ρV2 * (v2b - get(v2_avg, b, 0.0))
        end
    end
end

function update_duals_branch_rt_socp!(dualsR, dualsT, partitions, zR, zT, cut_edges, ρR, ρT)
    for (z, sub) in partitions
        @inbounds for e in cut_edges
            l = get(sub["edge_to_local"], e, nothing)
            if l !== nothing
                R, T = _rt_from_exprs_socp(sub["model"], l)
                dualsR[z][e] = get(dualsR[z], e, 0.0) + ρR * (value(R) - zR[e])
                dualsT[z][e] = get(dualsT[z], e, 0.0) + ρT * (value(T) - zT[e])
            end
        end
    end
end

function update_duals_branch_pq_socp!(dualsP, dualsQ, partitions, zP, zQ, cut_edges, ρP, ρQ)
    for (z, sub) in partitions
        @inbounds for e in cut_edges
            l = get(sub["edge_to_local"], e, nothing)
            if l !== nothing
                P, Q = _socp_flow_from_exprs_existing(sub["model"], l)
                dualsP[z][e] = get(dualsP[z], e, 0.0) + ρP * (value(P) - zP[e])
                dualsQ[z][e] = get(dualsQ[z], e, 0.0) + ρQ * (value(Q) - zQ[e])
            end
        end
    end
end

# ---- Resíduos / tolerâncias (★ evitar value.(...)) ----
function primal_residual_bus_socp(partitions, v2_avg, shared_buses)
    acc = 0.0
    for (_, sub) in partitions
        local_shared = get(sub, "shared_buses_local", shared_buses)
        @inbounds for b in local_shared
            acc += (value(sub["model"][:v_magnitude_sq][b, 1]) - get(v2_avg, b, 0.0))^2
        end
    end
    return sqrt(acc)
end

dual_residual_bus_socp(v2_avg, v2_prev, shared_buses) =
    sqrt(sum((get(v2_avg, b, 0.0) - get(v2_prev, b, 0.0))^2 for b in shared_buses))

function eps_primal_bus_socp(partitions, v2_avg, shared_buses, ϵ_abs, ϵ_rel)
    dim = 0
    for (_, sub) in partitions
        local_shared = get(sub, "shared_buses_local", shared_buses)
        dim += length(local_shared)
    end
    z_norm_sq = 0.0
    @inbounds for b in shared_buses
        z_norm_sq += get(v2_avg, b, 0.0)^2
    end
    return sqrt(max(1, dim)) * ϵ_abs + ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_bus_socp(v2_avg, v2_prev, shared_buses, ϵ_abs, ϵ_rel)
    dim = max(1, length(shared_buses))
    s = dual_residual_bus_socp(v2_avg, v2_prev, shared_buses)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
end

function primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
    acc = 0.0
    for (_, sub) in partitions
        @inbounds for e in cut_edges
            l = get(sub["edge_to_local"], e, nothing)
            if l !== nothing
                R, T = _rt_from_exprs_socp(sub["model"], l)
                acc += (value(R) - zR[e])^2 + (value(T) - zT[e])^2
            end
        end
    end
    return sqrt(acc)
end

dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges) =
    sqrt(
        sum(
            (get(zR, e, 0.0) - get(zR_prev, e, 0.0))^2 +
            (get(zT, e, 0.0) - get(zT_prev, e, 0.0))^2 for e in cut_edges
        ),
    )

function eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, ϵ_abs, ϵ_rel)
    dim = 0
    for (_, sub) in partitions
        dim += 2 * length(_zone_cut_edges(sub, cut_edges))
    end
    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(zR, e, 0.0)^2 + get(zT, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * ϵ_abs + ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, ϵ_abs, ϵ_rel)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
end

function primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
    acc = 0.0
    for (_, sub) in partitions
        @inbounds for e in cut_edges
            l = get(sub["edge_to_local"], e, nothing)
            if l !== nothing
                P, Q = _socp_flow_from_exprs_existing(sub["model"], l)
                acc += (value(P) - zP[e])^2 + (value(Q) - zQ[e])^2
            end
        end
    end
    return sqrt(acc)
end

dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges) =
    sqrt(
        sum(
            (get(zP, e, 0.0) - get(zP_prev, e, 0.0))^2 +
            (get(zQ, e, 0.0) - get(zQ_prev, e, 0.0))^2 for e in cut_edges
        ),
    )

function eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
    dim = 0
    for (_, sub) in partitions
        dim += 2 * length(_zone_cut_edges(sub, cut_edges))
    end
    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(zP, e, 0.0)^2 + get(zQ, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * ϵ_abs + ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
end

function adaptive_rho!(ρ::Base.RefValue{Float64}, r::Float64, s::Float64;
    μ_incr = 2.0, μ_decr = 2.0, τ_incr = 1.05, τ_decr = 0.98)
    if r > μ_incr * s
        ρ[] *= τ_incr
    elseif s > μ_decr * r
        ρ[] *= τ_decr
    end
end

# ---- Plots (wrappers) ----
function save_socp_v2rt_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
    rhoV2_hist, rhoR_hist, rhoT_hist; stem = "admm_socp_v2rt")
    plot_admm_convergence(r_hist, s_hist, eps_r_hist, eps_s_hist,
        "plot_convergencia_$(stem).html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoV2_hist),
                    y = rhoV2_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρ(V²)",
                );
                p.layout.title = "Evolução de ρ(V²)";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρ(V²)";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_v2_$(stem).html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot(); p(;
                    x = 1:length(rhoR_hist),
                    y = rhoR_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρR",
                ); p.layout.title = "Evolução de ρR"; p.layout.xaxis.title = "Iterações"; p.layout.yaxis.title = "ρR"; p.layout.yaxis.type = "log"; p
            )
        )(),
        "plot_rhoR_$(stem).html")
    return PlotlyLight.save(
        (
            () -> (
                p = Plot(); p(;
                    x = 1:length(rhoT_hist),
                    y = rhoT_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρT",
                ); p.layout.title = "Evolução de ρT"; p.layout.xaxis.title = "Iterações"; p.layout.yaxis.title = "ρT"; p.layout.yaxis.type = "log"; p
            )
        )(),
        "plot_rhoT_$(stem).html")
end

function save_socp_branch_pq_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
    rhoP_hist, rhoQ_hist; stem = "admm_socp_branch")
    plot_admm_convergence(r_hist, s_hist, eps_r_hist, eps_s_hist,
        "plot_convergencia_$(stem).html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot(); p(;
                    x = 1:length(rhoP_hist),
                    y = rhoP_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρP",
                ); p.layout.title = "Evolução de ρP"; p.layout.xaxis.title = "Iterações"; p.layout.yaxis.title = "ρP"; p.layout.yaxis.type = "log"; p
            )
        )(),
        "plot_rhoP_$(stem).html")
    return PlotlyLight.save(
        (
            () -> (
                p = Plot(); p(;
                    x = 1:length(rhoQ_hist),
                    y = rhoQ_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρQ",
                ); p.layout.title = "Evolução de ρQ"; p.layout.xaxis.title = "Iterações"; p.layout.yaxis.title = "ρQ"; p.layout.yaxis.type = "log"; p
            )
        )(),
        "plot_rhoQ_$(stem).html")
end

function save_socp_hybrid_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
    rhoV2_hist, rhoP_hist, rhoQ_hist, rhoR_hist, rhoT_hist;
    stem = "admm_socp_hybrid")
    plot_admm_convergence(r_hist, s_hist, eps_r_hist, eps_s_hist,
        "plot_convergencia_$(stem).html")
    for (lab, hist) in [
        ("ρ(V²)", rhoV2_hist),
        ("ρP", rhoP_hist),
        ("ρQ", rhoQ_hist),
        ("ρR", rhoR_hist),
        ("ρT", rhoT_hist),
    ]
        fname = replace(lab, "/" => "") |> x -> lowercase(replace(x, "ρ" => "rho"))
        PlotlyLight.save(
            (
                () -> (
                    p = Plot(); p(;
                        x = 1:length(hist),
                        y = hist,
                        type = :scatter,
                        mode = "lines+markers",
                        name = lab,
                    ); p.layout.title = "Evolução de $(lab)"; p.layout.xaxis.title = "Iterações"; p.layout.yaxis.title =
                        lab; p.layout.yaxis.type = "log"; p
                )
            )(),
            "plot_$(fname)_$(stem).html")
    end
end

# ================== RUNNERS ==================

function run_admm_socp(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == SOCP "ADMM-SOCP: base_model deve ser SOCP."

    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    # estado
    ρV2 = Ref(config.rho_v)
    ρR = Ref(config.rho_flow_p)
    ρT = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    # duais / consensos
    duals_v2 =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsR =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsT =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    v2_avg = Dict{Int, Float64}()
    v2_prev = Dict{Int, Float64}()
    zR = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zR_prev = Dict{Int, Float64}()
    zT = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zT_prev = Dict{Int, Float64}()

    # histórico
    r_hist = Float64[]
    s_hist = Float64[]
    eps_r_hist = Float64[]
    eps_s_hist = Float64[]
    rhoV2_hist = Float64[]
    rhoR_hist = Float64[]
    rhoT_hist = Float64[]

    optimizer_factory = resolve_socp_optimizer(config)

    # PRÉ-BUILD (★: preparar shared_buses_local, zone_edges, silent)
    for (zone_id, sub) in partitions
        m = build_model(sub, SOCP; optimizer = optimizer_factory)
        _set_model_silent!(m)
        set_socp_defaults!(m; eps_abs = 1e-6, eps_rel = 1e-6, max_iters = 200_000)

        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["v2_var"] = m[:v_magnitude_sq]
        sub["zone_edges"] = _zone_cut_edges(sub, cut_edges)
        sub["has_solution"] = false

        busset = _bus_set(sub)
        sub["shared_buses_local"] = [b for b in shared_buses if b in busset]

        # warm start
        @inbounds for I in keys(sub["v2_var"])
            set_start_value(sub["v2_var"][I], 1.0)
        end
    end

    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98

    for k in 1:config.max_iter
        println("\n[ADMM-SOCP (V²+R,T)] Iteração $k")

        # 1) subproblemas locais
        for (zone_id, sub) in partitions
            m = sub["model"]
            # (warm-start comentado – ligue se quiser)
            # if sub["has_solution"] === true
            #     V2_old = JuMP.value.(sub["v2_var"])
            #     @inbounds for I in keys(sub["v2_var"])
            #         set_start_value(sub["v2_var"][I], V2_old[I])
            #     end
            # end

            pen = add_penalty_bus_socp!(
                m,
                sub,
                shared_buses,
                v2_avg,
                duals_v2[zone_id],
                ρV2[],
            )
            pen += add_penalty_branch_rt_socp!(
                m,
                sub,
                sub["zone_edges"],
                zR,
                zT,
                dualsR[zone_id],
                dualsT[zone_id],
                ρR[],
                ρT[],
            )
            @objective(m, Min, sub["obj_base"] + pen)
            optimize!(m)
            sub["has_solution"] = true

            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] =
                compute_real_generation_cost(m, sub["objective_params"])
        end

        # 2) consensos
        v2_prev = isempty(v2_avg) ? v2_avg : copy(v2_avg)
        zR_prev = copy(zR)
        zT_prev = copy(zT)
        compute_v2_consensus!(v2_avg, partitions, shared_buses)
        compute_rt_consensus!(zR, zT, partitions, cut_edges)

        # 3) duais
        update_duals_bus_socp!(duals_v2, partitions, v2_avg, shared_buses, ρV2[])
        update_duals_branch_rt_socp!(
            dualsR,
            dualsT,
            partitions,
            zR,
            zT,
            cut_edges,
            ρR[],
            ρT[],
        )

        # 4) resíduos + tolerâncias
        r_bus = primal_residual_bus_socp(partitions, v2_avg, shared_buses)
        r_rt = primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
        r = sqrt(r_bus^2 + r_rt^2)

        s_bus = dual_residual_bus_socp(v2_avg, v2_prev, shared_buses)
        s_rt = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
        s = sqrt(s_bus^2 + s_rt^2)

        push!(r_hist, r)
        push!(s_hist, s)

        ε_pri = begin
            εb = eps_primal_bus_socp(partitions, v2_avg, shared_buses, ϵ_abs, ϵ_rel)
            εe = eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εe^2)
        end
        ε_du = begin
            εb = eps_dual_bus_socp(v2_avg, v2_prev, shared_buses, ϵ_abs, ϵ_rel)
            εe = eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εe^2)
        end
        push!(eps_r_hist, ε_pri)
        push!(eps_s_hist, ε_du)

        println("‖r‖ = $r | ε_pri = $ε_pri")
        println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoV2_hist, ρV2[])
        push!(rhoR_hist, ρR[])
        push!(rhoT_hist, ρT[])

        if r < ε_pri && s < ε_du
            println("Stopping criterion met.")
            break
        end

        # 6) atualização conjunta de ρ (V², R, T)
        if r > μ_incr * s
            ρV2[] *= τ_incr
            ρR[] *= τ_incr
            ρT[] *= τ_incr
        elseif s > μ_decr * r
            ρV2[] *= τ_decr
            ρR[] *= τ_decr
            ρT[] *= τ_decr
        end
    end

    save_socp_v2rt_plots(
        r_hist,
        s_hist,
        eps_r_hist,
        eps_s_hist,
        rhoV2_hist,
        rhoR_hist,
        rhoT_hist;
        stem = "admm_socp_no_theta",
    )

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-SOCP (V²+R,T)] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρ(V²)_final = $(ρV2[]) | ρR_final = $(ρR[]) | ρT_final = $(ρT[])")

    return total_real_cost, partitions
end

function run_admm_socp_branch(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == SOCP "ADMM-SOCP-branch: base_model deve ser SOCP."

    partitions, _ = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    dualsP =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ_prev = Dict{Int, Float64}()

    r_hist = Float64[]
    s_hist = Float64[]
    eps_r_hist = Float64[]
    eps_s_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]

    optimizer_factory = resolve_socp_optimizer(config)

    for (zone_id, sub) in partitions
        m = build_model(sub, SOCP; optimizer = optimizer_factory)
        _set_model_silent!(m)
        set_socp_defaults!(m; eps_abs = 1e-5, eps_rel = 1e-5, max_iters = 150_000,
            acceleration_lookback = 10, acceleration_interval = 10,
            alpha = 1.8, normalize = true, scale = 1.0)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["zone_edges"] = _zone_cut_edges(sub, cut_edges)
        sub["has_solution"] = false
    end

    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98

    for k in 1:config.max_iter
        println("\n[ADMM-SOCP-branch (P,Q)] Iteração $k")

        for (zone_id, sub) in partitions
            m = sub["model"]
            pen = add_penalty_branch_pq_socp!(
                m,
                sub,
                sub["zone_edges"],
                zP,
                zQ,
                dualsP[zone_id],
                dualsQ[zone_id],
                ρP[],
                ρQ[],
            )
            @objective(m, Min, sub["obj_base"] + pen)
            optimize!(m)
            sub["has_solution"] = true

            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] =
                compute_real_generation_cost(m, sub["objective_params"])
        end

        zP_prev = copy(zP)
        zQ_prev = copy(zQ)
        compute_pq_consensus!(zP, zQ, partitions, cut_edges)
        update_duals_branch_pq_socp!(
            dualsP,
            dualsQ,
            partitions,
            zP,
            zQ,
            cut_edges,
            ρP[],
            ρQ[],
        )

        r = primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
        s = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
        push!(r_hist, r)
        push!(s_hist, s)

        ε_pri = eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
        ε_du = eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
        push!(eps_r_hist, ε_pri)
        push!(eps_s_hist, ε_du)

        println("‖r‖ = $r | ε_pri = $ε_pri")
        println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoP_hist, ρP[])
        push!(rhoQ_hist, ρQ[])

        if r < ε_pri && s < ε_du
            println("Stopping criterion met.")
            break
        end

        if r > μ_incr * s
            ρP[] *= τ_incr
            ρQ[] *= τ_incr
        elseif s > μ_decr * r
            ρP[] *= τ_decr
            ρQ[] *= τ_decr
        end
    end

    save_socp_branch_pq_plots(
        r_hist,
        s_hist,
        eps_r_hist,
        eps_s_hist,
        rhoP_hist,
        rhoQ_hist;
        stem = "admm_socp_branch",
    )

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-SOCP-branch (P,Q)] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρP_final = $(ρP[]) | ρQ_final = $(ρQ[])")

    return total_real_cost, partitions
end

function run_admm_socp_hybrid(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == SOCP "ADMM-SOCP-hybrid: base_model deve ser SOCP."
    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    ρV2 = Ref(config.rho_v)
    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ρR = Ref(config.rho_flow_p)
    ρT = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    duals_v2 =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsP =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsR =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsT =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    v2_avg = Dict{Int, Float64}()
    v2_prev = Dict{Int, Float64}()
    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ_prev = Dict{Int, Float64}()
    zR = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zR_prev = Dict{Int, Float64}()
    zT = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zT_prev = Dict{Int, Float64}()

    r_hist = Float64[]
    s_hist = Float64[]
    eps_r_hist = Float64[]
    eps_s_hist = Float64[]
    rhoV2_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]
    rhoR_hist = Float64[]
    rhoT_hist = Float64[]

    optimizer_factory = resolve_socp_optimizer(config)

    # PRÉ-BUILD (★)
    for (zone_id, sub) in partitions
        m = build_model(sub, SOCP; optimizer = optimizer_factory)
        _set_model_silent!(m)
        set_socp_defaults!(m; eps_abs = 1e-5, eps_rel = 1e-5, max_iters = 150_000,
            acceleration_lookback = 10, acceleration_interval = 10,
            alpha = 1.8, normalize = true, scale = 1.0)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["v2_var"] = m[:v_magnitude_sq]
        sub["zone_edges"] = _zone_cut_edges(sub, cut_edges)
        sub["has_solution"] = false

        busset = _bus_set(sub)
        sub["shared_buses_local"] = [b for b in shared_buses if b in busset]

        @inbounds for I in keys(sub["v2_var"])
            set_start_value(sub["v2_var"][I], 1.0)
        end
    end

    for k in 1:config.max_iter
        println("\n[ADMM-SOCP-hybrid] Iteração $k")

        for (zone_id, sub) in partitions
            m = sub["model"]

            pen = @expression(m, 0.0)
            pen += add_penalty_bus_socp!(
                m,
                sub,
                shared_buses,
                v2_avg,
                duals_v2[zone_id],
                ρV2[],
            )
            pen += add_penalty_branch_pq_socp!(
                m,
                sub,
                sub["zone_edges"],
                zP,
                zQ,
                dualsP[zone_id],
                dualsQ[zone_id],
                ρP[],
                ρQ[],
            )
            pen += add_penalty_branch_rt_socp!(
                m,
                sub,
                sub["zone_edges"],
                zR,
                zT,
                dualsR[zone_id],
                dualsT[zone_id],
                ρR[],
                ρT[],
            )
            @objective(m, Min, sub["obj_base"] + pen)
            optimize!(m)
            sub["has_solution"] = true

            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] =
                compute_real_generation_cost(m, sub["objective_params"])
        end

        v2_prev = isempty(v2_avg) ? v2_avg : copy(v2_avg)
        zP_prev = copy(zP)
        zQ_prev = copy(zQ)
        zR_prev = copy(zR)
        zT_prev = copy(zT)

        compute_v2_consensus!(v2_avg, partitions, shared_buses)
        compute_pq_consensus!(zP, zQ, partitions, cut_edges)
        compute_rt_consensus!(zR, zT, partitions, cut_edges)

        update_duals_bus_socp!(duals_v2, partitions, v2_avg, shared_buses, ρV2[])
        update_duals_branch_pq_socp!(
            dualsP,
            dualsQ,
            partitions,
            zP,
            zQ,
            cut_edges,
            ρP[],
            ρQ[],
        )
        update_duals_branch_rt_socp!(
            dualsR,
            dualsT,
            partitions,
            zR,
            zT,
            cut_edges,
            ρR[],
            ρT[],
        )

        r_bus = primal_residual_bus_socp(partitions, v2_avg, shared_buses)
        r_pq = primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
        r_rt = primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
        r = sqrt(r_bus^2 + r_pq^2 + r_rt^2)

        s_bus = dual_residual_bus_socp(v2_avg, v2_prev, shared_buses)
        s_pq = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
        s_rt = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
        s = sqrt(s_bus^2 + s_pq^2 + s_rt^2)

        push!(r_hist, r)
        push!(s_hist, s)

        ε_pri = begin
            εb = eps_primal_bus_socp(partitions, v2_avg, shared_buses, ϵ_abs, ϵ_rel)
            εpq = eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
            εrt = eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εpq^2 + εrt^2)
        end
        ε_du = begin
            εb = eps_dual_bus_socp(v2_avg, v2_prev, shared_buses, ϵ_abs, ϵ_rel)
            εpq = eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
            εrt = eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εpq^2 + εrt^2)
        end
        push!(eps_r_hist, ε_pri)
        push!(eps_s_hist, ε_du)

        println("‖r‖ = $r | ε_pri = $ε_pri")
        println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoV2_hist, ρV2[])
        push!(rhoP_hist, ρP[])
        push!(rhoQ_hist, ρQ[])
        push!(rhoR_hist, ρR[])
        push!(rhoT_hist, ρT[])

        if r < ε_pri && s < ε_du
            println("Stopping criterion met.")
            break
        end

        adaptive_rho!(ρV2, r_bus, s_bus)
        adaptive_rho!(ρP, r_pq, s_pq)
        adaptive_rho!(ρQ, r_pq, s_pq)
        adaptive_rho!(ρR, r_rt, s_rt)
        adaptive_rho!(ρT, r_rt, s_rt)
    end

    save_socp_hybrid_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
        rhoV2_hist, rhoP_hist, rhoQ_hist, rhoR_hist, rhoT_hist; stem = "admm_socp_hybrid")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-SOCP-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println(
        "ρ(V²)_final = $(ρV2[]) | ρP_final = $(ρP[]) | ρQ_final = $(ρQ[]) | ρR_final = $(ρR[]) | ρT_final = $(ρT[])",
    )

    return total_real_cost, partitions
end
