# ================= SOCP Utils =================
# ---- Otimizador (default SCS) ----
function resolve_socp_optimizer(config)
    if hasproperty(config, :optimizer) && getproperty(config, :optimizer) !== nothing
        return getproperty(config, :optimizer)
    else
        return SCS.Optimizer
    end
end

# setters "seguros" de atributos (não quebram se o solver não reconhecer)
safe_set_optimizer_attr!(m, key, val) = try
    set_optimizer_attribute(m, key, val)
catch
    nothing
end

"Configura tolerâncias e aceleradores razoáveis do SCS."
function set_socp_defaults!(m; eps_abs=1e-6, eps_rel=1e-6, max_iters=200_000,
                            acceleration_lookback=10, acceleration_interval=10,
                            alpha=1.8, normalize=true, scale=1.0)
    safe_set_optimizer_attr!(m, "eps_abs", eps_abs)
    safe_set_optimizer_attr!(m, "eps_rel", eps_rel)
    safe_set_optimizer_attr!(m, "max_iters", max_iters)
    safe_set_optimizer_attr!(m, "acceleration_lookback", acceleration_lookback)
    safe_set_optimizer_attr!(m, "acceleration_interval", acceleration_interval)
    safe_set_optimizer_attr!(m, "alpha", alpha)
    safe_set_optimizer_attr!(m, "normalize", normalize)
    safe_set_optimizer_attr!(m, "scale", scale)
end

# ---- Penalidades ADMM ----
# BUS (V^2)
# BUS (V^2)
function add_penalty_bus_socp!(m, sub, shared_buses::Vector{Int},
                               v2_avg::Dict{Int,Float64}, duals_v2_zone::Dict{Int,Float64},
                               ρV2::Float64)
    bus_ids = sub["bus"]["bus_i"]
    pen = @expression(m, 0.0)   # <-- anônima
    for b in shared_buses
        if b in bus_ids
            V2  = m[:v_magnitude_sq][b, 1]
            V2̄ = get(v2_avg, b, 1.0)
            λ   = get(duals_v2_zone, b, 0.0)
            pen += λ * (V2 - V2̄) + (ρV2 / 2) * (V2 - V2̄)^2
        end
    end
    return pen
end

# BRANCH (R,T)
function add_penalty_branch_rt_socp!(m, sub, zone_edges::Vector{Int},
                                     zR::Dict{Int,Float64}, zT::Dict{Int,Float64},
                                     dualsR_zone::Dict{Int,Float64}, dualsT_zone::Dict{Int,Float64},
                                     ρR::Float64, ρT::Float64)
    pen = @expression(m, 0.0)   # <-- anônima
    for e in zone_edges
        l_local = findfirst(==(e), sub["branch_local_indices"])
        @assert l_local !== nothing "edge $e não está na partição"
        Rloc, Tloc = _rt_from_exprs_socp(m, l_local)
        λR = get(dualsR_zone, e, 0.0); λT = get(dualsT_zone, e, 0.0)
        pen += λR * (Rloc - zR[e]) + (ρR / 2) * (Rloc - zR[e])^2
        pen += λT * (Tloc - zT[e]) + (ρT / 2) * (Tloc - zT[e])^2
    end
    return pen
end

# BRANCH (P,Q)
function add_penalty_branch_pq_socp!(m, sub, zone_edges::Vector{Int},
                                     zP::Dict{Int,Float64}, zQ::Dict{Int,Float64},
                                     dualsP_zone::Dict{Int,Float64}, dualsQ_zone::Dict{Int,Float64},
                                     ρP::Float64, ρQ::Float64)
    pen = @expression(m, 0.0)   # <-- anônima
    for e in zone_edges
        l_local = findfirst(==(e), sub["branch_local_indices"])
        @assert l_local !== nothing "edge $e não está na partição"
        P_loc, Q_loc = _socp_flow_from_exprs_existing(m, l_local)
        λP = get(dualsP_zone, e, 0.0); λQ = get(dualsQ_zone, e, 0.0)
        pen += λP * (P_loc - zP[e]) + (ρP / 2) * (P_loc - zP[e])^2
        pen += λQ * (Q_loc - zQ[e]) + (ρQ / 2) * (Q_loc - zQ[e])^2
    end
    return pen
end

# ---- Consensos ----
compute_v2_consensus!(v2_avg, partitions, shared_buses) = begin
    for b in shared_buses
        acc=0.0; n=0
        for (_, sub) in partitions
            if b in sub["bus"]["bus_i"]
                acc += value(sub["model"][:v_magnitude_sq][b,1]); n+=1
            end
        end
        if n>0; v2_avg[b] = acc/n; end
    end
end

function compute_rt_consensus!(zR::Dict{Int,Float64}, zT::Dict{Int,Float64},
                               partitions, cut_edges::Vector{Int})
    for e in cut_edges
        accR=0.0; accT=0.0; n=0
        for (_, sub) in partitions
            l = findfirst(==(e), sub["branch_local_indices"])
            if l !== nothing
                R,T = _rt_from_exprs_socp(sub["model"], l)
                accR += value(R); accT += value(T); n+=1
            end
        end
        if n>0; zR[e]=accR/n; zT[e]=accT/n; end
    end
end

function compute_pq_consensus!(zP::Dict{Int,Float64}, zQ::Dict{Int,Float64},
                               partitions, cut_edges::Vector{Int})
    for e in cut_edges
        accP=0.0; accQ=0.0; n=0
        for (_, sub) in partitions
            l = findfirst(==(e), sub["branch_local_indices"])
            if l !== nothing
                P,Q = _socp_flow_from_exprs_existing(sub["model"], l)
                accP += value(P); accQ += value(Q); n+=1
            end
        end
        if n>0; zP[e]=accP/n; zQ[e]=accQ/n; end
    end
end

# ---- Duais ----
update_duals_bus_socp!(duals_v2, partitions, v2_avg, shared_buses, ρV2) = begin
    for (z, sub) in partitions
        m = sub["model"]
        for b in shared_buses
            if b in sub["bus"]["bus_i"]
                V2loc = value(m[:v_magnitude_sq][b,1])
                duals_v2[z][b] = get(duals_v2[z], b, 0.0) + ρV2*(V2loc - v2_avg[b])
            end
        end
    end
end

update_duals_branch_rt_socp!(dualsR, dualsT, partitions, zR, zT, cut_edges, ρR, ρT) = begin
    for (z, sub) in partitions
        for e in _zone_cut_edges(sub, cut_edges)
            l = findfirst(==(e), sub["branch_local_indices"])
            R,T = _rt_from_exprs_socp(sub["model"], l)
            dualsR[z][e] = get(dualsR[z], e, 0.0) + ρR*(value(R) - zR[e])
            dualsT[z][e] = get(dualsT[z], e, 0.0) + ρT*(value(T) - zT[e])
        end
    end
end

update_duals_branch_pq_socp!(dualsP, dualsQ, partitions, zP, zQ, cut_edges, ρP, ρQ) = begin
    for (z, sub) in partitions
        for e in _zone_cut_edges(sub, cut_edges)
            l = findfirst(==(e), sub["branch_local_indices"])
            P,Q = _socp_flow_from_exprs_existing(sub["model"], l)
            dualsP[z][e] = get(dualsP[z], e, 0.0) + ρP*(value(P) - zP[e])
            dualsQ[z][e] = get(dualsQ[z], e, 0.0) + ρQ*(value(Q) - zQ[e])
        end
    end
end

# ---- Resíduos e tolerâncias ----
# BUS (V^2)
function primal_residual_bus_socp(partitions, v2_avg, shared_buses)
    acc=0.0
    for (_, sub) in partitions
        for b in shared_buses
            if b in sub["bus"]["bus_i"]
                acc += (value(sub["model"][:v_magnitude_sq][b,1]) - v2_avg[b])^2
            end
        end
    end
    sqrt(acc)
end
dual_residual_bus_socp(v2_avg, v2_prev, shared_buses) =
    sqrt(sum((v2_avg[b] - get(v2_prev, b, 0.0))^2 for b in shared_buses))

function eps_primal_bus_socp(partitions, v2_avg, shared_buses, ϵ_abs, ϵ_rel)
    dim = sum(b in partitions[z]["bus"]["bus_i"] for z in keys(partitions), b in shared_buses)
    z_norm = sqrt(sum(v2_avg[b]^2 for b in shared_buses))
    sqrt(max(1,dim))*ϵ_abs + ϵ_rel*max(z_norm, 1e-12)
end
function eps_dual_bus_socp(v2_avg, v2_prev, shared_buses, ϵ_abs, ϵ_rel)
    dim = length(shared_buses)
    s = dual_residual_bus_socp(v2_avg, v2_prev, shared_buses)
    sqrt(max(1,dim))*ϵ_abs + ϵ_rel*s
end

# BRANCH (R,T)
function primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
    acc=0.0
    for (_, sub) in partitions
        for e in _zone_cut_edges(sub, cut_edges)
            l = findfirst(==(e), sub["branch_local_indices"])
            R,T = _rt_from_exprs_socp(sub["model"], l)
            acc += (value(R)-zR[e])^2 + (value(T)-zT[e])^2
        end
    end
    sqrt(acc)
end
dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges) =
    sqrt(sum((zR[e]-get(zR_prev,e,0.0))^2 + (zT[e]-get(zT_prev,e,0.0))^2 for e in cut_edges))
function eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, ϵ_abs, ϵ_rel)
    dim = 2 * sum(length(_zone_cut_edges(sub, cut_edges)) for (_,sub) in partitions)
    z_norm = sqrt(sum(zR[e]^2 + zT[e]^2 for e in cut_edges))
    sqrt(max(1,dim))*ϵ_abs + ϵ_rel*max(z_norm, 1e-12)
end
function eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, ϵ_abs, ϵ_rel)
    dim = 2 * length(cut_edges)
    s = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
    sqrt(max(1,dim))*ϵ_abs + ϵ_rel*s
end

# BRANCH (P,Q)
function primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
    acc=0.0
    for (_, sub) in partitions
        for e in _zone_cut_edges(sub, cut_edges)
            l = findfirst(==(e), sub["branch_local_indices"])
            P,Q = _socp_flow_from_exprs_existing(sub["model"], l)
            acc += (value(P)-zP[e])^2 + (value(Q)-zQ[e])^2
        end
    end
    sqrt(acc)
end
dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges) =
    sqrt(sum((zP[e]-get(zP_prev,e,0.0))^2 + (zQ[e]-get(zQ_prev,e,0.0))^2 for e in cut_edges))
function eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
    dim = 2 * sum(length(_zone_cut_edges(sub, cut_edges)) for (_,sub) in partitions)
    z_norm = sqrt(sum(zP[e]^2 + zQ[e]^2 for e in cut_edges))
    sqrt(max(1,dim))*ϵ_abs + ϵ_rel*max(z_norm, 1e-12)
end
function eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
    dim = 2 * length(cut_edges)
    s = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
    sqrt(max(1,dim))*ϵ_abs + ϵ_rel*s
end

function adaptive_rho!(ρ::Base.RefValue{Float64}, r::Float64, s::Float64;
                       μ_incr=2.0, μ_decr=2.0, τ_incr=1.05, τ_decr=0.98)
    if r > μ_incr*s
        ρ[] *= τ_incr
    elseif s > μ_decr*r
        ρ[] *= τ_decr
    end
end

# ---- Plots (wrappers) ----
function save_socp_v2rt_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
                              rhoV2_hist, rhoR_hist, rhoT_hist; stem="admm_socp_v2rt")
    plot_admm_convergence(r_hist, s_hist, eps_r_hist, eps_s_hist,
        "plot_convergencia_$(stem).html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoV2_hist), y=rhoV2_hist, type=:scatter, mode="lines+markers", name="ρ(V²)");
        p.layout.title="Evolução de ρ(V²)"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρ(V²)"; p.layout.yaxis.type="log"; p))(),
        "plot_rho_v2_$(stem).html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoR_hist), y=rhoR_hist, type=:scatter, mode="lines+markers", name="ρR");
        p.layout.title="Evolução de ρR"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρR"; p.layout.yaxis.type="log"; p))(),
        "plot_rhoR_$(stem).html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoT_hist), y=rhoT_hist, type=:scatter, mode="lines+markers", name="ρT");
        p.layout.title="Evolução de ρT"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρT"; p.layout.yaxis.type="log"; p))(),
        "plot_rhoT_$(stem).html")
end

function save_socp_branch_pq_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
                                   rhoP_hist, rhoQ_hist; stem="admm_socp_branch")
    plot_admm_convergence(r_hist, s_hist, eps_r_hist, eps_s_hist,
        "plot_convergencia_$(stem).html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoP_hist), y=rhoP_hist, type=:scatter, mode="lines+markers", name="ρP");
        p.layout.title="Evolução de ρP"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρP"; p.layout.yaxis.type="log"; p))(),
        "plot_rhoP_$(stem).html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoQ_hist), y=rhoQ_hist, type=:scatter, mode="lines+markers", name="ρQ");
        p.layout.title="Evolução de ρQ"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρQ"; p.layout.yaxis.type="log"; p))(),
        "plot_rhoQ_$(stem).html")
end

function save_socp_hybrid_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
                                rhoV2_hist, rhoP_hist, rhoQ_hist, rhoR_hist, rhoT_hist;
                                stem="admm_socp_hybrid")
    plot_admm_convergence(r_hist, s_hist, eps_r_hist, eps_s_hist,
        "plot_convergencia_$(stem).html")
    for (lab, hist) in [("ρ(V²)",rhoV2_hist), ("ρP",rhoP_hist), ("ρQ",rhoQ_hist), ("ρR",rhoR_hist), ("ρT",rhoT_hist)]
        fname = replace(lab, "/" => "") |> x->lowercase(replace(x, "ρ"=>"rho"))
        PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(hist), y=hist, type=:scatter, mode="lines+markers", name=lab);
            p.layout.title="Evolução de $(lab)"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title=lab; p.layout.yaxis.type="log"; p))(),
            "plot_$(fname)_$(stem).html")
    end
end

function run_admm_socp(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == SOCP "ADMM-SOCP: base_model deve ser SOCP."

    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    # estado
    ρV2 = Ref(config.rho_v)
    ρR  = Ref(config.rho_flow_p)
    ρT  = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    duals_v2 = Dict{Int,Dict{Int,Float64}}(z=>Dict{Int,Float64}() for z in keys(partitions))
    dualsR   = Dict{Int,Dict{Int,Float64}}(z=>Dict{Int,Float64}() for z in keys(partitions))
    dualsT   = Dict{Int,Dict{Int,Float64}}(z=>Dict{Int,Float64}() for z in keys(partitions))

    v2_avg = Dict{Int,Float64}(); v2_prev = Dict{Int,Float64}()
    zR = Dict{Int,Float64}(e=>0.0 for e in cut_edges); zR_prev = Dict{Int,Float64}()
    zT = Dict{Int,Float64}(e=>0.0 for e in cut_edges); zT_prev = Dict{Int,Float64}()

    # histórico
    r_hist = Float64[]; s_hist = Float64[]
    eps_r_hist = Float64[]; eps_s_hist = Float64[]
    rhoV2_hist = Float64[]; rhoR_hist = Float64[]; rhoT_hist = Float64[]

    optimizer_factory = resolve_socp_optimizer(config)

    # parâmetros de ajuste de ρ
    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98

    for k in 1:config.max_iter
        println("\n[ADMM-SOCP (V²+R,T)] Iteração $k")

        # 1) resolve locais (SCS)
        for (zone_id, sub) in partitions
            m = build_model(sub, SOCP; optimizer = optimizer_factory)
            set_socp_defaults!(m; eps_abs=1e-6, eps_rel=1e-6, max_iters=200_000)

            @expression(m, admm_pen, 0.0)
            admm_pen += add_penalty_bus_socp!(m, sub, shared_buses, v2_avg, duals_v2[zone_id], ρV2[])
            zone_edges = _zone_cut_edges(sub, cut_edges)
            admm_pen += add_penalty_branch_rt_socp!(m, sub, zone_edges, zR, zT, dualsR[zone_id], dualsT[zone_id], ρR[], ρT[])
            @objective(m, Min, m[:generation_cost] + admm_pen)

            optimize!(m)
            sub["model"] = m
            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] = compute_real_generation_cost(m, sub["objective_params"])
        end

        # 2) consensos
        v2_prev = isempty(v2_avg) ? v2_avg : deepcopy(v2_avg)
        zR_prev = deepcopy(zR); zT_prev = deepcopy(zT)
        compute_v2_consensus!(v2_avg, partitions, shared_buses)
        compute_rt_consensus!(zR, zT, partitions, cut_edges)

        # 3) duais
        update_duals_bus_socp!(duals_v2, partitions, v2_avg, shared_buses, ρV2[])
        update_duals_branch_rt_socp!(dualsR, dualsT, partitions, zR, zT, cut_edges, ρR[], ρT[])

        # 4) resíduos + tolerâncias (agregado V² + R,T)
        r_bus = primal_residual_bus_socp(partitions, v2_avg, shared_buses)
        r_rt  = primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
        r     = sqrt(r_bus^2 + r_rt^2)

        s_bus = dual_residual_bus_socp(v2_avg, v2_prev, shared_buses)
        s_rt  = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
        s     = sqrt(s_bus^2 + s_rt^2)

        push!(r_hist, r); push!(s_hist, s)

        ε_pri = begin
            εb = eps_primal_bus_socp(partitions, v2_avg, shared_buses, ϵ_abs, ϵ_rel)
            εe = eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εe^2)
        end
        ε_du  = begin
            εb = eps_dual_bus_socp(v2_avg, v2_prev, shared_buses, ϵ_abs, ϵ_rel)
            εe = eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εe^2)
        end
        push!(eps_r_hist, ε_pri); push!(eps_s_hist, ε_du)

        println("‖r‖ = $r | ε_pri = $ε_pri")
        println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoV2_hist, ρV2[]); push!(rhoR_hist, ρR[]); push!(rhoT_hist, ρT[])

        if r < ε_pri && s < ε_du
            println("Stopping criterion met.")
            break
        end

        # 5) update ρ (mesma regra para os três)
        if r > μ_incr*s
            ρV2[] *= τ_incr; ρR[] *= τ_incr; ρT[] *= τ_incr
        elseif s > μ_decr*r
            ρV2[] *= τ_decr; ρR[] *= τ_decr; ρT[] *= τ_decr
        end
    end

    save_socp_v2rt_plots(r_hist, s_hist, eps_r_hist, eps_s_hist, rhoV2_hist, rhoR_hist, rhoT_hist; stem="admm_socp_no_theta")

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

    ρP = Ref(config.rho_flow_p); ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    dualsP = Dict{Int,Dict{Int,Float64}}(z=>Dict{Int,Float64}() for z in keys(partitions))
    dualsQ = Dict{Int,Dict{Int,Float64}}(z=>Dict{Int,Float64}() for z in keys(partitions))
    zP = Dict{Int,Float64}(e=>0.0 for e in cut_edges); zP_prev = Dict{Int,Float64}()
    zQ = Dict{Int,Float64}(e=>0.0 for e in cut_edges); zQ_prev = Dict{Int,Float64}()

    r_hist=Float64[]; s_hist=Float64[]
    eps_r_hist=Float64[]; eps_s_hist=Float64[]
    rhoP_hist=Float64[]; rhoQ_hist=Float64[]

    optimizer_factory = resolve_socp_optimizer(config)

    # (mantive teus parâmetros mais agressivos)
    μ_incr, μ_decr = 10.0, 10.0
    τ_incr, τ_decr = 2.0, 0.5

    for k in 1:config.max_iter
        println("\n[ADMM-SOCP-branch (P,Q)] Iteração $k")

        # 1) resolve locais
        for (zone_id, sub) in partitions
            m = build_model(sub, SOCP; optimizer = optimizer_factory)
            set_socp_defaults!(m; eps_abs=1e-5, eps_rel=1e-5, max_iters=150_000,
                                  acceleration_lookback=10, acceleration_interval=10,
                                  alpha=1.8, normalize=true, scale=1.0)

            zone_edges = _zone_cut_edges(sub, cut_edges)
            admm_pen = add_penalty_branch_pq_socp!(m, sub, zone_edges, zP, zQ, dualsP[zone_id], dualsQ[zone_id], ρP[], ρQ[])
            @objective(m, Min, m[:generation_cost] + admm_pen)
            optimize!(m)

            sub["model"] = m
            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] = compute_real_generation_cost(m, sub["objective_params"])
        end

        # 2) consensos (P,Q)
        zP_prev = deepcopy(zP); zQ_prev = deepcopy(zQ)
        compute_pq_consensus!(zP, zQ, partitions, cut_edges)

        # 3) duais
        update_duals_branch_pq_socp!(dualsP, dualsQ, partitions, zP, zQ, cut_edges, ρP[], ρQ[])

        # 4) resíduos + tolerâncias
        r = primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
        s = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
        push!(r_hist, r); push!(s_hist, s)

        ε_pri = eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
        ε_du  = eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
        push!(eps_r_hist, ε_pri); push!(eps_s_hist, ε_du)

        println("‖r‖ = $r | ε_pri = $ε_pri")
        println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoP_hist, ρP[]); push!(rhoQ_hist, ρQ[])

        if r < ε_pri && s < ε_du
            println("Stopping criterion met.")
            break
        end

        # 5) ajuste adaptativo
        if r > μ_incr*s
            ρP[] *= τ_incr; ρQ[] *= τ_incr
        elseif s > μ_decr*r
            ρP[] *= τ_decr; ρQ[] *= τ_decr
        end
    end

    save_socp_branch_pq_plots(r_hist, s_hist, eps_r_hist, eps_s_hist, rhoP_hist, rhoQ_hist; stem="admm_socp_branch")

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

    # estado
    ρV2 = Ref(config.rho_v)
    ρP  = Ref(config.rho_flow_p)
    ρQ  = Ref(config.rho_flow_q)
    ρR  = Ref(config.rho_flow_p)
    ρT  = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    duals_v2 = Dict{Int,Dict{Int,Float64}}(z=>Dict{Int,Float64}() for z in keys(partitions))
    dualsP   = Dict{Int,Dict{Int,Float64}}(z=>Dict{Int,Float64}() for z in keys(partitions))
    dualsQ   = Dict{Int,Dict{Int,Float64}}(z=>Dict{Int,Float64}() for z in keys(partitions))
    dualsR   = Dict{Int,Dict{Int,Float64}}(z=>Dict{Int,Float64}() for z in keys(partitions))
    dualsT   = Dict{Int,Dict{Int,Float64}}(z=>Dict{Int,Float64}() for z in keys(partitions))

    v2_avg = Dict{Int,Float64}(); v2_prev = Dict{Int,Float64}()
    zP = Dict{Int,Float64}(e=>0.0 for e in cut_edges); zP_prev = Dict{Int,Float64}()
    zQ = Dict{Int,Float64}(e=>0.0 for e in cut_edges); zQ_prev = Dict{Int,Float64}()
    zR = Dict{Int,Float64}(e=>0.0 for e in cut_edges); zR_prev = Dict{Int,Float64}()
    zT = Dict{Int,Float64}(e=>0.0 for e in cut_edges); zT_prev = Dict{Int,Float64}()

    # histórico
    r_hist=Float64[]; s_hist=Float64[]
    eps_r_hist=Float64[]; eps_s_hist=Float64[]
    rhoV2_hist=Float64[]; rhoP_hist=Float64[]; rhoQ_hist=Float64[]; rhoR_hist=Float64[]; rhoT_hist=Float64[]

    optimizer_factory = resolve_socp_optimizer(config)

    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98

    for k in 1:config.max_iter
        println("\n[ADMM-SOCP-hybrid] Iteração $k")

        # 1) resolve locais
        for (zone_id, sub) in partitions
            m = build_model(sub, SOCP; optimizer = optimizer_factory)
            set_socp_defaults!(m; eps_abs=1e-5, eps_rel=1e-5, max_iters=150_000,
                                  acceleration_lookback=10, acceleration_interval=10,
                                  alpha=1.8, normalize=true, scale=1.0)

            @expression(m, admm_pen, 0.0)
            # BUS (V^2)
            admm_pen += add_penalty_bus_socp!(m, sub, shared_buses, v2_avg, duals_v2[zone_id], ρV2[])
            # BRANCH (P,Q,R,T)
            zone_edges = _zone_cut_edges(sub, cut_edges)
            admm_pen += add_penalty_branch_pq_socp!(m, sub, zone_edges, zP, zQ, dualsP[zone_id], dualsQ[zone_id], ρP[], ρQ[])
            admm_pen += add_penalty_branch_rt_socp!(m, sub, zone_edges, zR, zT, dualsR[zone_id], dualsT[zone_id], ρR[], ρT[])

            @objective(m, Min, m[:generation_cost] + admm_pen)
            optimize!(m)

            sub["model"] = m
            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] = compute_real_generation_cost(m, sub["objective_params"])
        end

        # 2) consensos
        v2_prev = isempty(v2_avg) ? v2_avg : deepcopy(v2_avg)
        zP_prev = deepcopy(zP); zQ_prev = deepcopy(zQ)
        zR_prev = deepcopy(zR); zT_prev = deepcopy(zT)

        compute_v2_consensus!(v2_avg, partitions, shared_buses)
        compute_pq_consensus!(zP, zQ, partitions, cut_edges)
        compute_rt_consensus!(zR, zT, partitions, cut_edges)

        # 3) duais
        update_duals_bus_socp!(duals_v2, partitions, v2_avg, shared_buses, ρV2[])
        update_duals_branch_pq_socp!(dualsP, dualsQ, partitions, zP, zQ, cut_edges, ρP[], ρQ[])
        update_duals_branch_rt_socp!(dualsR, dualsT, partitions, zR, zT, cut_edges, ρR[], ρT[])

        # 4) resíduos + tolerâncias (agregados)
        r_bus = primal_residual_bus_socp(partitions, v2_avg, shared_buses)
        r_pq  = primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
        r_rt  = primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
        r = sqrt(r_bus^2 + r_pq^2 + r_rt^2)

        s_bus = dual_residual_bus_socp(v2_avg, v2_prev, shared_buses)
        s_pq  = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
        s_rt  = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
        s = sqrt(s_bus^2 + s_pq^2 + s_rt^2)

        push!(r_hist, r); push!(s_hist, s)

        ε_pri = begin
            εb = eps_primal_bus_socp(partitions, v2_avg, shared_buses, ϵ_abs, ϵ_rel)
            εpq = eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
            εrt = eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εpq^2 + εrt^2)
        end
        ε_du  = begin
            εb = eps_dual_bus_socp(v2_avg, v2_prev, shared_buses, ϵ_abs, ϵ_rel)
            εpq = eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
            εrt = eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εpq^2 + εrt^2)
        end
        push!(eps_r_hist, ε_pri); push!(eps_s_hist, ε_du)

        println("‖r‖ = $r | ε_pri = $ε_pri")
        println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoV2_hist, ρV2[]); push!(rhoP_hist, ρP[]); push!(rhoQ_hist, ρQ[]); push!(rhoR_hist, ρR[]); push!(rhoT_hist, ρT[])

        if r < ε_pri && s < ε_du
            println("Stopping criterion met.")
            break
        end

        # 5) ρ adaptativo por bloco
        adaptive_rho!(ρV2, r_bus, s_bus)
        adaptive_rho!(ρP,  r_pq,  s_pq)
        adaptive_rho!(ρQ,  r_pq,  s_pq)
        adaptive_rho!(ρR,  r_rt,  s_rt)
        adaptive_rho!(ρT,  r_rt,  s_rt)
    end

    save_socp_hybrid_plots(r_hist, s_hist, eps_r_hist, eps_s_hist, rhoV2_hist, rhoP_hist, rhoQ_hist, rhoR_hist, rhoT_hist; stem="admm_socp_hybrid")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-SOCP-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρ(V²)_final = $(ρV2[]) | ρP_final = $(ρP[]) | ρQ_final = $(ρQ[]) | ρR_final = $(ρR[]) | ρT_final = $(ρT[])")

    return total_real_cost, partitions
end

