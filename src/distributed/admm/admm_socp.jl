# ---- ADMM-SOCP: BUS (V²) + BRANCH (R,T) com estrutura DC-like ----
function run_admm_socp(config::DistributedOPFModel, sys::ADMMSystem)
    @assert config.base_model == SOCP "ADMM-SOCP: base_model deve ser SOCP."

    t0_ns = time_ns()
    iter_done = 0

    partitions = build_socp_partitions(sys)        # Dict{Int,PartitionSOCP}
    shared_buses = sys.shared_buses
    cut_edges = sys.cut_edges

    # estados de ρ
    ρV2 = Ref(config.rho_v)
    ρR = Ref(config.rho_flow_p)
    ρT = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)

    # construtor de modelos
    optimizer_factory = resolve_socp_optimizer(config)

    # mapeamento denso (BUS)
    shared_buses_sorted, b2pos =  prepare_dense_maps(shared_buses)

    # build por zona (padrão DC/AC: usa SystemData da zona)
    for (z, sub) in partitions
        sys_zone = sys.parts[z]

        m = build_model(sys_zone, SOCP; optimizer = optimizer_factory)

        sub.model = m
        sub.obj_base = m[:generation_cost]
        sub.v2_var = m[:v_magnitude_sq]
        sub.objective_params = sys_zone.objective
        _set_model_silent!(m)

        # filtra arestas da zona
        sub.zone_edges = _zone_cut_edges(sub, cut_edges)

        # shared locais + buffers
        prepare_socp_hotpaths!(sub, shared_buses_sorted)

        # cache quadráticos
        _build_quad_bus_expr_socp!(sub, ρV2[])
        _build_quad_branch_rt_expr!(sub, ρR[], ρT[], sub.zone_edges)

        # objetivo "estático" = base + quads (lin será somado por iteração)
        sub.obj_static =
            (sub.obj_base) + sub.quad_bus_expr_socp + sub.quad_rt_expr_socp
        @objective(m, Min, sub.obj_static)
        optimize!(m)
    end

    # vetores densos (BUS)
    v2bar_vec = ones(Float64, length(shared_buses_sorted))
    v2bar_prev_vec = similar(v2bar_vec)

    # duais densos por zona (BUS)
    duals_v2_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )

    # consensos/duais (edges RT)
    zR = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zT = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zR_prev = Dict{Int, Float64}()
    zT_prev = Dict{Int, Float64}()
    dualsR =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsT =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    # históricos
    r_hist = Float64[]
    s_hist = Float64[]
    eps_r_hist = Float64[]
    eps_s_hist = Float64[]
    rhoV2_hist = Float64[]
    rhoR_hist = Float64[]
    rhoT_hist = Float64[]

    zone_ids = collect(keys(partitions))
    prev_ρV2, prev_ρR, prev_ρT = ρV2[], ρR[], ρT[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-SOCP (V²+R,T)] Iteração $k")

        # se ρ mudou, só reconstroi os quadráticos
        if ρV2[] != prev_ρV2
            for (_, sub) in partitions
                _build_quad_bus_expr_socp!(sub, ρV2[])
                sub.obj_static =
                    (sub.obj_base) + sub.quad_bus_expr_socp + sub.quad_rt_expr_socp
            end
            prev_ρV2 = ρV2[]
        end
        if ρR[] != prev_ρR || ρT[] != prev_ρT
            for (_, sub) in partitions
                _build_quad_branch_rt_expr!(sub, ρR[], ρT[], sub.zone_edges)
                sub.obj_static =
                    (sub.obj_base) + sub.quad_bus_expr_socp + sub.quad_rt_expr_socp
            end
            prev_ρR, prev_ρT = ρR[], ρT[]
        end

        α_used = nothing
        v2_target = v2bar_vec

        # resolve local (paralelizável)
        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        @threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            m = sub.model::JuMP.Model

            λz = (λ_zone_buf[z] .= duals_v2_dense[z])

            lin_bus = _build_linear_bus_expr_socp(
                sub,
                shared_buses_sorted,
                b2pos,
                v2_target,
                λz,
                ρV2[],
            )
            lin_rt = _build_linear_branch_rt_expr(
                sub,
                sub.zone_edges,
                zR,
                zT,
                dualsR[z],
                dualsT[z],
                ρR[],
                ρT[],
            )
            @objective(m, Min, sub.obj_static + (lin_bus + lin_rt))
            optimize!(m)

            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, sub.objective_params)
            results[i] = (gen_cost_expr = gen_cost_expr, real_cost = real_cost)
        end

        for (i, z) in enumerate(zone_ids)
            sub = partitions[z]
            out = results[i]
            sub.generation_cost_model_expr = out.gen_cost_expr
            sub.real_generation_cost = out.real_cost
        end

        # ===== consenso + duais =====
        v2bar_prev_vec .= v2bar_vec
        zR_prev = copy(zR)
        zT_prev = copy(zT)

        compute_v2_consensus_vec!(v2bar_vec, partitions, shared_buses_sorted, b2pos)
        compute_rt_consensus!(zR, zT, partitions, cut_edges)

        for z in zone_ids
            update_duals_v2_vec!(duals_v2_dense[z], partitions[z], v2bar_vec,
                v2bar_prev_vec,
                ρV2[], α_used, shared_buses_sorted, b2pos)
        end
        for (z, sub) in partitions
            m = sub.model::JuMP.Model
            @inbounds for e in cut_edges
                l = get(sub.edge_to_local, e, nothing)
                l === nothing && continue
                R, T = _rt_from_exprs_socp(m, l)
                dualsR[z][e] = get(dualsR[z], e, 0.0) + ρR[] * (value(R) - zR[e])
                dualsT[z][e] = get(dualsT[z], e, 0.0) + ρT[] * (value(T) - zT[e])
            end
        end

        # resíduos / tolerâncias
        r_bus =
            primal_residual_bus_socp_vec(partitions, v2bar_vec, shared_buses_sorted, b2pos)
        r_rt = primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
        r = hypot(r_bus, r_rt)

        s_bus = dual_residual_bus_socp_vec(v2bar_vec, v2bar_prev_vec)
        s_rt = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
        s = hypot(s_bus, s_rt)

        push!(r_hist, r)
        push!(s_hist, s)

        ε_pri = begin
            εb = eps_primal_bus_socp_vec(
                partitions,
                v2bar_vec,
                shared_buses_sorted,
                b2pos,
                ϵ_abs,
                ϵ_rel,
            )
            εrt = eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, ϵ_abs, ϵ_rel)
            hypot(εb, εrt)
        end
        ε_du = begin
            εb = eps_dual_bus_socp_vec(v2bar_vec, v2bar_prev_vec, ϵ_abs, ϵ_rel)
            εrt = eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, ϵ_abs, ϵ_rel)
            hypot(εb, εrt)
        end
        push!(eps_r_hist, ε_pri)
        push!(eps_s_hist, ε_du)

        _FAST_USE_LOGS && println("‖r‖ = $r | ε_pri = $ε_pri")
        _FAST_USE_LOGS && println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoV2_hist, ρV2[])
        push!(rhoR_hist, ρR[])
        push!(rhoT_hist, ρT[])

        iter_done = k
        if r < ε_pri && s < ε_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        # atualização adaptativa de ρ
        if adaptive
            # --- BUS (V²) ---
            adaptive_rho!(ρV2, r_bus, s_bus;
                μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)

            # --- BRANCH (R, T) ---
            adaptive_rho!(ρR, r_rt / sqrt(2), s_rt / sqrt(2);
                μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
            adaptive_rho!(ρT, r_rt / sqrt(2), s_rt / sqrt(2);
                μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
        end
    end

    save_socp_v2rt_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
        rhoV2_hist, rhoR_hist, rhoT_hist; stem = "admm_socp_no_theta")

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    println("\n[ADMM-SOCP (V²+R,T)] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρ(V²)_final = $(ρV2[]) | ρR_final = $(ρR[]) | ρT_final = $(ρT[])")
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (SOCP V²+RT): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end

# ---- ADMM-SOCP-branch (P,Q) ----
function run_admm_socp_branch(config::DistributedOPFModel, sys::ADMMSystem)
    @assert config.base_model == SOCP "ADMM-SOCP-branch: base_model deve ser SOCP."

    t0_ns = time_ns()
    iter_done = 0

    partitions = build_socp_partitions(sys)  # Dict{Int,PartitionSOCP}
    cut_edges = sys.cut_edges

    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_socp_optimizer(config)

    # build por zona (padrão DC/AC: usa SystemData da zona)
    for (z, sub) in partitions
        sys_zone = sys.parts[z]

        m = build_model(sys_zone, SOCP; optimizer = optimizer_factory)
        sub.model = m
        sub.obj_base = m[:generation_cost]
        sub.objective_params = sys_zone.objective
        _set_model_silent!(m)
        sub.zone_edges = _zone_cut_edges(sub, cut_edges)

        # 🔧 ridge (apenas se Gurobi) — antes de definir obj_static
        if is_gurobi(m)
            _soft_anchor_branch_expr!(
                sub;
                κV2 = 1e-6,
                κPQ = 1e-8,
                use_bus = true,
                use_pq = true,
            )
        else
            sub.soft_anchor_expr = zero(JuMP.QuadExpr)
        end

        # quadrático do ADMM (P,Q)
        _build_quad_branch_pq_expr!(sub, ρP[], ρQ[], sub.zone_edges)

        # 🔧 obj_static inclui SEMPRE o ridge (mesmo se zero)
        sub.obj_static =
            (sub.obj_base) + sub.quad_pq_expr_socp + sub.soft_anchor_expr
        @objective(m, Min, sub.obj_static)
        optimize!(m)
    end

    # consensos/duais
    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ_prev = Dict{Int, Float64}()
    dualsP =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    r_hist = Float64[]
    s_hist = Float64[]
    eps_r_hist = Float64[]
    eps_s_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]

    prev_ρP, prev_ρQ = ρP[], ρQ[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-SOCP-branch (P,Q)] Iteração $k")

        if ρP[] != prev_ρP || ρQ[] != prev_ρQ
            for (_, sub) in partitions
                _build_quad_branch_pq_expr!(sub, ρP[], ρQ[], sub.zone_edges)
                sub.obj_static =
                    (sub.obj_base) + sub.quad_pq_expr_socp + sub.soft_anchor_expr
            end
            prev_ρP, prev_ρQ = ρP[], ρQ[]
        end

        # locais (correto por zona)
        for (z, sub) in partitions
            m = sub.model::JuMP.Model
            lin = _build_linear_branch_pq_expr(sub, sub.zone_edges, zP, zQ,
                dualsP[z], dualsQ[z], ρP[], ρQ[])
            @objective(m, Min, sub.obj_static + lin)
            optimize!(m)

            sub.generation_cost_model_expr = value(m[:generation_cost])
            sub.real_generation_cost =
                compute_real_generation_cost(m, sub.objective_params)
        end

        zP_prev = copy(zP)
        zQ_prev = copy(zQ)
        compute_pq_consensus!(zP, zQ, partitions, cut_edges)

        for (z, sub) in partitions
            m = sub.model::JuMP.Model
            @inbounds for e in cut_edges
                l = get(sub.edge_to_local, e, nothing)
                l === nothing && continue
                P, Q = _socp_flow_from_exprs_existing(m, l)
                dualsP[z][e] = get(dualsP[z], e, 0.0) + ρP[] * (value(P) - zP[e])
                dualsQ[z][e] = get(dualsQ[z], e, 0.0) + ρQ[] * (value(Q) - zQ[e])
            end
        end

        r = primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
        s = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
        push!(r_hist, r)
        push!(s_hist, s)

        ε_pri = eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
        ε_du = eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
        push!(eps_r_hist, ε_pri)
        push!(eps_s_hist, ε_du)

        _FAST_USE_LOGS && println("‖r‖ = $r | ε_pri = $ε_pri")
        _FAST_USE_LOGS && println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoP_hist, ρP[])
        push!(rhoQ_hist, ρQ[])

        iter_done = k
        if r < ε_pri && s < ε_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        # === Atualização totalmente independente dos ρ (P e Q) ===
        if adaptive
            adaptive_rho!(ρP, r / sqrt(2), s / sqrt(2);
                μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
            adaptive_rho!(ρQ, r / sqrt(2), s / sqrt(2);
                μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
        end
    end

    save_socp_branch_pq_plots(r_hist, s_hist, eps_r_hist, eps_s_hist, rhoP_hist, rhoQ_hist;
        stem = "admm_socp_branch")

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))
    println("\n[ADMM-SOCP-branch (P,Q)] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρP_final = $(ρP[]) | ρQ_final = $(ρQ[])")
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (SOCP P/Q): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end

# ---- ADMM-SOCP-hybrid (V² + P/Q + R/T) ----
function run_admm_socp_hybrid(config::DistributedOPFModel, sys::ADMMSystem)
    @assert config.base_model == SOCP "ADMM-SOCP-hybrid: base_model deve ser SOCP."

    t0_ns = time_ns()
    iter_done = 0

    partitions = build_socp_partitions(sys)  # Dict{Int,PartitionSOCP}
    shared_buses = sys.shared_buses
    cut_edges = sys.cut_edges

    ρV2 = Ref(config.rho_v)
    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ρR = Ref(config.rho_flow_p)
    ρT = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_socp_optimizer(config)

    shared_buses_sorted, b2pos =  prepare_dense_maps(shared_buses)

    for (z, sub) in partitions
        sys_zone = sys.parts[z]

        m = build_model(sys_zone, SOCP; optimizer = optimizer_factory)
        sub.model = m
        sub.obj_base = m[:generation_cost]
        sub.v2_var = m[:v_magnitude_sq]
        sub.objective_params = sys_zone.objective
        _set_model_silent!(m)

        sub.zone_edges = _zone_cut_edges(sub, cut_edges)
        prepare_socp_hotpaths!(sub, shared_buses_sorted)

        _build_quad_bus_expr_socp!(sub, ρV2[])
        _build_quad_branch_pq_expr!(sub, ρP[], ρQ[], sub.zone_edges)
        _build_quad_branch_rt_expr!(sub, ρR[], ρT[], sub.zone_edges)
        sub.obj_static =
            (sub.obj_base) + sub.quad_bus_expr_socp +
            sub.quad_pq_expr_socp + sub.quad_rt_expr_socp
        @objective(m, Min, sub.obj_static)
        optimize!(m)
    end

    # BUS (denso) + EDGES (dict)
    v2bar_vec = ones(Float64, length(shared_buses_sorted))
    v2bar_prev_vec = similar(v2bar_vec)
    duals_v2_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )

    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ_prev = Dict{Int, Float64}()
    zR = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zR_prev = Dict{Int, Float64}()
    zT = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zT_prev = Dict{Int, Float64}()

    dualsP =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsR =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsT =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    r_hist = Float64[]
    s_hist = Float64[]
    eps_r_hist = Float64[]
    eps_s_hist = Float64[]
    rhoV2_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]
    rhoR_hist = Float64[]
    rhoT_hist = Float64[]

    prev_ρV2, prev_ρP, prev_ρQ, prev_ρR, prev_ρT = ρV2[], ρP[], ρQ[], ρR[], ρT[]
    zone_ids = collect(keys(partitions))

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-SOCP-hybrid] Iteração $k")

        # rebuild de quads caso ρ mude
        if ρV2[] != prev_ρV2 || ρP[] != prev_ρP || ρQ[] != prev_ρQ || ρR[] != prev_ρR ||
           ρT[] != prev_ρT
            for (_, sub) in partitions
                (ρV2[] != prev_ρV2) && _build_quad_bus_expr_socp!(sub, ρV2[])
                (ρP[] != prev_ρP || ρQ[] != prev_ρQ) &&
                    _build_quad_branch_pq_expr!(sub, ρP[], ρQ[], sub.zone_edges)
                (ρR[] != prev_ρR || ρT[] != prev_ρT) &&
                    _build_quad_branch_rt_expr!(sub, ρR[], ρT[], sub.zone_edges)
                sub.obj_static =
                    (sub.obj_base) + sub.quad_bus_expr_socp +
                    sub.quad_pq_expr_socp + sub.quad_rt_expr_socp
            end
            prev_ρV2, prev_ρP, prev_ρQ, prev_ρR, prev_ρT = ρV2[], ρP[], ρQ[], ρR[], ρT[]
        end

        α_used = nothing
        v2_target = v2bar_vec

        nz = length(zone_ids)
        @threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            m = sub.model::JuMP.Model
            λz = (λ_zone_buf[z] .= duals_v2_dense[z])

            lin_bus = _build_linear_bus_expr_socp(
                sub,
                shared_buses_sorted,
                b2pos,
                v2_target,
                λz,
                ρV2[],
            )
            lin_pq = _build_linear_branch_pq_expr(
                sub,
                sub.zone_edges,
                zP,
                zQ,
                dualsP[z],
                dualsQ[z],
                ρP[],
                ρQ[],
            )
            lin_rt = _build_linear_branch_rt_expr(
                sub,
                sub.zone_edges,
                zR,
                zT,
                dualsR[z],
                dualsT[z],
                ρR[],
                ρT[],
            )
            @objective(m, Min, sub.obj_static + (lin_bus + lin_pq + lin_rt))
            optimize!(m)

            sub.generation_cost_model_expr = value(m[:generation_cost])
            sub.real_generation_cost =
                compute_real_generation_cost(m, sub.objective_params)
        end

        # consenso + duais reais
        v2bar_prev_vec .= v2bar_vec
        zP_prev = copy(zP)
        zQ_prev = copy(zQ)
        zR_prev = copy(zR)
        zT_prev = copy(zT)

        compute_v2_consensus_vec!(v2bar_vec, partitions, shared_buses_sorted, b2pos)
        compute_pq_consensus!(zP, zQ, partitions, cut_edges)
        compute_rt_consensus!(zR, zT, partitions, cut_edges)

        for z in zone_ids
            update_duals_v2_vec!(duals_v2_dense[z], partitions[z], v2bar_vec,
                v2bar_prev_vec,
                ρV2[], α_used, shared_buses_sorted, b2pos)
        end
        for (z, sub) in partitions
            m = sub.model::JuMP.Model
            @inbounds for e in cut_edges
                l = get(sub.edge_to_local, e, nothing)
                l === nothing && continue
                P, Q = _socp_flow_from_exprs_existing(m, l)
                R, T = _rt_from_exprs_socp(m, l)
                dualsP[z][e] = get(dualsP[z], e, 0.0) + ρP[] * (value(P) - zP[e])
                dualsQ[z][e] = get(dualsQ[z], e, 0.0) + ρQ[] * (value(Q) - zQ[e])
                dualsR[z][e] = get(dualsR[z], e, 0.0) + ρR[] * (value(R) - zR[e])
                dualsT[z][e] = get(dualsT[z], e, 0.0) + ρT[] * (value(T) - zT[e])
            end
        end

        # resíduos / tolerâncias (agregados)
        r_bus =
            primal_residual_bus_socp_vec(partitions, v2bar_vec, shared_buses_sorted, b2pos)
        r_pq = primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
        r_rt = primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
        r = sqrt(r_bus^2 + r_pq^2 + r_rt^2)

        s_bus = dual_residual_bus_socp_vec(v2bar_vec, v2bar_prev_vec)
        s_pq = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
        s_rt = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
        s = sqrt(s_bus^2 + s_pq^2 + s_rt^2)

        push!(r_hist, r)
        push!(s_hist, s)

        ε_pri = begin
            εb = eps_primal_bus_socp_vec(
                partitions,
                v2bar_vec,
                shared_buses_sorted,
                b2pos,
                ϵ_abs,
                ϵ_rel,
            )
            εpq = eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
            εrt = eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εpq^2 + εrt^2)
        end
        ε_du = begin
            εb = eps_dual_bus_socp_vec(v2bar_vec, v2bar_prev_vec, ϵ_abs, ϵ_rel)
            εpq = eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
            εrt = eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εpq^2 + εrt^2)
        end
        push!(eps_r_hist, ε_pri)
        push!(eps_s_hist, ε_du)

        _FAST_USE_LOGS && println("‖r‖ = $r | ε_pri = $ε_pri")
        _FAST_USE_LOGS && println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoV2_hist, ρV2[])
        push!(rhoP_hist, ρP[])
        push!(rhoQ_hist, ρQ[])
        push!(rhoR_hist, ρR[])
        push!(rhoT_hist, ρT[])

        iter_done = k
        if r < ε_pri && s < ε_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        # === Atualização totalmente independente dos ρ (todos os blocos) ===
        if adaptive
            # BUS
            adaptive_rho!(ρV2, r_bus, s_bus;
                μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)

            # BRANCH (P, Q)
            adaptive_rho!(ρP, r_pq / sqrt(2), s_pq / sqrt(2);
                μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
            adaptive_rho!(ρQ, r_pq / sqrt(2), s_pq / sqrt(2);
                μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)

            # BRANCH (R, T)
            adaptive_rho!(ρR, r_rt / sqrt(2), s_rt / sqrt(2);
                μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
            adaptive_rho!(ρT, r_rt / sqrt(2), s_rt / sqrt(2);
                μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
        end
    end

    save_socp_hybrid_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
        rhoV2_hist, rhoP_hist, rhoQ_hist, rhoR_hist, rhoT_hist; stem = "admm_socp_hybrid")

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))
    println("\n[ADMM-SOCP-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println(
        "ρ(V²)_final = $(ρV2[]) | ρP_final = $(ρP[]) | ρQ_final = $(ρQ[]) | ρR_final = $(ρR[]) | ρT_final = $(ρT[])",
    )
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (SOCP híbrido): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end