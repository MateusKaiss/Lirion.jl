function init_socp_v2rt_partition!(
    sub::PartitionSOCP,
    sys_zone,
    optimizer_factory,
    rhoV2::Float64,
    rhoR::Float64,
    rhoT::Float64,
    cut_edges,
)
    prepare_hotpaths!(sub)
    sub.zone_edges = _zone_cut_edges(sub, cut_edges)

    m = build_model(sys_zone, SOCP; optimizer = optimizer_factory)
    sub.model = m
    sub.obj_base = m[:generation_cost]
    sub.v2_var = m[:v_magnitude_sq]
    sub.objective_params = sys_zone.objective
    sub.has_solution = false
    _set_model_silent!(m)

    for I in keys(sub.v2_var)
        set_start_value(sub.v2_var[I], 1.0)
    end

    _build_quad_bus_expr_socp!(sub, rhoV2)
    _build_quad_branch_rt_expr!(sub, rhoR, rhoT, sub.zone_edges)

    sub.obj_static = (sub.obj_base) + sub.quad_bus_expr_socp + sub.quad_rt_expr_socp
    @objective(m, Min, sub.obj_static)
    optimize!(m)

    return nothing
end

# BRANCH (P,Q)
function init_socp_pq_partition!(
    sub::PartitionSOCP,
    sys_zone,
    optimizer_factory,
    rhoP::Float64,
    rhoQ::Float64,
    cut_edges,
)
    prepare_hotpaths!(sub)
    sub.zone_edges = _zone_cut_edges(sub, cut_edges)

    m = build_model(sys_zone, SOCP; optimizer = optimizer_factory)
    sub.model = m
    sub.obj_base = m[:generation_cost]
    sub.objective_params = sys_zone.objective
    sub.has_solution = false
    _set_model_silent!(m)

    if is_gurobi(m)
        _soft_anchor_branch_expr!(
            sub;
            κV2 = 1e-6,
            κPQ = 1e-8,
            use_bus = true,
            use_pq  = true,
        )
    else
        sub.soft_anchor_expr = zero(JuMP.QuadExpr)
    end

    _build_quad_branch_pq_expr!(sub, rhoP, rhoQ, sub.zone_edges)
    sub.obj_static = (sub.obj_base) + sub.quad_pq_expr_socp + sub.soft_anchor_expr
    @objective(m, Min, sub.obj_static)
    optimize!(m)

    return nothing
end

# HYBRID (V² + P/Q + R/T)
function init_socp_hybrid_partition!(
    sub::PartitionSOCP,
    sys_zone,
    optimizer_factory,
    rhoV2::Float64,
    rhoP::Float64,
    rhoQ::Float64,
    rhoR::Float64,
    rhoT::Float64,
    cut_edges,
)
    prepare_hotpaths!(sub)
    sub.zone_edges = _zone_cut_edges(sub, cut_edges)

    m = build_model(sys_zone, SOCP; optimizer = optimizer_factory)
    sub.model = m
    sub.obj_base = m[:generation_cost]
    sub.v2_var = m[:v_magnitude_sq]
    sub.objective_params = sys_zone.objective
    sub.has_solution = false
    _set_model_silent!(m)

    for I in keys(sub.v2_var)
        set_start_value(sub.v2_var[I], 1.0)
    end

    _build_quad_bus_expr_socp!(sub, rhoV2)
    _build_quad_branch_pq_expr!(sub, rhoP, rhoQ, sub.zone_edges)
    _build_quad_branch_rt_expr!(sub, rhoR, rhoT, sub.zone_edges)

    sub.obj_static =
        (sub.obj_base) + sub.quad_bus_expr_socp +
        sub.quad_pq_expr_socp + sub.quad_rt_expr_socp

    @objective(m, Min, sub.obj_static)
    optimize!(m)

    return nothing
end

function run_admm_socp(config::DistributedOPFModel, sys::ADMMSystem)
    t0_ns = time_ns()

    @assert config.base_model == SOCP "ADMM-SOCP: base_model deve ser SOCP."

    # ---------------------------
    # 2) Inputs e parâmetros
    # ---------------------------
    shared_buses = sys.shared_buses
    cut_edges    = sys.cut_edges
    max_iter     = config.max_iter

    rhoV2 = Ref(config.rho_v)
    rhoR  = Ref(config.rho_flow_p)
    rhoT  = Ref(config.rho_flow_q)

    eps_abs = config.eps_primal
    eps_rel = config.eps_dual
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_socp_optimizer(config)

    # ---------------------------
    # 3) Build estrutural
    # ---------------------------
    partitions = build_socp_partitions(sys)
    zone_ids = collect(keys(partitions))
    n_zones = length(zone_ids)

    # hotpaths + edges antes
    for z in zone_ids
        sub = partitions[z]
        prepare_hotpaths!(sub)
        sub.zone_edges = _zone_cut_edges(sub, cut_edges)
    end

    shared_buses_sorted, shared_buses_idx = sorted_ids_and_pos(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    # ---------------------------
    # 4) Alocações (buffers + históricos)
    # ---------------------------
    # BUS (denso)
    v2bar      = ones(Float64, n_shared_buses)
    v2bar_prev = similar(v2bar)

    duals_v2 = Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)
    lambda_zone_buf = Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)

    # RT (dict)
    zR = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zT = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zR_prev = Dict{Int, Float64}()
    zT_prev = Dict{Int, Float64}()

    dualsR = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    dualsT = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    r_hist = Float64[]
    s_hist = Float64[]
    eps_r_hist = Float64[]
    eps_s_hist = Float64[]
    rhoV2_hist = Float64[]
    rhoR_hist  = Float64[]
    rhoT_hist  = Float64[]

    prev_rhoV2 = rhoV2[]
    prev_rhoR  = rhoR[]
    prev_rhoT  = rhoT[]

    results = Vector{NamedTuple{(:gen_cost_expr, :real_cost),
                                Tuple{Float64, Float64}}}(undef, n_zones)

    # ---------------------------
    # 5) Init por zona
    # ---------------------------
    for z in zone_ids
        sub = partitions[z]
        init_socp_v2rt_partition!(
            sub,
            sys.parts[z],
            optimizer_factory,
            rhoV2[],
            rhoR[],
            rhoT[],
            cut_edges,
        )
    end

    # ---------------------------
    # 6) Loop ADMM
    # ---------------------------
    iter_done = 0
    for k in 1:max_iter
        iter_done = k
        _FAST_USE_LOGS && println("\n[ADMM-SOCP (V²+R,T)] Iteração $k")

        if rhoV2[] != prev_rhoV2
            for sub in values(partitions)
                _build_quad_bus_expr_socp!(sub, rhoV2[])
                sub.obj_static = sub.obj_base + sub.quad_bus_expr_socp + sub.quad_rt_expr_socp
            end
            prev_rhoV2 = rhoV2[]
        end

        if (rhoR[] != prev_rhoR) || (rhoT[] != prev_rhoT)
            for sub in values(partitions)
                _build_quad_branch_rt_expr!(sub, rhoR[], rhoT[], sub.zone_edges)
                sub.obj_static = sub.obj_base + sub.quad_bus_expr_socp + sub.quad_rt_expr_socp
            end
            prev_rhoR = rhoR[]
            prev_rhoT = rhoT[]
        end

        α_used = nothing
        v2_target = v2bar

        Threads.@threads for i in 1:n_zones
            z = zone_ids[i]
            sub = partitions[z]
            m = sub.model::JuMP.Model

            lambdaz = (lambda_zone_buf[z] .= duals_v2[z])

            lin_bus = _build_linear_bus_expr_socp(
                sub, shared_buses_sorted, shared_buses_idx, v2_target, lambdaz, rhoV2[]
            )
            lin_rt = _build_linear_branch_rt_expr(
                sub, sub.zone_edges, zR, zT, dualsR[z], dualsT[z], rhoR[], rhoT[]
            )

            @objective(m, Min, sub.obj_static + (lin_bus + lin_rt))
            optimize!(m)

            results[i] = (
                gen_cost_expr = value(m[:generation_cost]),
                real_cost     = compute_real_generation_cost(m, sub.objective_params),
            )
        end

        for (i, z) in enumerate(zone_ids)
            sub = partitions[z]
            out = results[i]
            sub.has_solution = true
            sub.generation_cost_model_expr = out.gen_cost_expr
            sub.real_generation_cost       = out.real_cost
        end

        v2bar_prev .= v2bar
        zR_prev = copy(zR)
        zT_prev = copy(zT)

        compute_v2_consensus_vec!(v2bar, partitions, shared_buses_sorted, shared_buses_idx)
        compute_rt_consensus!(zR, zT, partitions, cut_edges)

        for z in zone_ids
            update_duals_v2_vec!(
                duals_v2[z],
                partitions[z],
                v2bar,
                v2bar_prev,
                rhoV2[],
                α_used,
                shared_buses_sorted,
                shared_buses_idx,
            )
        end

        # RT duals (por extração)
        for z in zone_ids
            sub = partitions[z]
            m = sub.model::JuMP.Model
            @inbounds for e in sub.zone_edges
                l = get(sub.edge_to_local, e, nothing)
                l === nothing && continue
                R, T = _rt_from_exprs_socp(m, l)
                dualsR[z][e] = get(dualsR[z], e, 0.0) + rhoR[] * (value(R) - zR[e])
                dualsT[z][e] = get(dualsT[z], e, 0.0) + rhoT[] * (value(T) - zT[e])
            end
        end

        r_bus = primal_residual_bus_socp_vec(partitions, v2bar, shared_buses_sorted, shared_buses_idx)
        r_rt  = primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
        r = hypot(r_bus, r_rt)

        s_bus = dual_residual_bus_socp_vec(v2bar, v2bar_prev)
        s_rt  = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
        s = hypot(s_bus, s_rt)

        eps_r = begin
            epsb  = eps_primal_bus_socp_vec(partitions, v2bar, shared_buses_sorted, shared_buses_idx, eps_abs, eps_rel)
            epsrt = eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, eps_abs, eps_rel)
            hypot(epsb, epsrt)
        end
        eps_s = begin
            epsb  = eps_dual_bus_socp_vec(v2bar, v2bar_prev, eps_abs, eps_rel)
            epsrt = eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, eps_abs, eps_rel)
            hypot(epsb, epsrt)
        end

        push!(r_hist, r)
        push!(s_hist, s)
        push!(eps_r_hist, eps_r)
        push!(eps_s_hist, eps_s)

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        push!(rhoV2_hist, rhoV2[])
        push!(rhoR_hist,  rhoR[])
        push!(rhoT_hist,  rhoT[])

        if r < eps_r && s < eps_s
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_dc!(rhoV2, r_bus, s_bus; mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
            adaptive_rho_dc!(rhoR,  r_rt / sqrt(2), s_rt / sqrt(2); mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
            adaptive_rho_dc!(rhoT,  r_rt / sqrt(2), s_rt / sqrt(2); mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
        end
    end

    # ---------------------------
    # 7) Pós-processo
    # ---------------------------
    save_socp_v2rt_plots(
        r_hist, s_hist, eps_r_hist, eps_s_hist,
        rhoV2_hist, rhoR_hist, rhoT_hist;
        stem = "admm_socp_no_theta",
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-SOCP (V²+R,T)] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho(V²) inicial: $(config.rho_v) | rho(V²) final: $(rhoV2[])")
    println("rhoR    inicial: $(config.rho_flow_p) | rhoR final: $(rhoR[])")
    println("rhoT    inicial: $(config.rho_flow_q) | rhoT final: $(rhoT[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (SOCP V²+RT): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end
function run_admm_socp_branch(config::DistributedOPFModel, sys::ADMMSystem)
    t0_ns = time_ns()

    @assert config.base_model == SOCP "ADMM-SOCP-branch: base_model deve ser SOCP."

    # ---------------------------
    # 2) Inputs e parâmetros
    # ---------------------------
    cut_edges = sys.cut_edges
    max_iter  = config.max_iter

    rhoP = Ref(config.rho_flow_p)
    rhoQ = Ref(config.rho_flow_q)

    eps_abs = config.eps_primal
    eps_rel = config.eps_dual
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_socp_optimizer(config)

    # ---------------------------
    # 3) Build estrutural
    # ---------------------------
    partitions = build_socp_partitions(sys)
    zone_ids = collect(keys(partitions))
    n_zones = length(zone_ids)

    # hotpaths + zone_edges antes
    for z in zone_ids
        sub = partitions[z]
        prepare_hotpaths!(sub)
        sub.zone_edges = _zone_cut_edges(sub, cut_edges)
    end

    # ---------------------------
    # 4) Alocações (buffers + históricos)
    # ---------------------------
    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    dualsP = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    dualsQ = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    r_hist = Float64[]
    s_hist = Float64[]
    eps_r_hist = Float64[]
    eps_s_hist = Float64[]
    rhoP_hist  = Float64[]
    rhoQ_hist  = Float64[]

    prev_rhoP = rhoP[]
    prev_rhoQ = rhoQ[]

    # ---------------------------
    # 5) Init por zona
    # ---------------------------
    for z in zone_ids
        sub = partitions[z]
        init_socp_pq_partition!(
            sub,
            sys.parts[z],
            optimizer_factory,
            rhoP[],
            rhoQ[],
            cut_edges,
        )
    end

    # ---------------------------
    # 6) Loop ADMM
    # ---------------------------
    iter_done = 0
    for k in 1:max_iter
        iter_done = k
        _FAST_USE_LOGS && println("\n[ADMM-SOCP-branch (P,Q)] Iteração $k")

        if (rhoP[] != prev_rhoP) || (rhoQ[] != prev_rhoQ)
            for sub in values(partitions)
                _build_quad_branch_pq_expr!(sub, rhoP[], rhoQ[], sub.zone_edges)
                sub.obj_static = sub.obj_base + sub.quad_pq_expr_socp + sub.soft_anchor_expr
            end
            prev_rhoP = rhoP[]
            prev_rhoQ = rhoQ[]
        end

        for z in zone_ids
            sub = partitions[z]
            m = sub.model::JuMP.Model

            lin = _build_linear_branch_pq_expr(
                sub, sub.zone_edges, zP, zQ, dualsP[z], dualsQ[z], rhoP[], rhoQ[]
            )

            @objective(m, Min, sub.obj_static + lin)
            optimize!(m)

            sub.has_solution = true
            sub.generation_cost_model_expr = value(m[:generation_cost])
            sub.real_generation_cost       = compute_real_generation_cost(m, sub.objective_params)
        end

        zP_prev = copy(zP)
        zQ_prev = copy(zQ)

        compute_pq_consensus!(zP, zQ, partitions, cut_edges)

        for z in zone_ids
            sub = partitions[z]
            m = sub.model::JuMP.Model
            @inbounds for e in sub.zone_edges
                l = get(sub.edge_to_local, e, nothing)
                l === nothing && continue
                P, Q = _socp_flow_from_exprs_existing(m, l)
                dualsP[z][e] = get(dualsP[z], e, 0.0) + rhoP[] * (value(P) - zP[e])
                dualsQ[z][e] = get(dualsQ[z], e, 0.0) + rhoQ[] * (value(Q) - zQ[e])
            end
        end

        r = primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
        s = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)

        eps_r = eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, eps_abs, eps_rel)
        eps_s = eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, eps_abs, eps_rel)

        push!(r_hist, r)
        push!(s_hist, s)
        push!(eps_r_hist, eps_r)
        push!(eps_s_hist, eps_s)

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        push!(rhoP_hist, rhoP[])
        push!(rhoQ_hist, rhoQ[])

        if r < eps_r && s < eps_s
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_dc!(rhoP, r / sqrt(2), s / sqrt(2); mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
            adaptive_rho_dc!(rhoQ, r / sqrt(2), s / sqrt(2); mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
        end
    end

    # ---------------------------
    # 7) Pós-processo
    # ---------------------------
    save_socp_branch_pq_plots(
        r_hist, s_hist, eps_r_hist, eps_s_hist, rhoP_hist, rhoQ_hist;
        stem = "admm_socp_branch",
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-SOCP-branch (P,Q)] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rhoP inicial: $(config.rho_flow_p) | rhoP final: $(rhoP[])")
    println("rhoQ inicial: $(config.rho_flow_q) | rhoQ final: $(rhoQ[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (SOCP P/Q): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end


function run_admm_socp_hybrid(config::DistributedOPFModel, sys::ADMMSystem)
    t0_ns = time_ns()

    @assert config.base_model == SOCP "ADMM-SOCP-hybrid: base_model deve ser SOCP."

    # ---------------------------
    # 2) Inputs e parâmetros
    # ---------------------------
    shared_buses = sys.shared_buses
    cut_edges    = sys.cut_edges
    max_iter     = config.max_iter

    rhoV2 = Ref(config.rho_v)
    rhoP  = Ref(config.rho_flow_p)
    rhoQ  = Ref(config.rho_flow_q)

    # reuso (como você já fazia)
    rhoR  = Ref(config.rho_flow_p)
    rhoT  = Ref(config.rho_flow_q)

    eps_abs = config.eps_primal
    eps_rel = config.eps_dual
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_socp_optimizer(config)

    # ---------------------------
    # 3) Build estrutural
    # ---------------------------
    partitions = build_socp_partitions(sys)
    zone_ids = collect(keys(partitions))
    n_zones = length(zone_ids)

    for z in zone_ids
        sub = partitions[z]
        prepare_hotpaths!(sub)
        sub.zone_edges = _zone_cut_edges(sub, cut_edges)
    end

    shared_buses_sorted, shared_buses_idx = sorted_ids_and_pos(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    # ---------------------------
    # 4) Alocações (buffers + históricos)
    # ---------------------------
    v2bar      = ones(Float64, n_shared_buses)
    v2bar_prev = similar(v2bar)

    duals_v2 = Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)
    lambda_zone_buf = Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)

    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zR = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zT = Dict{Int, Float64}(e => 0.0 for e in cut_edges)

    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()
    zR_prev = Dict{Int, Float64}()
    zT_prev = Dict{Int, Float64}()

    dualsP = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    dualsQ = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    dualsR = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    dualsT = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    r_hist = Float64[]
    s_hist = Float64[]
    eps_r_hist = Float64[]
    eps_s_hist = Float64[]

    rhoV2_hist = Float64[]
    rhoP_hist  = Float64[]
    rhoQ_hist  = Float64[]
    rhoR_hist  = Float64[]
    rhoT_hist  = Float64[]

    prev_rhoV2 = rhoV2[]
    prev_rhoP  = rhoP[]
    prev_rhoQ  = rhoQ[]
    prev_rhoR  = rhoR[]
    prev_rhoT  = rhoT[]

    # ---------------------------
    # 5) Init por zona
    # ---------------------------
    for z in zone_ids
        sub = partitions[z]
        init_socp_hybrid_partition!(
            sub,
            sys.parts[z],
            optimizer_factory,
            rhoV2[],
            rhoP[],
            rhoQ[],
            rhoR[],
            rhoT[],
            cut_edges,
        )
    end

    # ---------------------------
    # 6) Loop ADMM
    # ---------------------------
    iter_done = 0
    for k in 1:max_iter
        iter_done = k
        _FAST_USE_LOGS && println("\n[ADMM-SOCP-hybrid] Iteração $k")

        if (rhoV2[] != prev_rhoV2) || (rhoP[] != prev_rhoP) || (rhoQ[] != prev_rhoQ) ||
           (rhoR[] != prev_rhoR)  || (rhoT[] != prev_rhoT)

            for sub in values(partitions)
                (rhoV2[] != prev_rhoV2) && _build_quad_bus_expr_socp!(sub, rhoV2[])
                ((rhoP[] != prev_rhoP) || (rhoQ[] != prev_rhoQ)) && _build_quad_branch_pq_expr!(sub, rhoP[], rhoQ[], sub.zone_edges)
                ((rhoR[] != prev_rhoR) || (rhoT[] != prev_rhoT)) && _build_quad_branch_rt_expr!(sub, rhoR[], rhoT[], sub.zone_edges)

                sub.obj_static = sub.obj_base + sub.quad_bus_expr_socp + sub.quad_pq_expr_socp + sub.quad_rt_expr_socp
            end

            prev_rhoV2 = rhoV2[]
            prev_rhoP  = rhoP[]
            prev_rhoQ  = rhoQ[]
            prev_rhoR  = rhoR[]
            prev_rhoT  = rhoT[]
        end

        α_used = nothing
        v2_target = v2bar

        Threads.@threads for i in 1:n_zones
            z = zone_ids[i]
            sub = partitions[z]
            m = sub.model::JuMP.Model

            lambdaz = (lambda_zone_buf[z] .= duals_v2[z])

            lin_bus = _build_linear_bus_expr_socp(sub, shared_buses_sorted, shared_buses_idx, v2_target, lambdaz, rhoV2[])
            lin_pq  = _build_linear_branch_pq_expr(sub, sub.zone_edges, zP, zQ, dualsP[z], dualsQ[z], rhoP[], rhoQ[])
            lin_rt  = _build_linear_branch_rt_expr(sub, sub.zone_edges, zR, zT, dualsR[z], dualsT[z], rhoR[], rhoT[])

            @objective(m, Min, sub.obj_static + (lin_bus + lin_pq + lin_rt))
            optimize!(m)

            sub.has_solution = true
            sub.generation_cost_model_expr = value(m[:generation_cost])
            sub.real_generation_cost       = compute_real_generation_cost(m, sub.objective_params)
        end

        v2bar_prev .= v2bar
        zP_prev = copy(zP); zQ_prev = copy(zQ); zR_prev = copy(zR); zT_prev = copy(zT)

        compute_v2_consensus_vec!(v2bar, partitions, shared_buses_sorted, shared_buses_idx)
        compute_pq_consensus!(zP, zQ, partitions, cut_edges)
        compute_rt_consensus!(zR, zT, partitions, cut_edges)

        for z in zone_ids
            update_duals_v2_vec!(
                duals_v2[z],
                partitions[z],
                v2bar,
                v2bar_prev,
                rhoV2[],
                α_used,
                shared_buses_sorted,
                shared_buses_idx,
            )
        end

        for z in zone_ids
            sub = partitions[z]
            m = sub.model::JuMP.Model
            @inbounds for e in sub.zone_edges
                l = get(sub.edge_to_local, e, nothing)
                l === nothing && continue

                P, Q = _socp_flow_from_exprs_existing(m, l)
                R, T = _rt_from_exprs_socp(m, l)

                dualsP[z][e] = get(dualsP[z], e, 0.0) + rhoP[] * (value(P) - zP[e])
                dualsQ[z][e] = get(dualsQ[z], e, 0.0) + rhoQ[] * (value(Q) - zQ[e])
                dualsR[z][e] = get(dualsR[z], e, 0.0) + rhoR[] * (value(R) - zR[e])
                dualsT[z][e] = get(dualsT[z], e, 0.0) + rhoT[] * (value(T) - zT[e])
            end
        end

        r_bus = primal_residual_bus_socp_vec(partitions, v2bar, shared_buses_sorted, shared_buses_idx)
        r_pq  = primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
        r_rt  = primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
        r = sqrt(r_bus^2 + r_pq^2 + r_rt^2)

        s_bus = dual_residual_bus_socp_vec(v2bar, v2bar_prev)
        s_pq  = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
        s_rt  = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
        s = sqrt(s_bus^2 + s_pq^2 + s_rt^2)

        eps_r = begin
            epsb  = eps_primal_bus_socp_vec(partitions, v2bar, shared_buses_sorted, shared_buses_idx, eps_abs, eps_rel)
            epspq = eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, eps_abs, eps_rel)
            epsrt = eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, eps_abs, eps_rel)
            sqrt(epsb^2 + epspq^2 + epsrt^2)
        end
        eps_s = begin
            epsb  = eps_dual_bus_socp_vec(v2bar, v2bar_prev, eps_abs, eps_rel)
            epspq = eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, eps_abs, eps_rel)
            epsrt = eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, eps_abs, eps_rel)
            sqrt(epsb^2 + epspq^2 + epsrt^2)
        end

        push!(r_hist, r)
        push!(s_hist, s)
        push!(eps_r_hist, eps_r)
        push!(eps_s_hist, eps_s)

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        push!(rhoV2_hist, rhoV2[])
        push!(rhoP_hist,  rhoP[])
        push!(rhoQ_hist,  rhoQ[])
        push!(rhoR_hist,  rhoR[])
        push!(rhoT_hist,  rhoT[])

        if r < eps_r && s < eps_s
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_dc!(rhoV2, r_bus, s_bus; mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)

            adaptive_rho_dc!(rhoP, r_pq / sqrt(2), s_pq / sqrt(2); mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
            adaptive_rho_dc!(rhoQ, r_pq / sqrt(2), s_pq / sqrt(2); mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)

            adaptive_rho_dc!(rhoR, r_rt / sqrt(2), s_rt / sqrt(2); mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
            adaptive_rho_dc!(rhoT, r_rt / sqrt(2), s_rt / sqrt(2); mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
        end
    end

    # ---------------------------
    # 7) Pós-processo
    # ---------------------------
    save_socp_hybrid_plots(
        r_hist, s_hist, eps_r_hist, eps_s_hist,
        rhoV2_hist, rhoP_hist, rhoQ_hist, rhoR_hist, rhoT_hist;
        stem = "admm_socp_hybrid",
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-SOCP-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho(V²) inicial: $(config.rho_v) | rho(V²) final: $(rhoV2[])")
    println("rhoP    inicial: $(config.rho_flow_p) | rhoP final: $(rhoP[])")
    println("rhoQ    inicial: $(config.rho_flow_q) | rhoQ final: $(rhoQ[])")
    println("rhoR    inicial: $(config.rho_flow_p) | rhoR final: $(rhoR[])")
    println("rhoT    inicial: $(config.rho_flow_q) | rhoT final: $(rhoT[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (SOCP híbrido): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end
