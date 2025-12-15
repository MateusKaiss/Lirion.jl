###########################
# HELPERS GENÉRICOS AC    #
###########################

function resolve_optimizer_ac(config)
    if hasproperty(config, :optimizer) && getproperty(config, :optimizer) !== nothing
        return getproperty(config, :optimizer)
    else
        return optimizer_with_attributes(
            Ipopt.Optimizer,

            # ---------------------------
            # Saída silenciosa
            # ---------------------------
            "print_level" => 0,
            "sb" => "yes",
            "print_timing_statistics" => "no",
            "print_user_options" => "no",
            "file_print_level" => 0,

            # ---------------------------
            # Critérios de parada
            # ---------------------------
            "tol" => 1e-7,
            "acceptable_tol" => 1e-5,
            "acceptable_iter" => 10,
            "dual_inf_tol" => 1e-4,
            "constr_viol_tol" => 1e-4,
            "compl_inf_tol" => 1e-4,
            "acceptable_dual_inf_tol" => 1e-2,
            "acceptable_constr_viol_tol" => 1e-3,
            "acceptable_compl_inf_tol" => 1e-3,

            # ---------------------------
            # NLP & Bounds
            # ---------------------------
            "fixed_variable_treatment" => "make_parameter_nodual",
            "bound_relax_factor" => 1e-6,
            "check_derivatives_for_naninf" => "no",
            "jacobian_approximation" => "exact",
            "gradient_approximation" => "exact",
            "kappa_d" => 1e-5,

            # ---------------------------
            # Escalonamento
            # ---------------------------
            "nlp_scaling_max_gradient" => 100.0,
            "nlp_scaling_min_value" => 1e-8,

            # ---------------------------
            # Inicialização
            # ---------------------------
            "bound_push" => 1e-8,
            "bound_frac" => 1e-6,
            "slack_bound_push" => 1e-8,
            "slack_bound_frac" => 1e-6,

            # ---------------------------
            # Linear Solver & Refinement
            # ---------------------------
            "linear_solver" => "mumps",
            "linear_scaling_on_demand" => "yes",
            "fast_step_computation" => "no",
            "min_refinement_steps" => 1,
            "max_refinement_steps" => 8,
            "neg_curv_test_tol" => 1e-12,
            "neg_curv_test_reg" => "yes",
            "max_hessian_perturbation" => 1e20,
            "min_hessian_perturbation" => 1e-20,
            "mumps_mem_percent" => 500,
            "mumps_pivtol" => 1e-6,
            "mumps_pivtolmax" => 0.1,
            "mumps_pivot_order" => 7,
            "hessian_approximation" => "exact",
            "limited_memory_max_history" => 6,
        )
    end
end

@inline function _ac_flow_from_exprs_existing(
    part::PartitionAC,
    l_local::Int;
    t::Int = 1,
)
    m = part.model::JuMP.Model
    return m[:p_flow_branches][l_local, t], m[:q_flow_branches][l_local, t]
end

@inline function _val_or_start_ac(m::JuMP.Model, x::JuMP.VariableRef)
    if JuMP.has_values(m)
        return JuMP.value(x)
    else
        sv = JuMP.start_value(x)
        return sv === nothing ? 0.0 : sv
    end
end

@inline _theta_val_ac(
    m::JuMP.Model,
    θ::JuMP.Containers.DenseAxisArray,
    b::Int,
    t::Int,
) = _val_or_start_ac(m, θ[b, t])

@inline _vmag_val(
    m::JuMP.Model,
    V::JuMP.Containers.DenseAxisArray,
    b::Int,
    t::Int,
) = _val_or_start_ac(m, V[b, t])


@inline function snapshot_theta_ac!(part::PartitionAC)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.theta_buf
    bus_ids = part.bus.id
    @inbounds for i in eachindex(bus_ids)
        buf[i] = _theta_val_ac(m, θ, bus_ids[i], 1)
    end
    return buf
end

@inline function snapshot_vmag!(part::PartitionAC)
    m = part.model::JuMP.Model
    V = part.vmag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.vmag_buf
    bus_ids = part.bus.id
    @inbounds for i in eachindex(bus_ids)
        buf[i] = _vmag_val(m, V, bus_ids[i], 1)
    end
    return buf
end

@inline function snapshot_theta_ac_shared!(part::PartitionAC)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.theta_shared_buf
    pos = part.shared_pos_local
    bus_ids = part.bus.id
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] = _theta_val_ac(m, θ, bus_ids[bidx], 1)
    end
    return buf
end

@inline function snapshot_vmag_shared!(part::PartitionAC)
    m = part.model::JuMP.Model
    V = part.vmag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.vmag_shared_buf
    pos = part.shared_pos_local
    bus_ids = part.bus.id
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] = _vmag_val(m, V, bus_ids[bidx], 1)
    end
    return buf
end


function _build_quad_bus_expr_ac!(part::PartitionAC, ρθ::Float64, ρV::Float64)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    V = part.vmag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    q = zero(JuMP.QuadExpr)
    for b in part.shared_buses_local
        q += 0.5 * ρθ * θ[b, 1]^2 + 0.5 * ρV * V[b, 1]^2
    end
    part.quad_bus_expr_ac = q
    part.rho_bus_theta_cached = ρθ
    part.rho_bus_v_cached = ρV
    return q
end

function _build_linear_bus_expr_ac!(
    part::PartitionAC,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    λθ_zone::Vector{Float64},
    λV_zone::Vector{Float64},
    ρθ::Float64,
    ρV::Float64,
)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    V = part.vmag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    lin = zero(JuMP.AffExpr)
    for b in part.shared_buses_local
        j = b2pos[b]
        lin += (λθ_zone[j] - ρθ * θ̄_vec[j]) * θ[b, 1]
        lin += (λV_zone[j] - ρV * V̄_vec[j]) * V[b, 1]
    end
    return lin
end

function _build_quad_branch_expr_ac!(
    part::PartitionAC,
    ρP::Float64,
    ρQ::Float64,
    zone_edges::Vector{Int},
)
    m = part.model::JuMP.Model
    q = zero(JuMP.QuadExpr)
    for e in zone_edges
        l_local = get(part.edge_to_local, e, nothing)
        l_local === nothing && continue
        P_loc, Q_loc = _ac_flow_from_exprs_existing(part, l_local)
        q += 0.5 * ρP * P_loc^2 + 0.5 * ρQ * Q_loc^2
    end
    part.quad_branch_expr_ac = q
    part.rho_branch_p_cached = ρP
    part.rho_branch_q_cached = ρQ
    return q
end

function _build_linear_branch_expr_ac!(
    part::PartitionAC,
    zone_edges::Vector{Int},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    dualsP_zone::Dict{Int, Float64},
    dualsQ_zone::Dict{Int, Float64},
    ρP::Float64,
    ρQ::Float64,
)
    m = part.model::JuMP.Model
    lin = zero(JuMP.AffExpr)
    for e in zone_edges
        l_local = get(part.edge_to_local, e, nothing)
        l_local === nothing && continue
        P_loc, Q_loc = _ac_flow_from_exprs_existing(part, l_local)
        λP = get(dualsP_zone, e, 0.0)
        λQ = get(dualsQ_zone, e, 0.0)
        zPe = get(zP, e, 0.0)
        zQe = get(zQ, e, 0.0)
        lin += λP * (P_loc - zPe) + λQ * (Q_loc - zQe) -
               ρP * zPe * P_loc - ρQ * zQe * Q_loc
    end
    return lin
end

function compute_ac_bus_consensus_vec!(
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    partitions::Dict{Int, PartitionAC},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    fill!(θ̄_vec, 0.0)
    fill!(V̄_vec, 0.0)
    sums_sin = zeros(Float64, length(shared_buses_sorted))
    sums_cos = zeros(Float64, length(shared_buses_sorted))
    counts = zeros(Int, length(shared_buses_sorted))

    for (_, part) in partitions
        θs = snapshot_theta_ac_shared!(part)
        Vs = snapshot_vmag_shared!(part)
        lbs = part.shared_buses_local
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            s, c = sincos(θs[i])
            sums_sin[j] += s
            sums_cos[j] += c
            V̄_vec[j] += Vs[i]
            counts[j] += 1
        end
    end

    @inbounds for j in eachindex(θ̄_vec)
        if counts[j] != 0
            θ̄_vec[j] = atan(sums_sin[j], sums_cos[j])
            V̄_vec[j] /= counts[j]
        end
    end
end

function compute_ac_branch_consensus!(
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    partitions::Dict{Int, PartitionAC},
    cut_edges::Vector{Int},
)
    sumsP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    sumsQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    cntP = Dict{Int, Int}(e => 0 for e in cut_edges)
    cntQ = Dict{Int, Int}(e => 0 for e in cut_edges)

    for (_, part) in partitions
        m = part.model::JuMP.Model
        for e in part.zone_edges
            l_local = get(part.edge_to_local, e, nothing)
            l_local === nothing && continue
            P_loc, Q_loc = _ac_flow_from_exprs_existing(part, l_local)
            pval = JuMP.value(P_loc)
            qval = JuMP.value(Q_loc)
            if isfinite(pval)
                sumsP[e] += pval
                cntP[e] += 1
            end
            if isfinite(qval)
                sumsQ[e] += qval
                cntQ[e] += 1
            end
        end
    end

    for e in cut_edges
        if cntP[e] > 0
            zP[e] = sumsP[e] / cntP[e]
        end
        if cntQ[e] > 0
            zQ[e] = sumsQ[e] / cntQ[e]
        end
    end
    return nothing
end

@inline function update_duals_bus_ac_vec!(
    dualsθ_vec_zone::Vector{Float64},
    dualsV_vec_zone::Vector{Float64},
    part::PartitionAC,
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    ρθ::Float64,
    ρV::Float64,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int};
    α::Union{Nothing, Float64} = nothing,
    θ̄_prev::Union{Nothing, Vector{Float64}} = nothing,
    V̄_prev::Union{Nothing, Vector{Float64}} = nothing,
)
    θs = snapshot_theta_ac_shared!(part)
    Vs = snapshot_vmag_shared!(part)
    lbs = part.shared_buses_local

    if α === nothing
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            dθ = _wrap_angle_diff(θs[i] - θ̄_vec[j])
            dV = (Vs[i] - V̄_vec[j])
            dualsθ_vec_zone[j] += ρθ * dθ
            dualsV_vec_zone[j] += ρV * dV
        end
    else
        αv = something(α, 1.0)
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            θrel = _wrap_angle_diff(αv * θs[i] + (1 - αv) * θ̄_prev[j] - θ̄_vec[j])
            Vrel = αv * Vs[i] + (1 - αv) * V̄_prev[j] - V̄_vec[j]
            dualsθ_vec_zone[j] += ρθ * θrel
            dualsV_vec_zone[j] += ρV * Vrel
        end
    end
    return nothing
end

#############################
# DUAIS BRANCH AC           #
#############################

function update_duals_branch_ac!(
    dualsP::Dict{Int, Dict{Int, Float64}},
    dualsQ::Dict{Int, Dict{Int, Float64}},
    partitions::Dict{Int, PartitionAC},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ρP::Float64,
    ρQ::Float64;
    α::Union{Nothing, Float64} = nothing,
    zP_prev::Union{Nothing, Dict{Int, Float64}} = nothing,
    zQ_prev::Union{Nothing, Dict{Int, Float64}} = nothing,
)
    for (z, part) in partitions
        m = part.model::JuMP.Model
        λPz = get!(dualsP, z, Dict{Int, Float64}())
        λQz = get!(dualsQ, z, Dict{Int, Float64}())

        for e in part.zone_edges
            l_local = get(part.edge_to_local, e, nothing)
            l_local === nothing && continue
            P_loc, Q_loc = _ac_flow_from_exprs_existing(part, l_local)
            pval = JuMP.value(P_loc)
            qval = JuMP.value(Q_loc)

            if α === nothing
                if isfinite(pval)
                    λPz[e] = get(λPz, e, 0.0) + ρP * (pval - get(zP, e, 0.0))
                end
                if isfinite(qval)
                    λQz[e] = get(λQz, e, 0.0) + ρQ * (qval - get(zQ, e, 0.0))
                end
            else
                αv = something(α, 1.0)
                if isfinite(pval)
                    p_rel = αv * pval + (1 - αv) * get(zP_prev, e, get(zP, e, 0.0))
                    λPz[e] = get(λPz, e, 0.0) + ρP * (p_rel - get(zP, e, 0.0))
                end
                if isfinite(qval)
                    q_rel = αv * qval + (1 - αv) * get(zQ_prev, e, get(zQ, e, 0.0))
                    λQz[e] = get(λQz, e, 0.0) + ρQ * (q_rel - get(zQ, e, 0.0))
                end
            end
        end
    end
    return nothing
end

#############################
# RESÍDUOS BUS AC           #
#############################

function primal_residual_bus_ac_vec(
    partitions::Dict{Int, PartitionAC},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    acc = 0.0
    for (_, part) in partitions
        θbuf = snapshot_theta_ac!(part)
        Vbuf = snapshot_vmag!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            j = b2pos[b]
            acc +=
                _wrap_angle_diff(θbuf[pos[b]] - θ̄_vec[j])^2 +
                (Vbuf[pos[b]] - V̄_vec[j])^2
        end
    end
    return sqrt(acc)
end

function primal_residual_bus_ac_components(
    partitions::Dict{Int, PartitionAC},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    rθ = 0.0
    rV = 0.0
    for (_, part) in partitions
        θbuf = snapshot_theta_ac!(part)
        Vbuf = snapshot_vmag!(part)
        pos = part.bus_local_index
        @inbounds for b in part.shared_buses_local
            j = b2pos[b]
            rθ += _wrap_angle_diff(θbuf[pos[b]] - θ̄_vec[j])^2
            rV += (Vbuf[pos[b]] - V̄_vec[j])^2
        end
    end
    return sqrt(rθ), sqrt(rV)
end

function eps_primal_bus_ac_vec(
    partitions::Dict{Int, PartitionAC},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim = 0
    for (_, part) in partitions
        dim += 2 * length(part.shared_buses_local)
    end
    z_norm_sq = 0.0
    @inbounds for j in eachindex(shared_buses_sorted)
        z_norm_sq += θ̄_vec[j]^2 + V̄_vec[j]^2
    end
    return sqrt(max(1, dim)) * ϵ_abs +
           ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function dual_residual_bus_ac_vec(
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    θ̄_prev::Vector{Float64},
    V̄_prev::Vector{Float64},
)
    sθ = 0.0
    sV = 0.0
    @inbounds for j in eachindex(θ̄_vec)
        sθ += _wrap_angle_diff(θ̄_vec[j] - θ̄_prev[j])^2
        sV += (V̄_vec[j] - V̄_prev[j])^2
    end
    return sqrt(sθ + sV)
end

function dual_residual_bus_ac_components(
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    θ̄_prev::Vector{Float64},
    V̄_prev::Vector{Float64},
)
    sθ = 0.0
    sV = 0.0
    @inbounds for j in eachindex(θ̄_vec)
        sθ += _wrap_angle_diff(θ̄_vec[j] - θ̄_prev[j])^2
        sV += (V̄_vec[j] - V̄_prev[j])^2
    end
    return sqrt(sθ), sqrt(sV)
end

function eps_dual_bus_ac_vec(
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    θ̄_prev::Vector{Float64},
    V̄_prev::Vector{Float64},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim = 2 * max(1, length(θ̄_vec))
    s = dual_residual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
end

#############################
# ρ ADAPTATIVO AC           #
#############################

function adaptive_rho_ac!(
    ρ::Base.RefValue{Float64},
    r::Float64,
    s::Float64;
    μ_incr = 2.0,
    μ_decr = 2.0,
    τ_incr = 1.05,
    τ_decr = 0.98,
)
    if r > μ_incr * s
        ρ[] *= τ_incr
    elseif s > μ_decr * r
        ρ[] *= τ_decr
    end
end

#############################
# RESÍDUOS BRANCH AC        #
#############################

function primal_residual_branch_ac(
    partitions::Dict{Int, PartitionAC},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    acc = 0.0
    for (_, part) in partitions
        m = part.model::JuMP.Model
        for e in part.zone_edges
            l_local = get(part.edge_to_local, e, nothing)
            l_local === nothing && continue
            P_loc, Q_loc = _ac_flow_from_exprs_existing(part, l_local)
            p = JuMP.value(P_loc)
            q = JuMP.value(Q_loc)
            if isfinite(p)
                acc += (p - get(zP, e, 0.0))^2
            end
            if isfinite(q)
                acc += (q - get(zQ, e, 0.0))^2
            end
        end
    end
    return sqrt(acc)
end

function dual_residual_branch_ac(
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    zP_prev::Dict{Int, Float64},
    zQ_prev::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    s = 0.0
    for e in cut_edges
        s += (get(zP, e, 0.0) - get(zP_prev, e, 0.0))^2
        s += (get(zQ, e, 0.0) - get(zQ_prev, e, 0.0))^2
    end
    return sqrt(s)
end

function eps_primal_branch_ac(
    partitions::Dict{Int, PartitionAC},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim = 0
    for (_, part) in partitions
        dim += 2 * length(part.zone_edges)
    end
    z_norm_sq = 0.0
    for e in cut_edges
        z_norm_sq += get(zP, e, 0.0)^2 + get(zQ, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * ϵ_abs +
           ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_branch_ac(
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    zP_prev::Dict{Int, Float64},
    zQ_prev::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
end


#################################
# ADMM: BUS (AC) – usando sys   #
#################################

function run_admm_ac(config::DistributedOPFModel, sys::ADMMSystem)
    t0 = time_ns()
    @assert config.base_model == AC "ADMM-AC: base_model deve ser AC."

    partitions = build_ac_partitions(sys)
    shared_buses = sys.shared_buses

    ρθ = Ref(config.rho_theta)
    ρV = Ref(config.rho_v)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)

    for (_, part) in partitions
        prepare_hotpaths_ac!(part)
    end
    shared_buses_sorted, b2pos = prepare_dense_maps(shared_buses)

    dualsθ_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    dualsV_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )

    θ̄_vec = zeros(Float64, length(shared_buses_sorted))
    V̄_vec = ones(Float64, length(shared_buses_sorted))
    θ̄_prev = similar(θ̄_vec)
    V̄_prev = similar(V̄_vec)

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rhoθ_hist = Float64[]
    rhoV_hist = Float64[]
    θ_avg_hist = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_ac(config)

    # build inicial + termo quadrático em bus
    for (zone_id, part) in partitions
        sys_zone = sys.parts[zone_id]

        m = build_model(
            sys_zone,
            AC;
            optimizer = optimizer_factory,
        )

        part.model = m
        part.obj_base = m[:generation_cost]
        part.theta_var = m[:theta]
        part.vmag_var = m[:v_magnitude]
        part.has_solution = false

        _set_model_silent!(m)
        for I in keys(part.theta_var)
            set_start_value(part.theta_var[I], 0.0)
        end
        for I in keys(part.vmag_var)
            set_start_value(part.vmag_var[I], 1.0)
        end

        _build_quad_bus_expr_ac!(part, ρθ[], ρV[])
        @objective(m, Min, part.obj_base + part.quad_bus_expr_ac)
    end

    prev_ρθ, prev_ρV = ρθ[], ρV[]
    zone_ids = collect(keys(partitions))
    λθ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )
    λV_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    iter_done = 0

    for k in 1:config.max_iter
        iter_done = k
        _FAST_USE_LOGS && println("\n[ADMM-AC] Iteração $k")

        if ρθ[] != prev_ρθ || ρV[] != prev_ρV
            for (_, part) in partitions
                _build_quad_bus_expr_ac!(part, ρθ[], ρV[])
            end
            prev_ρθ, prev_ρV = ρθ[], ρV[]
        end

        # sem over-relaxation / FISTA
        α_used = nothing
        θ̄_target = θ̄_vec
        V̄_target = V̄_vec

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            λθ = (λθ_zone_buf[z] .= dualsθ_dense[z])
            λV = (λV_zone_buf[z] .= dualsV_dense[z])

            lin = _build_linear_bus_expr_ac!(
                part,
                shared_buses_sorted,
                b2pos,
                θ̄_target,
                V̄_target,
                λθ,
                λV,
                ρθ[],
                ρV[],
            )

            m = part.model::JuMP.Model
            @objective(m, Min, part.obj_base + part.quad_bus_expr_ac + lin)
            optimize!(m)

            results[i] = (
                expr_cost = value(m[:generation_cost]),
                real_cost = compute_real_generation_cost(m, part.objective_params),
            )
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.expr_cost
            part.real_generation_cost = out.real_cost
        end

        θ̄_prev .= θ̄_vec
        V̄_prev .= V̄_vec
        compute_ac_bus_consensus_vec!(
            θ̄_vec,
            V̄_vec,
            partitions,
            shared_buses_sorted,
            b2pos,
        )
        @inbounds for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_hist, b, Float64[]), θ̄_vec[j])
        end

        for z in zone_ids
            update_duals_bus_ac_vec!(
                dualsθ_dense[z],
                dualsV_dense[z],
                partitions[z],
                θ̄_vec,
                V̄_vec,
                ρθ[],
                ρV[],
                shared_buses_sorted,
                b2pos;
                α = α_used,
                θ̄_prev = θ̄_prev,
                V̄_prev = V̄_prev,
            )
        end

        r = primal_residual_bus_ac_vec(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
        )
        s = dual_residual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev)
        εp = eps_primal_bus_ac_vec(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
            ϵ_abs,
            ϵ_rel,
        )
        εd = eps_dual_bus_ac_vec(
            θ̄_vec,
            V̄_vec,
            θ̄_prev,
            V̄_prev,
            ϵ_abs,
            ϵ_rel,
        )

        push!(res_pri_hist, r)
        push!(res_du_hist, s)
        push!(eps_pri_hist, εp)
        push!(eps_du_hist, εd)

        _FAST_USE_LOGS && println("‖r‖ = $r | ε_pri = $εp")
        _FAST_USE_LOGS && println("‖s‖ = $s | ε_du  = $εd")

        push!(rhoθ_hist, ρθ[])
        push!(rhoV_hist, ρV[])

        rθ, rV = primal_residual_bus_ac_components(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
        )
        sθ, sV = dual_residual_bus_ac_components(
            θ̄_vec,
            V̄_vec,
            θ̄_prev,
            V̄_prev,
        )

        if r < εp && s < εd
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_ac!(
                ρθ,
                rθ,
                sθ;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
            adaptive_rho_ac!(
                ρV,
                rV,
                sV;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
        end
    end

    save_ac_bus_plots(
        res_pri_hist,
        res_du_hist,
        eps_pri_hist,
        eps_du_hist,
        rhoθ_hist,
        rhoV_hist,
        θ_avg_hist;
        stem = "admm_ac",
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    println("\n[ADMM-AC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ_final = $(ρθ[]) | ρV_final = $(ρV[])")
    elapsed_s = (time_ns() - t0) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (bus/AC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end

#################################
# ADMM: BRANCH (AC) – sys       #
#################################

function run_admm_ac_branch(
    config::DistributedOPFModel,
    sys::ADMMSystem;
    warm_start::Bool = false,
)
    t0 = time_ns()
    @assert config.base_model == AC "ADMM-AC-branch: base_model deve ser AC."

    partitions = build_ac_partitions(sys)
    cut_edges = sys.cut_edges

    for (_, part) in partitions
        part.zone_edges = _zone_cut_edges(part, cut_edges)
    end

    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)

    dualsP =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]

    optimizer_factory = resolve_optimizer_ac(config)

    # build inicial + termo quadrático em branch
    for (zone_id, part) in partitions
        sys_zone = sys.parts[zone_id]

        m = build_model(
            sys_zone,
            AC;
            optimizer = optimizer_factory,
        )
        part.model = m
        part.obj_base = m[:generation_cost]
        part.theta_var = m[:theta]
        part.vmag_var = m[:v_magnitude]
        part.has_solution = false

        _set_model_silent!(m)
        for I in keys(part.theta_var)
            set_start_value(part.theta_var[I], 0.0)
        end
        for I in keys(part.vmag_var)
            set_start_value(part.vmag_var[I], 1.0)
        end

        _build_quad_branch_expr_ac!(part, ρP[], ρQ[], part.zone_edges)
        @objective(m, Min, part.obj_base + part.quad_branch_expr_ac)
    end

    prev_ρP, prev_ρQ = ρP[], ρQ[]
    zone_ids = collect(keys(partitions))
    iter_done = 0

    for k in 1:config.max_iter
        iter_done = k
        _FAST_USE_LOGS && println("\n[ADMM-AC-branch] Iteração $k")

        if ρP[] != prev_ρP || ρQ[] != prev_ρQ
            for (_, part) in partitions
                _build_quad_branch_expr_ac!(part, ρP[], ρQ[], part.zone_edges)
            end
            prev_ρP, prev_ρQ = ρP[], ρQ[]
        end

        # sem over-relaxation
        α_used = nothing

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]
            lin = _build_linear_branch_expr_ac!(
                part,
                part.zone_edges,
                zP,
                zQ,
                dualsP[z],
                dualsQ[z],
                ρP[],
                ρQ[],
            )
            m = part.model::JuMP.Model
            @objective(m, Min, part.obj_base + part.quad_branch_expr_ac + lin)
            optimize!(m)

            results[i] = (
                expr_cost = value(m[:generation_cost]),
                real_cost = compute_real_generation_cost(m, part.objective_params),
            )
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.expr_cost
            part.real_generation_cost = out.real_cost
        end

        zP_prev = copy(zP)
        zQ_prev = copy(zQ)
        compute_ac_branch_consensus!(zP, zQ, partitions, cut_edges)

        update_duals_branch_ac!(
            dualsP,
            dualsQ,
            partitions,
            zP,
            zQ,
            cut_edges,
            ρP[],
            ρQ[];
            α = α_used,
            zP_prev = zP_prev,
            zQ_prev = zQ_prev,
        )

        r = primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
        s = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
        εp = eps_primal_branch_ac(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
        εd = eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)

        push!(res_pri_hist, r)
        push!(res_du_hist, s)
        push!(eps_pri_hist, εp)
        push!(eps_du_hist, εd)

        _FAST_USE_LOGS && println("‖r‖ = $r | ε_pri = $εp")
        _FAST_USE_LOGS && println("‖s‖ = $s | ε_du  = $εd")

        push!(rhoP_hist, ρP[])
        push!(rhoQ_hist, ρQ[])

        if r < εp && s < εd
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_ac!(
                ρP,
                r,
                s;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
            adaptive_rho_ac!(
                ρQ,
                r,
                s;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
        end
    end

    save_ac_branch_plots(
        res_pri_hist,
        res_du_hist,
        eps_pri_hist,
        eps_du_hist,
        rhoP_hist,
        rhoQ_hist;
        stem = "admm_ac_branch",
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    println("\n[ADMM-AC-branch] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρP_final = $(ρP[]) | ρQ_final = $(ρQ[])")
    elapsed_s = (time_ns() - t0) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (branch/AC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end

#################################
# ADMM: HÍBRIDO (AC) – sys      #
#################################

function run_admm_ac_hybrid(config::DistributedOPFModel, sys::ADMMSystem)
    t0 = time_ns()
    @assert config.base_model == AC "Hybrid AC: base_model deve ser AC."

    partitions = build_ac_partitions(sys)
    shared_buses = sys.shared_buses
    cut_edges = sys.cut_edges

    for (_, part) in partitions
        part.zone_edges = _zone_cut_edges(part, cut_edges)
        prepare_hotpaths_ac!(part)
    end

    ρθ = Ref(config.rho_theta)
    ρV = Ref(config.rho_v)
    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)

    adaptive_bus = adaptive
    adaptive_br = adaptive

    shared_buses_sorted, b2pos = prepare_dense_maps(shared_buses)

    dualsθ_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    dualsV_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    dualsP =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    θ̄_vec = zeros(Float64, length(shared_buses_sorted))
    V̄_vec = ones(Float64, length(shared_buses_sorted))
    θ̄_prev = similar(θ̄_vec)
    V̄_prev = similar(V̄_vec)

    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rhoθ_hist = Float64[]
    rhoV_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]
    θ_avg_hist = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_ac(config)

    # build inicial + termos quadráticos
    for (zone_id, part) in partitions
        sys_zone = sys.parts[zone_id]
        m = build_model(
            sys_zone,
            AC;
            optimizer = optimizer_factory,
        )
        part.model = m
        part.obj_base = m[:generation_cost]
        part.theta_var = m[:theta]
        part.vmag_var = m[:v_magnitude]
        part.has_solution = false

        _set_model_silent!(m)
        for I in keys(part.theta_var)
            set_start_value(part.theta_var[I], 0.0)
        end
        for I in keys(part.vmag_var)
            set_start_value(part.vmag_var[I], 1.0)
        end

        _build_quad_bus_expr_ac!(part, ρθ[], ρV[])
        _build_quad_branch_expr_ac!(part, ρP[], ρQ[], part.zone_edges)
        @objective(
            m,
            Min,
            part.obj_base + part.quad_bus_expr_ac + part.quad_branch_expr_ac,
        )
    end

    prev_ρθ, prev_ρV = ρθ[], ρV[]
    prev_ρP, prev_ρQ = ρP[], ρQ[]
    zone_ids = collect(keys(partitions))
    λθ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )
    λV_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    iter_done = 0

    for k in 1:config.max_iter
        iter_done = k
        _FAST_USE_LOGS && println("\n[ADMM-AC-hybrid] Iteração $k")

        if ρθ[] != prev_ρθ || ρV[] != prev_ρV
            for (_, part) in partitions
                _build_quad_bus_expr_ac!(part, ρθ[], ρV[])
            end
            prev_ρθ, prev_ρV = ρθ[], ρV[]
        end
        if ρP[] != prev_ρP || ρQ[] != prev_ρQ
            for (_, part) in partitions
                _build_quad_branch_expr_ac!(part, ρP[], ρQ[], part.zone_edges)
            end
            prev_ρP, prev_ρQ = ρP[], ρQ[]
        end

        # sem over-relaxation / aceleração
        α_used = nothing
        θ̄_target = θ̄_vec
        V̄_target = V̄_vec

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            λθ = (λθ_zone_buf[z] .= dualsθ_dense[z])
            λV = (λV_zone_buf[z] .= dualsV_dense[z])

            lin_bus = _build_linear_bus_expr_ac!(
                part,
                shared_buses_sorted,
                b2pos,
                θ̄_target,
                V̄_target,
                λθ,
                λV,
                ρθ[],
                ρV[],
            )
            lin_br = _build_linear_branch_expr_ac!(
                part,
                part.zone_edges,
                zP,
                zQ,
                dualsP[z],
                dualsQ[z],
                ρP[],
                ρQ[],
            )

            m = part.model::JuMP.Model
            @objective(
                m,
                Min,
                part.obj_base + part.quad_bus_expr_ac + part.quad_branch_expr_ac +
                (lin_bus + lin_br),
            )
            optimize!(m)

            results[i] = (
                expr_cost = value(m[:generation_cost]),
                real_cost = compute_real_generation_cost(m, part.objective_params),
            )
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.expr_cost
            part.real_generation_cost = out.real_cost
        end

        θ̄_prev .= θ̄_vec
        V̄_prev .= V̄_vec
        zP_prev = copy(zP)
        zQ_prev = copy(zQ)

        compute_ac_bus_consensus_vec!(
            θ̄_vec,
            V̄_vec,
            partitions,
            shared_buses_sorted,
            b2pos,
        )
        for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_hist, b, Float64[]), θ̄_vec[j])
        end
        compute_ac_branch_consensus!(zP, zQ, partitions, cut_edges)

        for z in zone_ids
            update_duals_bus_ac_vec!(
                dualsθ_dense[z],
                dualsV_dense[z],
                partitions[z],
                θ̄_vec,
                V̄_vec,
                ρθ[],
                ρV[],
                shared_buses_sorted,
                b2pos;
                α = α_used,
                θ̄_prev = θ̄_prev,
                V̄_prev = V̄_prev,
            )
        end
        update_duals_branch_ac!(
            dualsP,
            dualsQ,
            partitions,
            zP,
            zQ,
            cut_edges,
            ρP[],
            ρQ[];
            α = α_used,
            zP_prev = zP_prev,
            zQ_prev = zQ_prev,
        )

        r_bus = primal_residual_bus_ac_vec(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
        )
        r_br = primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
        r_tot = hypot(r_bus, r_br)

        s_bus = dual_residual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev)
        s_br = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
        s_tot = hypot(s_bus, s_br)

        εb_pri = eps_primal_bus_ac_vec(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
            ϵ_abs,
            ϵ_rel,
        )
        εe_pri = eps_primal_branch_ac(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
        ε_pri = hypot(εb_pri, εe_pri)

        εb_du = eps_dual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev, ϵ_abs, ϵ_rel)
        εe_du = eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
        ε_du = hypot(εb_du, εe_du)

        push!(res_pri_hist, r_tot)
        push!(res_du_hist, s_tot)
        push!(eps_pri_hist, ε_pri)
        push!(eps_du_hist, ε_du)

        _FAST_USE_LOGS && println("‖r‖ = $r_tot | ε_pri = $ε_pri")
        _FAST_USE_LOGS && println("‖s‖ = $s_tot | ε_du  = $ε_du")

        push!(rhoθ_hist, ρθ[])
        push!(rhoV_hist, ρV[])
        push!(rhoP_hist, ρP[])
        push!(rhoQ_hist, ρQ[])

        if r_tot < ε_pri && s_tot < ε_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive_bus
            rθ, rV = primal_residual_bus_ac_components(
                partitions,
                θ̄_vec,
                V̄_vec,
                shared_buses_sorted,
                b2pos,
            )
            sθ, sV = dual_residual_bus_ac_components(
                θ̄_vec,
                V̄_vec,
                θ̄_prev,
                V̄_prev,
            )
            adaptive_rho_ac!(
                ρθ,
                rθ,
                sθ;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
            adaptive_rho_ac!(
                ρV,
                rV,
                sV;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
        end
        if adaptive_br
            adaptive_rho_ac!(
                ρP,
                r_br,
                s_br;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
            adaptive_rho_ac!(
                ρQ,
                r_br,
                s_br;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
        end
    end

    plot_admm_convergence(
        res_pri_hist,
        res_du_hist,
        eps_pri_hist,
        eps_du_hist,
        "plot_convergencia_admm_ac_hybrid.html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoθ_hist),
                    y = rhoθ_hist,
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
        "plot_rho_theta_ac_hybrid.html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoV_hist),
                    y = rhoV_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρV",
                );
                p.layout.title = "Evolução de ρV";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρV";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_v_ac_hybrid.html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoP_hist),
                    y = rhoP_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρP",
                );
                p.layout.title = "Evolução de ρP";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρP";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_p_ac_hybrid.html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoQ_hist),
                    y = rhoQ_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρQ",
                );
                p.layout.title = "Evolução de ρQ";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρQ";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_q_ac_hybrid.html",
    )

    plot_theta_avg_heatmap(θ_avg_hist, "plot_theta_avg_heatmap_ac_hybrid.html")

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    println("\n[ADMM-AC-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println(
        "ρθ_final = $(ρθ[]) | ρV_final = $(ρV[]) | ρP_final = $(ρP[]) | ρQ_final = $(ρQ[])",
    )
    elapsed_s = (time_ns() - t0) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (hybrid/AC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end