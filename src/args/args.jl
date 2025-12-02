function finish_path(path::String)
    return abspath(path)
end

struct Args
    database_path::String
    model_type::Union{Type{<:OPFModel}, DistributedOPFModel}
    optimizer::Any
    write_outputs::Bool

    function Args(raw_args::Vector{String})
        s = ArgParseSettings()

        @add_arg_table! s begin
            "path"
            help = "Path to MATPOWER .m file or a CSV folder"
            arg_type = String
            required = true

            "--model"
            help = "Model: dc | ac | socp | distributed-admm-dc | distributed-admm-ac | distributed-admm-socp | distributed-admm-branch-dc | distributed-admm-branch-ac | distributed-admm-branch-socp | distributed-admm-hybrid-dc | distributed-admm-hybrid-ac | distributed-admm-hybrid-socp"
            arg_type = String
            required = true

            "--optimizer"
            help = "Optimizer function, e.g., OSQP.Optimizer (default if omitted), or Ipopt.Optimizer / SCS.Optimizer / Gurobi.Optimizer, etc."
            arg_type = String

            "--outputs"
            help = "Generate CSV output files"
            action = :store_true

            # ======= ADMM =======
            "--admm-rho"
            help = "Default ρ used if a specific rho (theta/v/flow) is not provided (default: 1000.0)"
            arg_type = Float64
            default = 1_000.0

            "--admm-rho-theta"
            help = "ρθ (consensus penalty for bus angles; ignored for SOCP-hybrid)"
            arg_type = Float64

            "--admm-rho-v"
            help = "ρV (consensus penalty for voltage magnitudes; AC usa |V|, SOCP usa V^2)"
            arg_type = Float64

            "--admm-rho-flow"
            help = "ρe (consensus penalty for branch flows; applies to P and Q if specifics are absent)"
            arg_type = Float64

            "--admm-rho-flow-p"
            help = "ρP (consensus penalty for active power flow in branch consensus; em SOCP-hybrid também usado para R=Re(W))"
            arg_type = Float64

            "--admm-rho-flow-q"
            help = "ρQ (consensus penalty for reactive power flow in branch consensus; em SOCP-hybrid também usado para T=Im(W))"
            arg_type = Float64

            "--admm-max-iter"
            help = "Max ADMM iterations"
            arg_type = Int
            default = 350

            "--admm-eps-primal"
            help = "Primal tolerance ε_primal"
            arg_type = Float64
            default = 1e-7

            "--admm-eps-dual"
            help = "Dual tolerance ε_dual"
            arg_type = Float64
            default = 1e-7

            "--admm-alpha"
            help = "Over-relaxation α (use only if 1<α<2). If omitted, over-relaxation is disabled."
            arg_type = Float64

            "--admm-warm-start"
            help = "Warm start dos subproblemas (pré-solve base e usa θ como start)"
            action = :store_true

            "--admm-adaptive"
            help = "Liga/desliga o ajuste adaptativo de ρ no ADMM: on | off (padrão: on)"
            arg_type = String
            default = "on"

            # ----- ADMM: parameters -----
            "--admm-mu-incr"
            help = "Fator μ_incr para comparar r vs s (default: 3.0)"
            arg_type = Float64
            default = 3.0

            "--admm-mu-decr"
            help = "Fator μ_decr para comparar s vs r (default: 3.0)"
            arg_type = Float64
            default = 3.0

            "--admm-tau-incr"
            help = "Fator multiplicativo τ_incr quando r >> s (default: 1.05)"
            arg_type = Float64
            default = 1.05

            "--admm-tau-decr"
            help = "Fator multiplicativo τ_decr quando s >> r (default: 0.98)"
            arg_type = Float64
            default = 0.98

            # ----- ADMM: Nesterov -----
            "--admm-accel"
            help = "Habilita aceleração tipo Nesterov no ADMM (padrão: desligado)"
            action = :store_true
        end

        parsed = parse_args(raw_args, s)

        path = finish_path(parsed["path"])
        if !(isfile(path) || isdir(path))
            error("Invalid path: $path is not a file or directory.")
        end

        model_key = lowercase(strip(parsed["model"]))
        base_map = Dict(
            "dc"   => DC,
            "ac"   => AC,
            "socp" => SOCP,
        )

        # ===== rhos =====
        rho_base    = parsed["admm-rho"]
        rho_theta_  = get(parsed, "admm-rho-theta", nothing)
        rho_v_      = get(parsed, "admm-rho-v", nothing)
        rho_flow_   = get(parsed, "admm-rho-flow", nothing)
        rho_flow_p_ = get(parsed, "admm-rho-flow-p", nothing)
        rho_flow_q_ = get(parsed, "admm-rho-flow-q", nothing)

        # fallback: specific > generic flow > base
        rho_theta_val  = isnothing(rho_theta_)  ? rho_base : rho_theta_
        rho_v_val      = isnothing(rho_v_)      ? rho_base : rho_v_
        rho_flow_p_val = isnothing(rho_flow_p_) ? (isnothing(rho_flow_) ? rho_base : rho_flow_) : rho_flow_p_
        rho_flow_q_val = isnothing(rho_flow_q_) ? (isnothing(rho_flow_) ? rho_base : rho_flow_) : rho_flow_q_

        # ===== alpha (over-relaxation) =====
        alpha_ = get(parsed, "admm-alpha", nothing)
        if alpha_ !== nothing && !(1.0 < alpha_ < 2.0)
            @warn "Ignoring --admm-alpha=$(alpha_) because it must satisfy 1<α<2."
            alpha_ = nothing
        end

        # ===== adaptive rho (on/off) =====
        adaptive_flag_str = lowercase(strip(parsed["admm-adaptive"]))
        adaptive_flag = adaptive_flag_str == "on"  ? true  : adaptive_flag_str == "off" ? false : ( @warn "Valor inválido para --admm-adaptive=$(parsed["admm-adaptive"]). Usando 'on'."; true )
        mu_incr_val = parsed["admm-mu-incr"]
        mu_decr_val = parsed["admm-mu-decr"]
        tau_incr_val = parsed["admm-tau-incr"]
        tau_decr_val = parsed["admm-tau-decr"]

        # ===== acceleration (on/off) =====
        accel_flag = get(parsed, "admm-accel", false)

        # ===== optimizer =====
        optimizer = if haskey(parsed, "optimizer") && !isnothing(parsed["optimizer"])
            try
                eval(Meta.parse(parsed["optimizer"]))
            catch
                @warn "Invalid optimizer: $(parsed["optimizer"]). Falling back to default (OSQP) via resolve_optimizer_dc."
                nothing
            end
        else
            nothing
        end
        local function _warn_irrelevant_rhos(
            method_sym::Symbol,
            base_str::AbstractString,
            model_key_local::AbstractString,
            rho_theta_orig::Union{Nothing, Float64},
            rho_v_orig::Union{Nothing, Float64},
            rho_flow_orig::Union{Nothing, Float64},
            rho_flow_p_orig::Union{Nothing, Float64},
            rho_flow_q_orig::Union{Nothing, Float64},
        )
            # θ em :admm(bus) para DC/AC/SOCP; em :admm_hybrid só para DC/AC
            uses_theta =
                (
                    method_sym == :admm &&
                    (base_str == "dc" || base_str == "ac" || base_str == "socp")
                ) ||
                (method_sym == :admm_hybrid && (base_str == "dc" || base_str == "ac"))
            # V (|V| em AC, V^2 em SOCP) em :admm e :admm_hybrid para AC/SOCP
            uses_v =
                (method_sym == :admm || method_sym == :admm_hybrid) &&
                (base_str == "ac" || base_str == "socp")
            # flow em :admm_branch e :admm_hybrid
            uses_flow = (method_sym == :admm_branch) || (method_sym == :admm_hybrid)

            if rho_theta_orig !== nothing && !uses_theta
                @warn "Flag --admm-rho-theta fornecida mas ignorada para '$model_key_local'."
            end
            if rho_v_orig !== nothing && !uses_v
                @warn "Flag --admm-rho-v fornecida mas ignorada para '$model_key_local'."
            end
            if !(uses_flow)
                if rho_flow_orig !== nothing
                    @warn "Flag --admm-rho-flow fornecida mas ignorada para '$model_key_local'."
                end
                if rho_flow_p_orig !== nothing
                    @warn "Flag --admm-rho-flow-p fornecida mas ignorada para '$model_key_local'."
                end
                if rho_flow_q_orig !== nothing
                    @warn "Flag --admm-rho-flow-q fornecida mas ignorada para '$model_key_local'."
                end
            end
            return nothing
        end

        model_type = begin
            if startswith(model_key, "distributed-")
                parts = split(model_key, '-')
                if length(parts) == 3
                    _, method_str, base_str = parts
                    variant_str = ""
                elseif length(parts) == 4
                    _, method_str, variant_str, base_str = parts
                else
                    error("Invalid distributed model key: $(parsed["model"])")
                end

                method_str  = String(method_str)
                base_str    = String(base_str)
                variant_str = String(variant_str)

                method_str == "admm" ||
                    error("Only ADMM is supported for distributed-* now.")
                haskey(base_map, base_str) ||
                    error("Invalid base model in distributed-* key: $base_str")

                method_sym =
                    variant_str == ""       ? :admm :
                    variant_str == "branch" ? :admm_branch :
                    variant_str == "hybrid" ? :admm_hybrid :
                    error("Unsupported distributed method variant: $variant_str")

                _warn_irrelevant_rhos(
                    method_sym, base_str, model_key,
                    rho_theta_, rho_v_, rho_flow_, rho_flow_p_, rho_flow_q_,
                )

                warm_start_flag = get(parsed, "admm-warm-start", false)

                local cfg
                try
                    cfg = DistributedOPFModel(
                        base_map[base_str],
                        method_sym,
                        parsed["admm-max-iter"];
                        rho = rho_base,
                        rho_theta = rho_theta_val,
                        rho_v = rho_v_val,
                        rho_flow_p = rho_flow_p_val,
                        rho_flow_q = rho_flow_q_val,
                        ε_primal = parsed["admm-eps-primal"],
                        ε_dual   = parsed["admm-eps-dual"],
                        alpha_overrelax = alpha_,
                        warm_start  = warm_start_flag,
                        adaptive_rho = adaptive_flag,
                        mu_incr = mu_incr_val,
                        mu_decr = mu_decr_val,
                        tau_incr = tau_incr_val,
                        tau_decr = tau_decr_val,
                        accel_nesterov = accel_flag,
                        optimizer = optimizer,
                        prefer_socp_solver = :auto,
                    )
                catch err
                    @warn "Esta versão de DistributedOPFModel não aceita alguns campos novos (warm_start/adaptive_rho/accel_nesterov). Prosseguindo sem eles." err
                    cfg = DistributedOPFModel(
                        base_map[base_str],
                        method_sym,
                        parsed["admm-max-iter"];
                        rho = rho_base,
                        rho_theta = rho_theta_val,
                        rho_v = rho_v_val,
                        rho_flow_p = rho_flow_p_val,
                        rho_flow_q = rho_flow_q_val,
                        ε_primal = parsed["admm-eps-primal"],
                        ε_dual   = parsed["admm-eps-dual"],
                        alpha_overrelax = alpha_,
                        optimizer = optimizer,
                        prefer_socp_solver = :auto,
                    )
                end
                cfg
            else
                haskey(base_map, model_key) || error("Invalid model: $model_key")
                base_map[String(model_key)]
            end
        end

        return new(
            path,
            model_type,
            optimizer,
            get(parsed, "outputs", false),
        )
    end
end