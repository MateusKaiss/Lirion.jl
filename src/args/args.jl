function finish_path(path::String)
    return abspath(path)
end

struct Args
    database_path::String
    model_type::Union{Type{<:OPFModel}, DistributedOPFModel}
    optimizer::Any
    plot_outputs::Bool
    write_lp::Bool
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
            help = "Optimizer function, e.g., Ipopt.Optimizer or SCS.Optimizer or Gurobi.Optimizer (default: Ipopt.Optimizer)"
            arg_type = String

            "--plot-outputs"
            help = "Generate plots"
            action = :store_true

            "--write-lp"
            help = "Write LP file exported by JuMP/optimizer (if supported)"
            action = :store_true

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
        end

        parsed = parse_args(raw_args, s)

        path = finish_path(parsed["path"])
        if !(isfile(path) || isdir(path))
            error("Invalid path: $path is not a file or directory.")
        end

        model_key = lowercase(strip(parsed["model"]))
        base_map = Dict(
            "dc" => DC,
            "ac" => AC,
            "socp" => SOCP,
        )

        # ===== rhos =====
        rho_base    = parsed["admm-rho"]
        rho_theta_  = get(parsed, "admm-rho-theta", nothing)
        rho_v_      = get(parsed, "admm-rho-v", nothing)
        rho_flow_   = get(parsed, "admm-rho-flow", nothing)
        rho_flow_p_ = get(parsed, "admm-rho-flow-p", nothing)
        rho_flow_q_ = get(parsed, "admm-rho-flow-q", nothing)

        # fallback: específicos > flow genérico > base
        rho_theta_val = isnothing(rho_theta_) ? rho_base : rho_theta_
        rho_v_val     = isnothing(rho_v_)     ? rho_base : rho_v_
        rho_flow_p_val =
            isnothing(rho_flow_p_) ? (isnothing(rho_flow_) ? rho_base : rho_flow_) :
            rho_flow_p_
        rho_flow_q_val =
            isnothing(rho_flow_q_) ? (isnothing(rho_flow_) ? rho_base : rho_flow_) :
            rho_flow_q_

        # avisos de flags irrelevantes/ignoradas
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
            # θ é usado em :admm(bus) para DC/AC/SOCP; em :admm_hybrid só para DC/AC (SOCP-hybrid NÃO usa θ)
            uses_theta =
                (method_sym == :admm  && (base_str == "dc" || base_str == "ac" || base_str == "socp")) ||
                (method_sym == :admm_hybrid && (base_str == "dc" || base_str == "ac"))
            # V (|V| em AC, V^2 em SOCP) é usado em :admm e :admm_hybrid para AC/SOCP
            uses_v =
                (method_sym == :admm || method_sym == :admm_hybrid) &&
                (base_str == "ac" || base_str == "socp")
            # flow (P/Q e afins) é usado em :admm_branch e :admm_hybrid para todos
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
                # formatos aceitos:
                #  - distributed-admm-<base>
                #  - distributed-admm-branch-<base>
                #  - distributed-admm-hybrid-<base>
                if length(parts) == 3
                    _, method_str, base_str = parts
                    variant_str = ""
                elseif length(parts) == 4
                    _, method_str, variant_str, base_str = parts
                else
                    error("Invalid distributed model key: $(parsed["model"])")
                end

                method_str = String(method_str)
                base_str   = String(base_str)
                variant_str = String(variant_str)

                method_str == "admm" || error("Only ADMM is supported for distributed-* now.")
                haskey(base_map, base_str) || error("Invalid base model in distributed-* key: $base_str")

                method_sym =
                    variant_str == ""        ? :admm :
                    variant_str == "branch"  ? :admm_branch :
                    variant_str == "hybrid"  ? :admm_hybrid :
                    error("Unsupported distributed method variant: $variant_str")

                _warn_irrelevant_rhos(
                    method_sym, base_str, model_key,
                    rho_theta_, rho_v_, rho_flow_, rho_flow_p_, rho_flow_q_,
                )

                DistributedOPFModel(
                    base_map[base_str],
                    method_sym,
                    parsed["admm-max-iter"];
                    rho         = rho_base,
                    rho_theta   = rho_theta_val,
                    rho_v       = rho_v_val,
                    rho_flow_p  = rho_flow_p_val,   # em SOCP-hybrid vale para P e R
                    rho_flow_q  = rho_flow_q_val,   # em SOCP-hybrid vale para Q e T
                    ε_primal    = parsed["admm-eps-primal"],
                    ε_dual      = parsed["admm-eps-dual"],
                )
            else
                haskey(base_map, model_key) || error("Invalid model: $model_key")
                base_map[String(model_key)]
            end
        end

        optimizer = if haskey(parsed, "optimizer") && !isnothing(parsed["optimizer"])
            try
                eval(Meta.parse(parsed["optimizer"]))
            catch
                @warn "Invalid optimizer: $(parsed["optimizer"]). Defaulting to Ipopt.Optimizer"
                Ipopt.Optimizer
            end
        else
            Ipopt.Optimizer
        end

        return new(
            path,
            model_type,
            optimizer,
            get(parsed, "plot-outputs", false),
            get(parsed, "write-lp", false),
            get(parsed, "outputs", false),
        )
    end
end