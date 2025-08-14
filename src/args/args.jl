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
            help = "Model: dc | ac | socp | distributed-admm-dc | distributed-admm-ac | distributed-admm-branch-dc | distributed-admm-hybrid-dc"
            arg_type = String
            required = true

            "--optimizer"
            help = "Optimizer function, e.g., Ipopt.Optimizer or Gurobi.Optimizer (default: Ipopt.Optimizer)"
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
            help = "Initial ρ (fallback if --admm-rho-theta/--admm-rho-v not provided)"
            arg_type = Float64
            default = 1_000.0

            "--admm-rho-theta"
            help = "Initial ρθ (consensus penalty for bus angles)"
            arg_type = Float64

            "--admm-rho-v"
            help = "Initial ρV (consensus penalty for voltage magnitudes)"
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

        rho_base = parsed["admm-rho"]
        rho_theta = get(parsed, "admm-rho-theta", nothing)
        rho_v     = get(parsed, "admm-rho-v", nothing)
        rho_theta = isnothing(rho_theta) ? rho_base : rho_theta
        rho_v     = isnothing(rho_v)     ? rho_base : rho_v

        model_type = begin
            if startswith(model_key, "distributed-")
                parts = split(model_key, '-')
                # accepted formats:
                #  - distributed-admm-dc
                #  - distributed-admm-ac
                #  - distributed-admm-branch-dc
                #  - distributed-admm-hybrid-dc
                if length(parts) == 3
                    _, method_str, base_str = parts
                    variant_str = ""  # ex.: "distributed-admm-dc" ou "distributed-admm-ac"
                elseif length(parts) == 4
                    _, method_str, variant_str, base_str = parts  # ex.: "distributed-admm-branch-dc"
                else
                    error("Invalid distributed model key: $(parsed["model"])")
                end

                method_str == "admm" ||
                    error("Only ADMM is supported for distributed-* now.")
                haskey(base_map, base_str) ||
                    error("Invalid base model in distributed-* key: $base_str")

                method_sym =
                    variant_str == ""       ? :admm :
                    variant_str == "branch" ? :admm_branch :
                    variant_str == "hybrid" ? :admm_hybrid :
                    error("Unsupported distributed method variant: $variant_str")

                if (method_sym == :admm_branch || method_sym == :admm_hybrid) && base_str != "dc"
                    error("distributed-admm-[branch|hybrid]-* currently supported only with DC")
                end

                DistributedOPFModel(
                    base_map[base_str],
                    method_sym,
                    parsed["admm-max-iter"],
                    rho_theta,             # << ρθ
                    rho_v,                 # << ρV
                    parsed["admm-eps-primal"],
                    parsed["admm-eps-dual"],
                )
            else
                haskey(base_map, model_key) || error("Invalid model: $model_key")
                base_map[model_key]
            end
        end

        # Optimizer (fallback: Ipopt.Optimizer)
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

# Examples of how DistributedOPFModel(DC, :admm, ...) was hardcoded before CLI parameters:
# DistributedOPFModel(DC, :admm, 350,       1000,       1e-7, 1e-7)       # default ADMM parameters
# DistributedOPFModel(DC, :admm, 700,       1000,       1e-7, 1e-7)       # more iterations
# DistributedOPFModel(DC, :admm, 350,       1000,       1e-7, 1e-7)       # duplicate default example
# DistributedOPFModel(DC, :admm, 500,       500_000,    1e-7, 1e-7)       # very large rho
# DistributedOPFModel(DC, :admm, 600,       1_000_000,  1e-7, 1e-7)       # extreme rho and iterations
