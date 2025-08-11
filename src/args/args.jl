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
            help = "Path to MATPOWER .m file or CSV folder"
            arg_type = String
            required = true

            "--model"
            help = "Model: dc | ac | socp | distributed-admm-dc | distributed-admm-ac | distributed-admm-socp"
            arg_type = String
            required = true

            "--optimizer"
            help = "Optimizer function (e.g., Ipopt.Optimizer, Gurobi.Optimizer)"
            arg_type = String

            "--plot-outputs"
            help = "Generate plots"
            action = :store_true

            "--write-lp"
            help = "Write LP file"
            action = :store_true

            "--outputs"
            help = "Generate CSV output files"
            action = :store_true

            # ======= ADMM =======
            "--admm-rho"
            help = "Initial ρ (only for distributed-admm-*)"
            arg_type = Float64
            default = 1_000.0

            "--admm-max-iter"
            help = "Max ADMM iterations (only for distributed-admm-*)"
            arg_type = Int
            default = 350

            "--admm-eps-primal"
            help = "Primal tolerance ε_primal (only for distributed-admm-*)"
            arg_type = Float64
            default = 1e-7

            "--admm-eps-dual"
            help = "Dual tolerance ε_dual (only for distributed-admm-*)"
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

        model_type = begin
            if startswith(model_key, "distributed-")
                parts = split(model_key, '-')
                if length(parts) != 3
                    error("Invalid distributed model key: $(parsed["model"])")
                end
                method_str = parts[2]  # "admm"
                base_str = parts[3]  # "dc" | "ac" | "socp"

                method_str == "admm" ||
                    error("Only ADMM is supported for distributed-* now.")
                haskey(base_map, base_str) ||
                    error("Invalid base model in distributed-* key: $base_str")

                DistributedOPFModel(
                    base_map[base_str],
                    :admm,
                    parsed["admm-max-iter"],
                    parsed["admm-rho"],
                    parsed["admm-eps-primal"],
                    parsed["admm-eps-dual"],
                )
            else
                haskey(base_map, model_key) || error("Invalid model: $model_key")
                base_map[model_key]
            end
        end

        # Optimizer (fallback Ipopt)
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
