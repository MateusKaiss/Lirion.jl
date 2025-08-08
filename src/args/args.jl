function finish_path(path::String)
    return abspath(path)
end

function parse_commandline(args)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "path"
        help = "Path to the case (.m or .csv folder)"
        arg_type = String
        required = true

        "--model"
        help = "Model type to use: DC, AC, or SOCP"
        arg_type = String
        required = true

        "--optimizer"
        help = "Optimizer function name (e.g., Ipopt.Optimizer)"
        arg_type = String
        required = true

        "--plot-outputs"
        help = "Generate plots"
        action = :store_true

        "--write-lp"
        help = "Write LP file"
        action = :store_true

        "--outputs"
        help = "Generate CSV output files"
        action = :store_true
    end

    parsed_args = parse_args(args, s)
    parsed_args["path"] = finish_path(parsed_args["path"])

    if !(isfile(parsed_args["path"]) || isdir(parsed_args["path"]))
        error("Invalid path: $(parsed_args["path"]) is not a file or directory.")
    end

    return parsed_args
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
            help = "Path to MATPOWER .m or CSV folder"

            "--model"
            help = "Model to use (dc, ac, socp)"
            arg_type = String
            required = true

            "--optimizer"
            help = "Optimizer (e.g., Ipopt.Optimizer, Gurobi.Optimizer)"
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
        end

        parsed = parse_args(raw_args, s)

        model_map = Dict(
            "dc" => DC,
            "ac" => AC,
            "socp" => SOCP,
            "distributed-admm-dc" =>
                DistributedOPFModel(DC, :admm, 350 , 1000, 1e-7, 1e-7),
                # DistributedOPFModel(DC, :admm, 350 , 500_000, 1e-7, 1e-7),
                # DistributedOPFModel(DC, :admm, 350 , 1_000_000, 1e-7, 1e-7),
            "distributed-admm-ac" =>
                DistributedOPFModel(AC, :admm, 100, 10.0, 1e-4, 1e-4),
            "distributed-admm-socp" =>
                DistributedOPFModel(SOCP, :admm, 100, 10.0, 1e-4, 1e-4),
        )

        model_key = lowercase(strip(parsed["model"]))
        model_type = get(model_map, model_key) do
            return error("Invalid model: $model_key")
        end

        optimizer = nothing
        if haskey(parsed, "optimizer")
            try
                optimizer = eval(Meta.parse(parsed["optimizer"]))
            catch
                @warn "Invalid optimizer: $(parsed["optimizer"]). Defaulting to Ipopt."
                optimizer = Ipopt.Optimizer
            end
        else
            optimizer = Ipopt.Optimizer
        end

        return new(
            abspath(parsed["path"]),
            model_type,
            optimizer,
            get(parsed, "plot-outputs", false),
            get(parsed, "write-lp", false),
            get(parsed, "outputs", false),
        )
    end
end
