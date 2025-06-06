using ArgParse

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
    model_type::Type{<:OPFModel}
    optimizer
    plot_outputs::Bool
    write_lp::Bool

    function Args(raw_args::Vector{String})
        parsed_args = parse_commandline(raw_args)

        model = lowercase(strip(parsed_args["model"]))
        model_type = model == "dc"   ? PowerOptimization.DC :
                     model == "ac"   ? PowerOptimization.AC :
                     model == "socp" ? PowerOptimization.SOCP :
                     error("Invalid model: $model")

        optimizer_name = strip(parsed_args["optimizer"])
        optimizer_expr = Meta.parse(optimizer_name)

        optimizer = try
            Base.eval(Main, optimizer_expr)
        catch
            error("Failed to evaluate optimizer: $optimizer_name. Make sure it's loaded in Main.")
        end

        return new(
            parsed_args["path"],
            model_type,
            optimizer,
            get(parsed_args, "plot-outputs", false),
            get(parsed_args, "write-lp", false),
        )
    end
end