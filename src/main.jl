# ============================================================
# src/main.jl
# ============================================================
@inline function _input_basedir(path::AbstractString)::String
    p = abspath(String(path))
    return isdir(p) ? p : dirname(p)  # file => folder; folder => itself
end

@inline function _plots_outdir(input_path::AbstractString)::String
    return joinpath(_input_basedir(input_path), "plots_outputs")
end

@inline function _plots_stem(model_type)::String
    if model_type isa DistributedOPFModel
        return "admm_" * string(model_type.base_model) * "_" * string(model_type.method)
    else
        # DC/AC/SOCP
        return "opf_" * lowercase(string(model_type))
    end
end

@inline function _maybe_build_plots!(
    write_plots::Bool,
    input_path::AbstractString,
    model_type,
    out::LirionRunResult,
)
    write_plots || return nothing

    diag = out.diagnostics
    if diag isa ADMMHistory
        outdir = _plots_outdir(input_path)
        stem = _plots_stem(model_type)
        mkpath(outdir)
        build_plots(diag; outdir = outdir, stem = stem)
    else
        @warn "--plots foi solicitado, mas não há ADMMHistory neste run (nenhum plot gerado)."
    end
    return nothing
end

function run_lirion(inputs::LirionInputs)::LirionRunResult
    @info "Case: $(inputs.database_path)"
    @info "Model: $(inputs.model_type)"
    @info "Optimizer (raw): $(inputs.optimizer)"

    if is_distributed_model(inputs.model_type)
        return run_distributed(inputs)
    else
        result_model = solve_model(inputs)
        sys = inputs.system
        data = sys.raw

        result_dict = build_result_dict(result_model, data)

        if inputs.write_outputs
            save_result_to_csv(result_dict, data, inputs.database_path)
        end

        return LirionRunResult(;
            result = result_dict,
            model = result_model,
            system = sys,
            diagnostics = NoDiagnostics(),
        )
    end
end

@inline function _legacy_return(out::LirionRunResult)
    if out.result isa AbstractDict
        return out.result, out.model, out.system
    end
    return out.result
end

function main(ARGS)
    args = Args(ARGS)
    inputs = build_inputs(args)

    out = run_lirion(inputs)

    _maybe_build_plots!(args.write_plots, args.database_path, args.model_type, out)

    return _legacy_return(out)
end

function julia_main()::Cint
    main(ARGS)
    return 0
end
