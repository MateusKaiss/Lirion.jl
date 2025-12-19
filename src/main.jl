# ============================================================
# src/main.jl
# ============================================================

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

        return LirionRunResult(
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
    diag = out.diagnostics
    if diag isa ADMMHistory
        stem = "admm_" * string(typeof(inputs.model_type))
        try
            if inputs.model_type isa DistributedOPFModel
                stem = "admm_" * string(inputs.model_type.base_model) * "_" * string(inputs.model_type.method)
            end
        catch
        end
        build_plots(diag; outdir = ".", stem = stem)
    end
    return _legacy_return(out)
end

function julia_main()::Cint
    main(ARGS)
    return 0
end