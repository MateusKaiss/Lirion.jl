# ============================================================
# create_partition — auto-zone assignment for distributed OPF
# ============================================================

function create_partition(::AbstractString; kwargs...)
    return _metis_not_loaded_error()
end

function create_partition(::CaseRaw; kwargs...)
    return _metis_not_loaded_error()
end

function _metis_not_loaded_error()
    return error(
        "create_partition requires Metis.jl, Graphs.jl and SimpleWeightedGraphs.jl.\n" *
        "Install with:  using Pkg; Pkg.add([\"Metis\", \"Graphs\", \"SimpleWeightedGraphs\"])\n" *
        "Then load:     using Metis, Graphs, SimpleWeightedGraphs",
    )
end

# ============================================================
# write_matpower — write a CaseRaw back to .m format
# ============================================================

function write_matpower(case::CaseRaw, path::AbstractString)
    open(path, "w") do io
        name = case.name
        println(io, "function mpc = $name")
        println(io, "mpc.version = '2';")
        println(io, "mpc.baseMVA = $(case.baseMVA);")
        println(io)
        _write_bus_table(io, case.bus)
        _write_gen_table(io, case.gen)
        _write_branch_table(io, case.branch)
        _write_gencost_table(io, case.gencost)
        if case.storage !== nothing
            _write_storage_table(io, case.storage)
        end
    end
    return path
end

@inline function _fmt(x::Float64)
    return isinteger(x) ? string(Int(x)) : string(x)
end

function _write_bus_table(io::IO, b::BusRaw)
    println(io, "%% bus data")
    println(io, "%\tbus_i\ttype\tPd\tQd\tGs\tBs\tarea\tVm\tVa\tbaseKV\tzone\tVmax\tVmin")
    println(io, "mpc.bus = [")
    for i in eachindex(b.bus_i)
        print(io, "\t", b.bus_i[i], "\t", b.bus_type[i], "\t")
        print(io, _fmt(b.pd[i]), "\t", _fmt(b.qd[i]), "\t")
        print(io, _fmt(b.gs[i]), "\t", _fmt(b.bs[i]), "\t")
        print(io, b.area[i], "\t", _fmt(b.vm[i]), "\t", _fmt(b.va[i]), "\t")
        print(io, _fmt(b.base_kv[i]), "\t", b.zone[i], "\t")
        println(io, _fmt(b.vmax[i]), "\t", _fmt(b.vmin[i]), ";")
    end
    println(io, "];")
    return println(io)
end

function _write_gen_table(io::IO, g::GenRaw)
    println(io, "%% generator data")
    println(io, "mpc.gen = [")
    for i in eachindex(g.gen_bus)
        print(io, "\t", g.gen_bus[i], "\t", _fmt(g.pg[i]), "\t", _fmt(g.qg[i]), "\t")
        print(io, _fmt(g.qmax[i]), "\t", _fmt(g.qmin[i]), "\t")
        print(io, _fmt(g.vg[i]), "\t", _fmt(g.mbase[i]), "\t", g.gen_status[i], "\t")
        print(io, _fmt(g.pmax[i]), "\t", _fmt(g.pmin[i]))
        if !isempty(g.pc1) && i <= length(g.pc1)
            print(io, "\t", _fmt(g.pc1[i]), "\t", _fmt(g.pc2[i]))
            print(io, "\t", _fmt(g.qc1min[i]), "\t", _fmt(g.qc1max[i]))
            print(io, "\t", _fmt(g.qc2min[i]), "\t", _fmt(g.qc2max[i]))
            print(io, "\t", _fmt(g.ramp_agc[i]), "\t", _fmt(g.ramp_10[i]))
            print(io, "\t", _fmt(g.ramp_30[i]), "\t", _fmt(g.ramp_q[i]))
            print(io, "\t", _fmt(g.apf[i]))
        end
        println(io, ";")
    end
    println(io, "];")
    return println(io)
end

function _write_branch_table(io::IO, br::BranchRaw)
    println(io, "%% branch data")
    println(io, "mpc.branch = [")
    for i in eachindex(br.f_bus)
        print(io, "\t", br.f_bus[i], "\t", br.t_bus[i], "\t")
        print(io, _fmt(br.br_r[i]), "\t", _fmt(br.br_x[i]), "\t", _fmt(br.br_b[i]), "\t")
        print(
            io,
            _fmt(br.rate_a[i]),
            "\t",
            _fmt(br.rate_b[i]),
            "\t",
            _fmt(br.rate_c[i]),
            "\t",
        )
        tap_v = isnan(br.tap[i]) ? 0.0 : br.tap[i]
        shift_v = isnan(br.shift[i]) ? 0.0 : br.shift[i]
        print(io, _fmt(tap_v), "\t", _fmt(shift_v), "\t")
        print(io, br.br_status[i], "\t")
        angmin_v = isnan(br.angmin[i]) ? -360.0 : br.angmin[i]
        angmax_v = isnan(br.angmax[i]) ? 360.0 : br.angmax[i]
        print(io, _fmt(angmin_v), "\t", _fmt(angmax_v))
        println(io, ";")
    end
    println(io, "];")
    return println(io)
end

function _write_gencost_table(io::IO, gc::GencostRaw)
    isempty(gc.model) && return
    println(io, "%% generator cost data")
    println(io, "mpc.gencost = [")
    for i in eachindex(gc.model)
        print(io, "\t", gc.model[i], "\t", _fmt(gc.startup[i]), "\t")
        print(io, _fmt(gc.shutdown[i]), "\t", gc.ncost[i])
        for c in gc.cost_coefficients[i]
            print(io, "\t", _fmt(c))
        end
        println(io, ";")
    end
    println(io, "];")
    return println(io)
end

function _write_storage_table(io::IO, s)
    println(io, "%% storage data")
    println(io, "mpc.storage = [")
    for i in eachindex(s.storage_bus)
        print(io, "\t", s.storage_bus[i])
        print(io, "\t", _fmt(s.ps[i]), "\t", _fmt(s.qs[i]))
        print(io, "\t", _fmt(s.energy[i]), "\t", _fmt(s.energy_rating[i]))
        print(io, "\t", _fmt(s.charge_rating[i]), "\t", _fmt(s.discharge_rating[i]))
        print(io, "\t", _fmt(s.charge_efficiency[i]), "\t", _fmt(s.discharge_efficiency[i]))
        print(io, "\t", _fmt(s.thermal_rating[i]))
        print(io, "\t", _fmt(s.qmin[i]), "\t", _fmt(s.qmax[i]))
        print(io, "\t", _fmt(s.r[i]), "\t", _fmt(s.x[i]))
        print(io, "\t", _fmt(s.p_loss[i]), "\t", _fmt(s.q_loss[i]))
        print(io, "\t", s.status[i])
        println(io, ";")
    end
    return println(io, "];")
end
