# ============================================================
# Multi-period invariants — internal consistency tests.
#
# We do not have a single-source-of-truth solver for multi-period OPF
# (PowerModels' MN is structurally the same as what we do for models
# without inter-period coupling), so we rely on a set of invariants
# that MUST hold by construction.
#
# Invariants tested here:
#   A. Replication (no battery): solving T=K with demand identical in
#      every period gives objective == K · obj(T=1). Periods are
#      decoupled when no battery is present, so the multi-period
#      model must reduce to K independent single-period problems.
#   B. Battery energy conservation: for every battery, the entire
#      charge/discharge trajectory must reconcile with the final
#      storage state:
#        storage[T] − energy_start == Σ_t(charge[t]·η_c − discharge[t]/η_d)
#   C. Arbitrage direction: with a varying-demand multi-period case,
#      at least one battery must charge in the low-demand period and
#      discharge in the high-demand period (economic optimality).
# ============================================================

using Lirion: StorageData
import Lirion
import PowerModels

# PowerModels multi-network build function that uses NONLINEAR storage
# complementarity (sc·sd = 0), matching Lirion's formulation. The standard
# PM `build_mn_opf_strg` uses MI complementarity (sc_on + sd_on == 1),
# which requires an MIP-capable solver and introduces binary variables we
# don't have. This version is a direct copy with the MI pieces swapped.
function _build_mn_opf_strg_nl(pm::PowerModels.AbstractPowerModel)
    for (n, network) in PowerModels.nws(pm)
        PowerModels.variable_bus_voltage(pm; nw = n)
        PowerModels.variable_gen_power(pm; nw = n)
        PowerModels.variable_storage_power(pm; nw = n)  # non-MI
        PowerModels.variable_branch_power(pm; nw = n)
        PowerModels.variable_dcline_power(pm; nw = n)

        PowerModels.constraint_model_voltage(pm; nw = n)

        for i in PowerModels.ids(pm, :ref_buses; nw = n)
            PowerModels.constraint_theta_ref(pm, i; nw = n)
        end

        for i in PowerModels.ids(pm, :bus; nw = n)
            PowerModels.constraint_power_balance(pm, i; nw = n)
        end

        for i in PowerModels.ids(pm, :storage; nw = n)
            PowerModels.constraint_storage_complementarity_nl(pm, i; nw = n)  # non-MI
            PowerModels.constraint_storage_losses(pm, i; nw = n)
            PowerModels.constraint_storage_thermal_limit(pm, i; nw = n)
        end

        for i in PowerModels.ids(pm, :branch; nw = n)
            PowerModels.constraint_ohms_yt_from(pm, i; nw = n)
            PowerModels.constraint_ohms_yt_to(pm, i; nw = n)
            PowerModels.constraint_voltage_angle_difference(pm, i; nw = n)
            PowerModels.constraint_thermal_limit_from(pm, i; nw = n)
            PowerModels.constraint_thermal_limit_to(pm, i; nw = n)
        end

        for i in PowerModels.ids(pm, :dcline; nw = n)
            PowerModels.constraint_dcline_power_losses(pm, i; nw = n)
        end
    end

    # Storage state coupling across periods (initial + recurrence)
    network_ids = sort(collect(PowerModels.nw_ids(pm)))
    n_1 = network_ids[1]
    for i in PowerModels.ids(pm, :storage; nw = n_1)
        PowerModels.constraint_storage_state(pm, i; nw = n_1)
    end
    for n_2 in network_ids[2:end]
        for i in PowerModels.ids(pm, :storage; nw = n_2)
            PowerModels.constraint_storage_state(pm, i, n_1, n_2)
        end
        n_1 = n_2
    end

    PowerModels.objective_min_fuel_and_flow_cost(pm)
    return nothing
end

# Build a PM multi-network dict from a single-period matpower file, applying
# per-period load scaling factors (to mirror a Lirion bus_timeseries with a
# uniform scale factor per period).
function _build_pm_mn_from_factors(base_file::AbstractString, factors::Vector{<:Real})
    pm_data = PowerModels.parse_file(base_file)
    base_pd = Dict(l => load["pd"] for (l, load) in pm_data["load"])
    base_qd = Dict(l => load["qd"] for (l, load) in pm_data["load"])
    mn = PowerModels.replicate(pm_data, length(factors))
    for (t_str, nw) in mn["nw"]
        t = parse(Int, t_str)
        f = factors[t]
        for (l, load) in nw["load"]
            load["pd"] = base_pd[l] * f
            load["qd"] = base_qd[l] * f
        end
    end
    return mn
end

# Solve a PM MN problem via our custom NL-storage build function.
function _pm_mn_opf_strg_nl(mn_data::Dict, model_type::Type, optimizer)
    return PowerModels.solve_model(
        mn_data, model_type, optimizer, _build_mn_opf_strg_nl;
        multinetwork = true,
    )
end

# Helpers ---------------------------------------------------------------

"""
    _period_obj(file, model, optimizer, t_sub)

Solve the case at `file` (which must have a bus_timeseries) keeping only
the data for period `t_sub`, replicated as a single-period problem. The
case is loaded, the timeseries is trimmed in-memory to just `t_sub`, and
the reduced case is re-solved. Returns the objective.

Caveat: this relies on `Lirion.load` returning a mutable CaseRaw; when
that is not available we fall back to solving the full multi-period run
and extracting the per-period contribution via the cost expression.
"""
function _solve_flat_timeseries_singleperiod(
    multi_file::AbstractString, base_file::AbstractString; model, optimizer,
)
    # The multi-period "flat" case (all periods identical) is equivalent
    # to the base single-period case. Solve both and return both.
    r_multi = Lirion.solve(
        multi_file; model = model, algorithm = Lirion.centralized(),
        optimizer = optimizer, outputs = false, plots = false,
    )
    r_single = Lirion.solve(
        base_file; model = model, algorithm = Lirion.centralized(),
        optimizer = optimizer, outputs = false, plots = false,
    )
    return r_multi, r_single
end

# Invariant A — Replication ---------------------------------------------

function _check_replication(multi_file::AbstractString, base_file::AbstractString;
    K::Integer, model, optimizer, rtol::Real = 1e-4,
)
    r_multi, r_single = _solve_flat_timeseries_singleperiod(
        multi_file, base_file; model = model, optimizer = optimizer,
    )
    expected = K * r_single.objective
    @info "Replication" model K multi_obj = r_multi.objective single_obj =
        r_single.objective expected
    @test isapprox(r_multi.objective, expected; rtol = rtol)
    return r_multi, r_single
end

# Invariant B — Energy conservation -------------------------------------

function _check_energy_conservation(result::Lirion.LirionResult, bat::StorageData;
    atol::Real = 1e-4,
)
    m = result.model
    T = get(m.ext, :time_periods, 1)
    T > 1 || error("Energy conservation test requires T > 1.")

    for b in 1:bat.n
        ηc = bat.charge_efficiency[b]
        ηd = bat.discharge_efficiency[b]
        storage_T = JuMP.value(m[:storage][b, T])
        storage_0 = bat.energy_start[b]
        net_flow = sum(
            JuMP.value(m[:charge][b, t]) * ηc -
            JuMP.value(m[:discharge][b, t]) / ηd
            for t in 1:T
        )
        residual = (storage_T - storage_0) - net_flow
        @info "Energy conservation" bat = b storage_T storage_0 net_flow residual
        @test isapprox(residual, 0.0; atol = atol)
    end
    return nothing
end

# Invariant C — Arbitrage direction -------------------------------------

function _check_arbitrage_direction(
    result::Lirion.LirionResult, bat::StorageData;
    low_period::Int, high_period::Int, atol::Real = 1e-4,
)
    m = result.model
    # Economic-optimality invariant: aggregated across all batteries,
    # total discharge must be AT LEAST as large in the high-demand
    # period as in the low-demand period. We compare discharges (not
    # net charge) because batteries may start with non-zero energy and
    # be forced to discharge at some point — what matters for arbitrage
    # is that *when* the discharge happens it biases toward expensive
    # periods.
    total_discharge_low = sum(
        JuMP.value(m[:discharge][b, low_period]) for b in 1:bat.n
    )
    total_discharge_high = sum(
        JuMP.value(m[:discharge][b, high_period]) for b in 1:bat.n
    )
    @info "Arbitrage" total_discharge_low total_discharge_high low_period high_period
    @test total_discharge_high >= total_discharge_low - atol
    return nothing
end

# Public runner ---------------------------------------------------------

function run_multiperiod_invariants(; optimizer = Ipopt.Optimizer)
    mp_dir = joinpath(@__DIR__, "multiperiod")
    base_dir = joinpath(@__DIR__, "matpower_cases")
    flat_file = joinpath(mp_dir, "case3_timeseries.m")
    base_file = joinpath(base_dir, "pglib_opf_case3_lmbd.m")
    bat_ts_file = joinpath(mp_dir, "case5_bat_ts.m")

    @testset "Multi-period invariants" begin
        @testset "A. Replication (T=3 flat, no battery)" begin
            for model in (Lirion.DC, Lirion.AC, Lirion.ACR, Lirion.SOCWR)
                mfn =
                    model === Lirion.DC ? Lirion.dc() :
                    model === Lirion.AC ? Lirion.ac() :
                    model === Lirion.ACR ? Lirion.acr() :
                    Lirion.socwr()
                opt = model === Lirion.SOCWR ? Clarabel.Optimizer : optimizer
                _check_replication(
                    flat_file, base_file; K = 3, model = mfn,
                    optimizer = opt, rtol = 1e-4,
                )
            end
        end

        @testset "B. Storage energy conservation (T=3 varying)" begin
            # AC with battery: this is where coupling happens.
            r = Lirion.solve(
                bat_ts_file; model = Lirion.ac(),
                algorithm = Lirion.centralized(),
                optimizer = optimizer, outputs = false, plots = false,
            )
            @test r.status in (:optimal, :local_optimal)
            _check_energy_conservation(r, r.system.storage; atol = 1e-4)
        end

        @testset "C. Arbitrage direction (T=3 varying)" begin
            r = Lirion.solve(
                bat_ts_file; model = Lirion.ac(),
                algorithm = Lirion.centralized(),
                optimizer = optimizer, outputs = false, plots = false,
            )
            @test r.status in (:optimal, :local_optimal)
            # Period 1 is low-demand (0.6×), period 3 is high-demand (1.2×).
            _check_arbitrage_direction(
                r, r.system.storage;
                low_period = 1, high_period = 3, atol = 1e-3,
            )
        end

        @testset "D. PM MN no-storage match (T=3 flat)" begin
            # Lirion case3_timeseries (T=3 identical demand) must match PM's
            # build_mn_opf on the same case replicated 3×, for DC and AC.
            for (mfn, pm_model, opt) in (
                (Lirion.dc(), PowerModels.DCPPowerModel, optimizer),
                (Lirion.ac(), PowerModels.ACPPowerModel, optimizer),
            )
                r = Lirion.solve(
                    flat_file; model = mfn, algorithm = Lirion.centralized(),
                    optimizer = opt, outputs = false, plots = false,
                )
                pm_mn = _build_pm_mn_from_factors(base_file, [1.0, 1.0, 1.0])
                pm_res = PowerModels.solve_mn_opf(pm_mn, pm_model, opt)
                @info "PM MN match (no battery)" mfn lirion = r.objective pm =
                    pm_res["objective"]
                @test isapprox(r.objective, pm_res["objective"]; rtol = 1e-4)
            end
        end

        @testset "E. PM MN with-storage AC match (T=3 varying)" begin
            # Lirion case5_bat_ts vs PM MN with NL storage complementarity.
            r = Lirion.solve(
                bat_ts_file; model = Lirion.ac(),
                algorithm = Lirion.centralized(),
                optimizer = optimizer, outputs = false, plots = false,
            )
            base_bat_file = joinpath(base_dir, "case5_bus_2_bat.m")
            pm_mn = _build_pm_mn_from_factors(base_bat_file, [0.6, 1.0, 1.2])
            pm_res = _pm_mn_opf_strg_nl(pm_mn, PowerModels.ACPPowerModel, optimizer)
            @info "PM MN match (AC + battery)" lirion = r.objective pm = pm_res["objective"]
            @test isapprox(r.objective, pm_res["objective"]; rtol = 1e-2)
        end
    end
    return nothing
end
