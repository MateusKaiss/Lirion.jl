# ============================================================
# Fast regression tests against frozen fixtures.
#
# Each TOML in `test/fixtures/centralized/` pins the Lirion objective
# for every (case, model) pair the regeneration script validated
# against PowerModels. These tests just re-solve with Lirion and assert
# the objective matches the fixture within a tight tolerance. They
# replace `test_objective.jl`'s PM-roundtrip as the default fast suite
# — no PowerModels, no Clarabel noise.
#
# Regenerate fixtures with:
#   julia --project=. test/scripts/regen_fixtures.jl
# ============================================================

using TOML

const FIXTURES_CENTRALIZED = joinpath(@__DIR__, "fixtures", "centralized")
const CASES_CENTRALIZED = joinpath(@__DIR__, "matpower_cases")

const _MODEL_CTOR = Dict{String, Any}(
    "dc" => Lirion.dc,
    "ac" => Lirion.ac,
    "acr" => Lirion.acr,
    "ivr" => Lirion.ivr,
    "lpac" => Lirion.lpac,
    "socwr" => Lirion.socwr,
    "socwrnl" => Lirion.socwrnl,
    "qc" => Lirion.qc,
    "sdpwrm" => Lirion.sdpwrm,
    "socbf" => Lirion.socbf,
    "bfa" => Lirion.bfa,
    "nfa" => Lirion.nfa,
)

function _optimizer_for_tag(tag::AbstractString)
    return tag in ("socwr", "socbf", "bfa", "sdpwrm") ?
           optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false) :
           optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
end

# Per-model tolerance. SDPWRM is numerically fragile with Clarabel
# (often terminates NUMERICAL_ERROR / ALMOST_OPTIMAL) so we accept a
# 2% drift from the stored objective to absorb run-to-run noise.
@inline _rtol_for_tag(tag::AbstractString, default::Real) =
    tag == "sdpwrm" ? max(default, 2e-2) : default

function run_centralized_fixture_tests(; rtol::Real = 1e-3)
    fixture_files = sort!(filter(endswith(".toml"), readdir(FIXTURES_CENTRALIZED)))

    @testset "Centralized vs fixtures" begin
        for fname in fixture_files
            fx = TOML.parsefile(joinpath(FIXTURES_CENTRALIZED, fname))
            case = fx["metadata"]["case"]
            file = joinpath(CASES_CENTRALIZED, case)
            isfile(file) || continue
            objectives = fx["objectives"]::Dict

            @testset "$case" begin
                for (tag, expected) in sort(collect(objectives))
                    haskey(_MODEL_CTOR, tag) || continue
                    ctor = _MODEL_CTOR[tag]
                    opt = _optimizer_for_tag(tag)

                    got = try
                        r = Lirion.solve(
                            file; model = ctor(),
                            algorithm = Lirion.centralized(),
                            optimizer = opt, outputs = false, plots = false,
                        )
                        Float64(r.objective)
                    catch e
                        @warn "Lirion.solve threw" case tag exc = typeof(e)
                        @test false   # visible failure
                        continue
                    end
                    tol = _rtol_for_tag(tag, rtol)
                    @test isapprox(got, Float64(expected); rtol = tol)
                end
            end
        end
    end
    return nothing
end
