# ============================================================
# ======  DummyBus exactness: Ybus equivalence proof  ========
# ============================================================


# Centralized π-branch 2×2 admittance block.
function _ybus_original(r, x, b_c, tap, shift)
    Tc = tap * cis(shift)
    Ys = 1 / (r + im * x)
    Yff = (Ys + im * b_c / 2) / abs2(Tc)
    Yft = -Ys / conj(Tc)
    Ytf = -Ys / Tc
    Ytt = Ys + im * b_c / 2
    return (Yff, Yft, Ytf, Ytt)
end

function _ybus_dummybus_reduced(r, x, b_c, tap, shift)
    Tc = tap * cis(shift)
    Ys = 1 / (r + im * x)
    yh = 2 * Ys                       # each half: 1/((r+jx)/2) = 2·Ys

    # 3-node Ybus over (f=1, D=2, t=3); pure series, NO charging on the halves.
    Yff = yh / abs2(Tc)               # half-A (transformer f→D)
    YfD = -yh / conj(Tc)
    YDf = -yh / Tc
    YDD = yh + yh                     # half-A + half-B meet at D
    YDt = -yh                         # half-B (plain D→t)
    YtD = -yh
    Ytt = yh
    Yft = zero(ComplexF64)
    Ytf = zero(ComplexF64)

    # Relocated line charging as bus shunts (add j·bs on the diagonal).
    Yff += im * (b_c / 2) / abs2(Tc)  # from-end, tap-scaled
    Ytt += im * (b_c / 2)             # to-end

    # Kron-eliminate node D (scalar pivot YDD).
    Yff_r = Yff - YfD * YDf / YDD
    Yft_r = Yft - YfD * YDt / YDD
    Ytf_r = Ytf - YtD * YDf / YDD
    Ytt_r = Ytt - YtD * YDt / YDD
    return (Yff_r, Yft_r, Ytf_r, Ytt_r)
end

function run_dummybus_ybus_equiv_tests(; atol = 1e-12)
    @testset "DummyBus Ybus equivalence (centralized π ≡ DummyBus − D)" begin
        cases = [
            ("plain line  tap=1, b_c≠0", 0.02, 0.06, 0.03, 1.0, 0.0),
            ("transformer tap≠1, shift≠0, b_c≠0", 0.01, 0.20, 0.04, 0.978, deg2rad(3.0)),
            ("transformer tap≠1, shift=0, b_c≠0", 0.01, 0.20, 0.04, 1.05, 0.0),
            ("zero charging b_c=0 (sanity)", 0.02, 0.06, 0.0, 1.0, 0.0),
        ]
        for (name, r, x, b_c, tap, shift) in cases
            Yo = _ybus_original(r, x, b_c, tap, shift)
            Yd = _ybus_dummybus_reduced(r, x, b_c, tap, shift)
            err = maximum(abs.(Yo .- Yd))
            @test err < atol
        end
    end
    return nothing
end
