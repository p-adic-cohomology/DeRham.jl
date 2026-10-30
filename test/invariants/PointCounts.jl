"""
Expands the zeta function `Z` (an element of `fraction_field(ZZ[T])`) as a
power series and reads off the first `k` point counts via the
log-derivative identity `log(Z) = sum_i #X(F_{q^i}) T^i / i`. Used only to
independently check `point_counts` against `zeta_function`.
"""
function point_counts_via_zeta_log(Z, k)
    num = numerator(Z)
    den = denominator(Z)
    R, Tser = power_series_ring(QQ, k + 10, "T"; model = :capped_absolute)
    numser = sum(coeff(num, i) * Tser^i for i = 0:degree(num))
    denser = sum(coeff(den, i) * Tser^i for i = 0:degree(den))
    logz = log(numser * inv(denser))
    return [ZZ(coeff(logz, i) * i) for i = 1:k]
end

@testset "point_counts and zeta_function" begin
    q_curve = 7
    m_curve = 1

    fermat_cubic_p11 = [-1771561, 0, 43923, 0, -363, 0, 1]
    q_cubic = 11
    m_cubic = 2

    p5_q4K3_sparse_fk_001 = [
        -476837158203125,
        95367431640625,
        15258789062500,
        -6103515625000,
        305175781250,
        61035156250,
        36621093750,
        -7324218750,
        -1220703125,
        439453125,
        0,
        0,
        -3515625,
        390625,
        93750,
        -18750,
        -1250,
        -250,
        200,
        -20,
        -5,
        1,
    ]
    q_k3 = 5
    m_k3 = 2

    curve_coeffs = [7, 1, 1]

    examples = [
        ("curve", curve_coeffs, q_curve, m_curve),
        ("cubic surface", fermat_cubic_p11, q_cubic, m_cubic),
        ("K3", p5_q4K3_sparse_fk_001, q_k3, m_k3),
    ]

    @testset "zeta_function log-derivative identity matches point_counts" begin
        for (name, coeffs, q, m) in examples
            Z = DeRham.zeta_function(coeffs, q, m)
            @test point_counts_via_zeta_log(Z, 3) == DeRham.point_counts(coeffs, q, m, 3)
        end
        @test point_counts_via_zeta_log(
            DeRham.zeta_function(curve_coeffs, q_curve, m_curve),
            3,
        ) == DeRham.point_counts(curve_coeffs, q_curve, m_curve, 3)
    end

    @testset "twisted-input round trip" begin
        for (name, coeffs, q, m) in examples
            tw = DeRham.tate_twist(coeffs, q)
            normtw = DeRham.normalized_tate_twist(coeffs, q, 1)

            pc_raw = DeRham.point_counts(coeffs, q, m, 2)
            @test DeRham.point_counts(tw, q, m, 2) == pc_raw
            @test DeRham.point_counts(normtw, q, m, 2) == pc_raw

            Z_raw = DeRham.zeta_function(coeffs, q, m)
            @test DeRham.zeta_function(tw, q, m) == Z_raw
            @test DeRham.zeta_function(normtw, q, m) == Z_raw
        end
    end

    @testset "non-smooth returns false" begin
        R5, (y1, y2, y3, y4) = polynomial_ring(GF(5), ["x1", "x2", "x3", "x4"])
        f_nonsmooth = (y1^2 + y2^2)^2 + y3^4 + y4^4
        @test DeRham.point_counts(f_nonsmooth, 2) == false
        @test DeRham.zeta_function(f_nonsmooth) == false
    end
end
