"""
A degree-`2g` raw L-polynomial satisfying the weight-1 functional equation,
built as `prod_i (1 + a_i T + q T^2)` from integer "traces" `as`, so its
roots have absolute value `sqrt(q)` without requiring an actual curve of
genus `g` (no smooth plane curve has genus 2, so the data-form tests below
use such synthetic coefficients for `g = 2, 3`).
"""
function synthetic_weil_coeffs(as::Vector{Int}, q)
    P, T = polynomial_ring(ZZ, "T")
    poly = prod(1 + ZZ(a) * T + ZZ(q) * T^2 for a in as)
    d = degree(poly)
    return [coeff(poly, d - i) for i = 0:d]
end

@testset "jacobian invariants" begin
    R, (x, y, z) = polynomial_ring(GF(7), ["x", "y", "z"])
    f_curve = y^2 * z - x^3 - x * z^2 - z^3
    q_curve = 7
    curve_coeffs = DeRham.zeta_coefficients(f_curve)

    coeffs_g2 = synthetic_weil_coeffs([3, -5], 11)
    q_g2 = 11
    coeffs_g3 = synthetic_weil_coeffs([1, -2, 3], 5)
    q_g3 = 5

    @testset "jacobian_zeta_coefficients(f) matches zeta_coefficients(f)" begin
        @test DeRham.jacobian_zeta_coefficients(f_curve) == curve_coeffs
    end

    @testset "genus 1: jacobian_point_counts matches point_counts(coeffs, q, 1, k) (J is E)" begin
        @test DeRham.jacobian_point_counts(curve_coeffs, q_curve, 3) ==
              DeRham.point_counts(curve_coeffs, q_curve, 1, 3)
        @test DeRham.jacobian_point_counts(f_curve, 3) ==
              DeRham.jacobian_point_counts(curve_coeffs, q_curve, 3)
    end

    @testset "genus 2 and 3: #J(F_q) matches L_polynomial(coeffs) evaluated at 1" begin
        P, T = polynomial_ring(ZZ, "T")
        for (coeffs, q) in [(coeffs_g2, q_g2), (coeffs_g3, q_g3)]
            L1 = evaluate(DeRham.L_polynomial(coeffs; ring = P), ZZ(1))
            @test DeRham.jacobian_point_counts(coeffs, q, 1)[1] == L1
        end
    end

    @testset "twisted-input round trip (point counts)" begin
        for (coeffs, q, g) in
            [(curve_coeffs, q_curve, 1), (coeffs_g2, q_g2, 2), (coeffs_g3, q_g3, 3)]
            tw = DeRham.tate_twist(coeffs, q)
            normtw = DeRham.normalized_tate_twist(coeffs, q, g)

            pc_raw = DeRham.jacobian_point_counts(coeffs, q, 2)
            @test DeRham.jacobian_point_counts(tw, q, 2) == pc_raw
            @test DeRham.jacobian_point_counts(normtw, q, 2) == pc_raw
        end
    end

    @testset "genus 1: jacobian_zeta_function matches zeta_function(coeffs, q, 1)" begin
        @test DeRham.jacobian_zeta_function(curve_coeffs, q_curve) ==
              DeRham.zeta_function(curve_coeffs, q_curve, 1)
        @test DeRham.jacobian_zeta_function(f_curve) ==
              DeRham.jacobian_zeta_function(curve_coeffs, q_curve)
    end

    @testset "deg P_i == binomial(2g, i), and P_0 = 1 - T" begin
        P, T = polynomial_ring(ZZ, "T")
        for (coeffs, q, g) in
            [(curve_coeffs, q_curve, 1), (coeffs_g2, q_g2, 2), (coeffs_g3, q_g3, 3)]
            P_list = DeRham._jacobian_wedge_charpolys(coeffs, q; ring = P)
            @test [degree(Pi) for Pi in P_list] == [binomial(2g, i) for i = 0:2g]
            @test P_list[1] == 1 - T
        end
    end

    @testset "zeta function series coefficients match jacobian_point_counts (genus 2)" begin
        Z = DeRham.jacobian_zeta_function(coeffs_g2, q_g2)
        num = numerator(Z)
        den = denominator(Z)
        R, Tser = power_series_ring(QQ, 13, "T"; model = :capped_absolute)
        numser = sum(coeff(num, i) * Tser^i for i = 0:degree(num))
        denser = sum(coeff(den, i) * Tser^i for i = 0:degree(den))
        logz = log(numser * inv(denser))
        series_counts = [ZZ(coeff(logz, i) * i) for i = 1:3]
        @test series_counts == DeRham.jacobian_point_counts(coeffs_g2, q_g2, 3)
    end

    @testset "twisted-input round trip (zeta function)" begin
        for (coeffs, q, g) in
            [(curve_coeffs, q_curve, 1), (coeffs_g2, q_g2, 2), (coeffs_g3, q_g3, 3)]
            tw = DeRham.tate_twist(coeffs, q)
            normtw = DeRham.normalized_tate_twist(coeffs, q, g)

            Z_raw = DeRham.jacobian_zeta_function(coeffs, q)
            @test DeRham.jacobian_zeta_function(tw, q) == Z_raw
            @test DeRham.jacobian_zeta_function(normtw, q) == Z_raw
        end
    end

    @testset "ArgumentError on a surface" begin
        R5, (x1, x2, x3, x4) = polynomial_ring(GF(5), ["x1", "x2", "x3", "x4"])
        fk3 = x1^4 + x2^4 + x3^4 + x4^4
        @test_throws ArgumentError DeRham.jacobian_zeta_coefficients(fk3)
        @test_throws ArgumentError DeRham.jacobian_point_counts(fk3, 2)
        @test_throws ArgumentError DeRham.jacobian_zeta_function(fk3)
    end

    @testset "non-smooth returns false" begin
        f_singular = y^2 * z - x^3
        @test DeRham.jacobian_point_counts(f_singular, 2) == false
        @test DeRham.jacobian_zeta_function(f_singular) == false
    end
end
