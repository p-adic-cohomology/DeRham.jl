@testset "L_polynomial and Tate twists" begin
    # Old LPolynomial(zeta_coeffs) formula, reimplemented here (the function
    # itself was deleted in favor of L_polynomial): T^(deg+1-i) * coeffs[i].
    function old_Lpolynomial(coeffs)
        P, T = polynomial_ring(ZZ, "T")
        deg = length(coeffs) - 1
        return sum(T^(deg + 1 - i) * ZZ(coeffs[i]) for i = 1:(deg+1))
    end

    # Old boat_shape_Lpoly(coeffs, deg, q) formula, reimplemented here (the
    # function itself was deleted in favor of normalized_tate_twist).
    function old_boat_shape_Lpoly(coeffs, deg, q)
        n = deg + 1
        coeffs_new = [ZZ(0) for i = 1:n]
        for i = 1:n
            coeffs_new[i] = div(q * ZZ(coeffs[i]), ZZ(q)^(n - i))
        end
        return coeffs_new
    end

    fermat_cubic_p11 = [-1771561, 0, 43923, 0, -363, 0, 1]
    q_cubic = 11

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

    R, (x, y, z) = polynomial_ring(GF(7), ["x", "y", "z"])
    f_curve = y^2 * z - x^3 - x * z^2 - z^3
    curve_coeffs = DeRham.zeta_coefficients(f_curve)
    q_curve = 7

    examples = [
        ("curve", curve_coeffs, q_curve),
        ("cubic surface", fermat_cubic_p11, q_cubic),
        ("K3", p5_q4K3_sparse_fk_001, q_k3),
    ]

    @testset "L_polynomial matches the old formula" begin
        for (name, coeffs, q) in examples
            @test DeRham.L_polynomial(coeffs) == old_Lpolynomial(coeffs)
        end
    end

    @testset "normalized_tate_twist(c, q, 1) matches the old boat_shape_Lpoly" begin
        for (name, coeffs, q) in examples
            deg = length(coeffs) - 1
            @test DeRham.normalized_tate_twist(coeffs, q, 1) ==
                  old_boat_shape_Lpoly(coeffs, deg, q)
        end
    end

    @testset "round trips" begin
        for (name, coeffs, q) in examples
            raw = ZZRingElem.(coeffs)
            @test DeRham.tate_untwist(DeRham.tate_twist(coeffs, q), q) == raw
            @test DeRham.tate_untwist(
                DeRham.normalized_tate_twist(coeffs, q, 1),
                q;
                s = 1,
            ) == raw
        end
    end

    @testset "normalized_tate_twist throws for non-integral s" begin
        for (name, coeffs, q) in examples
            @test_throws ArgumentError DeRham.normalized_tate_twist(coeffs, q, -1)
        end
    end

    @testset "detection helper maps all three input kinds back to raw" begin
        for (name, coeffs, q) in examples
            raw = ZZRingElem.(coeffs)
            twisted = DeRham.tate_twist(coeffs, q)
            normalized = DeRham.normalized_tate_twist(coeffs, q, 1)

            @test DeRham._normalize_raw_coefficients(coeffs, q) == raw
            @test DeRham._normalize_raw_coefficients(twisted, q) == raw
            @test DeRham._normalize_raw_coefficients(normalized, q) == raw
        end
    end
end
