function _reference_generalized_binomial(a::Rational{BigInt}, i::Integer)
    value = BigInt(1) // BigInt(1)
    for t = 0:(i-1)
        value *= (a - BigInt(t)) / BigInt(t + 1)
    end
    return value
end

function _reference_computeD_half_integral(N::Integer, k::Integer)
    a = -BigInt(k) // BigInt(2)
    return [
        sum(
            (-1)^(i+j) * _reference_generalized_binomial(a, i) * binomial(i, j) for
            i = j:(N-1);
            init = zero(a),
        ) for j = 0:(N-1)
    ]
end

@testset "Double cover frontend" begin
    for N = 0:6, m = 1:4
        @test DeRham.computeD_half_integral(N, 2m) == QQ.(DeRham.computeD(N, m))
    end

    for k in (1, 3, 5), N = 0:5
        expected = _reference_computeD_half_integral(N, k)
        @test DeRham.computeD_half_integral(N, k) == QQ.(expected)
    end

    Rcurve, (x, z) = polynomial_ring(GF(5), ["x", "z"])
    quartic = x^4 + z^4
    curve = DeRham.DoubleCover(quartic)
    P = DeRham.presentation(curve)
    @test DeRham.polynomial(curve) === quartic
    @test P == DeRham.DoubleCoverPresentation(1, 4)
    @test DeRham.weight(curve) == 1
    @test DeRham.check_input(curve) === nothing
    @test typeof(DeRham.frobenius_term_pole_order(curve, 1, 0)) == Int

    scale = QQ(1)
    for m = 1:6
        @test DeRham.reduction_prescale(curve, m) == scale
        @test DeRham.reduction_prescale(curve, m) isa Oscar.QQFieldElem
        if m < 6
            @test DeRham.reduction_prescale(curve, m + 1) / scale ==
                  inv(DeRham.reduction_scalar(P, m))
            scale *= inv(DeRham.reduction_scalar(P, m))
        end
    end
    @test [DeRham.reduction_prescale(curve, m) for m = 1:3] == [QQ(1), QQ(1, 2), QQ(3, 4)]

    for N = 1:5, m = 1:4
        terminal_scale = DeRham.reduction_prescale(
            curve,
            DeRham.frobenius_term_pole_order(curve, m, N - 1),
        )
        expected = [
            terminal_scale / DeRham.reduction_prescale(
                curve,
                DeRham.frobenius_term_pole_order(curve, m, j),
            ) for j = 0:(N-1)
        ]
        @test DeRham._frobenius_prescale_factors(curve, m, N) == expected
    end

    for N = 1:5, m = 1:4
        @test DeRham.frobenius_series(curve, N, m) ==
              DeRham.computeD_half_integral(N, DeRham.pole_exponent(P, m))
        for j = 0:(N-1)
            @test DeRham.frobenius_term_pole_order(curve, m, j) ==
                  (5 * (2m - 1 + 2j) + 1) ÷ 2
        end
    end

    @test length(DeRham.get_basis_of_cohomology(curve)) == 2
    genus_two = DeRham.DoubleCover(x^6 + z^6)
    @test length(DeRham.get_basis_of_cohomology(genus_two)) == 4

    Rsurface, (a, b, c) = polynomial_ring(GF(5), ["a", "b", "c"])
    sextic = DeRham.DoubleCover(a^6 + b^6 + c^6)
    @test length(DeRham.get_basis_of_cohomology(sextic)) == 21

    odd_branch = x^3 + z^3
    @test_throws ArgumentError DeRham.DoubleCover(odd_branch)
    mismatched_odd_frontend = DeRham.DoubleCover{typeof(odd_branch)}(
        odd_branch,
        DeRham.DoubleCoverPresentation(1, 4),
    )
    @test_throws ArgumentError DeRham.check_input(mismatched_odd_frontend)

    Rchar2, (u, v) = polynomial_ring(GF(2), ["u", "v"])
    char2_frontend = DeRham.DoubleCover(u^4 + v^4)
    @test_throws ArgumentError DeRham.check_input(char2_frontend)
    @test_throws ArgumentError DeRham.get_basis_of_cohomology(char2_frontend)

    Rchar0, (u0, v0) = polynomial_ring(ZZ, ["u0", "v0"])
    char0_frontend = DeRham.DoubleCover(u0^4 + v0^4)
    @test_throws r"finite field" DeRham.check_input(char0_frontend)

    inhomogeneous_branch = x^4 + x^2 + z^4
    @test_throws r"homogeneous" DeRham.check_input(DeRham.DoubleCover(inhomogeneous_branch))

    Rsingular, (s, t) = polynomial_ring(GF(3), ["s", "t"])
    singular_frontend = DeRham.DoubleCover(s^2 * t^2)
    @test DeRham.zeta_coefficients(
        singular_frontend;
        series_precision = 1,
        algorithm_precision = 0,
    ) == false
end

@testset "Double cover zeta suite" begin
    R3, (x3, z3) = polynomial_ring(GF(3), ["x", "z"])
    R7, (x7, z7) = polynomial_ring(GF(7), ["x", "z"])
    suite_cases = (
        ("g1_d4_a", DeRham.DoubleCover(x3^4 - x3*z3^3 + z3^4), 3, [1, 3, 3], 0, (2, 5)),
        (
            "g2_d6_a",
            DeRham.DoubleCover(4*x7^6 - 4*x7*z7^5 + 4*z7^6),
            7,
            [1, 1, 9, 7, 49],
            0,
            (2, 5),
        ),
        (
            "g3_d8_a",
            DeRham.DoubleCover(4*x7^8 - 4*x7*z7^7 + 4*z7^8),
            7,
            [1, 0, 0, 0, 0, 0, 343],
            1,
            (3, 7),
        ),
    )

    for (id, fe, p, expected_ascending, precision_increment, expected_pair) in suite_cases
        basis = DeRham.get_basis_of_cohomology(fe)
        polygon = DeRham._hodge_polygon(basis, DeRham.presentation(fe))
        relative_precision =
            DeRham.calculate_relative_precision(polygon, DeRham.weight(fe), p)
        N = 2 + precision_increment
        M = maximum(relative_precision) + 3 + precision_increment
        @test (N, M) == expected_pair

        for always_use_bigints in (false, true)
            result = DeRham.zeta_coefficients(
                fe;
                algorithm = :depthfirst,
                series_precision = N,
                algorithm_precision = M,
                always_use_bigints,
            )
            next_result = DeRham.zeta_coefficients(
                fe;
                algorithm = :depthfirst,
                series_precision = N + 1,
                algorithm_precision = M + 1,
                always_use_bigints,
            )
            @test result !== false
            @test result !== false && reverse(result) == expected_ascending
            @test result !== false && next_result == result
        end
    end

    R5, (x5, y5, z5) = polynomial_ring(GF(5), ["x", "y", "z"])
    double_sextic = DeRham.DoubleCover(
        x5^6 + 3*x5*y5^5 + x5*y5^2*z5^3 + y5^6 + 2*y5^3*z5^3 + 4*y5^2*z5^4 + z5^6,
    )
    sextic_basis = DeRham.get_basis_of_cohomology(double_sextic)
    sextic_polygon = DeRham._hodge_polygon(sextic_basis, DeRham.presentation(double_sextic))
    sextic_relative_precision =
        DeRham.calculate_relative_precision(sextic_polygon, DeRham.weight(double_sextic), 5)
    sextic_M = maximum(sextic_relative_precision) + 4
    @test length(sextic_basis) == 21
    @test sextic_M == 8
    sextic_result = DeRham.zeta_coefficients(
        double_sextic;
        algorithm = :depthfirst,
        series_precision = 3,
        algorithm_precision = sextic_M,
    )
    next_sextic_result = DeRham.zeta_coefficients(
        double_sextic;
        algorithm = :depthfirst,
        series_precision = 4,
        algorithm_precision = sextic_M + 1,
    )
    expected_sextic_ascending = [
        1,
        0,
        -10,
        75,
        -250,
        -3750,
        6250,
        31250,
        -468750,
        781250,
        9765625,
        -48828125,
        -97656250,
        1464843750,
        -2441406250,
        -12207031250,
        183105468750,
        305175781250,
        -2288818359375,
        7629394531250,
        0,
        -476837158203125,
    ]
    @test sextic_result !== false
    @test sextic_result !== false && reverse(sextic_result) == expected_sextic_ascending
    @test sextic_result !== false && next_sextic_result == sextic_result

    _, (u, v) = polynomial_ring(GF(5), ["u", "v"])
    genus_zero = DeRham.DoubleCover(u^2 + v^2)
    @test isempty(DeRham.get_basis_of_cohomology(genus_zero))
    @test DeRham.zeta_coefficients(
        genus_zero;
        algorithm = :depthfirst,
        series_precision = 1,
        algorithm_precision = 2,
    ) == ZZRingElem[ZZ(1)]
    empty_frobenius, empty_coefficients = DeRham.zeta_coefficients(
        genus_zero;
        algorithm = :depthfirst,
        series_precision = 1,
        algorithm_precision = 2,
        givefrobmat = true,
    )
    @test size(empty_frobenius) == (0, 0)
    @test empty_coefficients == ZZRingElem[ZZ(1)]
end

@testset "Double cover conic Frobenius scaling" begin
    R3, (x, y, z) = polynomial_ring(GF(3), ["x", "y", "z"])
    conics = ((x^2 + y^2 + z^2, -3), (x^2 + y^2 + 2*z^2, 3))

    for (F, expected_entry) in conics
        FM, _ = DeRham.zeta_coefficients(
            DeRham.DoubleCover(F);
            S = [0, 1],
            algorithm = :depthfirst,
            relative_precision = 2,
            series_precision = 1,
            algorithm_precision = 1,
            givefrobmat = true,
        )
        entry = lift(ZZ, FM[1, 1])
        @test mod(entry, 3^2) == mod(expected_entry, 3^2)
    end
end
