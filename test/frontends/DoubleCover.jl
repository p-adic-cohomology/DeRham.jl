function _reference_generalized_binomial(a::Rational{BigInt}, i::Integer)
    value = BigInt(1) // BigInt(1)
    for t in 0:(i-1)
        value *= (a - BigInt(t)) / BigInt(t + 1)
    end
    return value
end

function _reference_computeD_half_integral(N::Integer, k::Integer)
    a = -BigInt(k) // BigInt(2)
    return [
        sum(
            (-1)^(i+j) * _reference_generalized_binomial(a, i) * binomial(i, j) for i = j:(N-1);
            init = zero(a),
        ) for j = 0:(N-1)
    ]
end

@testset "Double cover frontend" begin
    for N in 0:6, m in 1:4
        @test DeRham.computeD_half_integral(N, 2m) == QQ.(DeRham.computeD(N, m))
    end

    for k in (1, 3, 5), N in 0:5
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

    scale = QQ(1)
    for m in 1:6
        @test DeRham.reduction_prescale(curve, m) == scale
        @test DeRham.reduction_prescale(curve, m) isa Oscar.QQFieldElem
        if m < 6
            @test DeRham.reduction_prescale(curve, m + 1) / scale ==
                  inv(DeRham.reduction_scalar(P, m))
            scale *= inv(DeRham.reduction_scalar(P, m))
        end
    end
    @test [DeRham.reduction_prescale(curve, m) for m in 1:3] == [QQ(1), QQ(1, 2), QQ(3, 4)]

    for N in 1:5, m in 1:4
        @test DeRham.frobenius_series(curve, N, m) ==
              DeRham.computeD_half_integral(N, DeRham.pole_exponent(P, m))
        for j in 0:(N-1)
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

    Rsingular, (s, t) = polynomial_ring(GF(3), ["s", "t"])
    singular_frontend = DeRham.DoubleCover(s^2 * t^2)
    @test DeRham.zeta_coefficients(
        singular_frontend;
        series_precision = 1,
        algorithm_precision = 0,
    ) == false
end
