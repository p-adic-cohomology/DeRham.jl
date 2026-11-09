@testset "Projective hypersurface frontend" begin
    Rcurve, (x, y, z) = polynomial_ring(GF(2), ["x", "y", "z"])
    cubic = x^2*y + y^2*z + z^2*x
    cubic_frontend = DeRham.ProjectiveHypersurface(cubic)
    @test DeRham.polynomial(cubic_frontend) === cubic
    @test DeRham.presentation(cubic_frontend) == DeRham.HypersurfacePresentation(2, 3)
    @test DeRham.weight(cubic_frontend) == 1
    @test DeRham.check_input(cubic_frontend) === nothing
    @test DeRham.check_smoothness(cubic)
    @test_throws ArgumentError DeRham.zeta_coefficients(cubic_frontend; series_precision = 3)
    @test DeRham.zeta_coefficients(cubic_frontend) == DeRham.zeta_coefficients(cubic)
    @test typeof(DeRham.frobenius_term_pole_order(cubic_frontend, 1, 0)) == Int

    Rsurface, (a, b, c, d) = polynomial_ring(GF(2), ["a", "b", "c", "d"])
    quartic = a^3*c + a^2*b^2 + a*b*d^2 + a*c^2*d + a*d^3 + b^4 + c^3*d
    quartic_frontend = DeRham.ProjectiveHypersurface(quartic)
    @test DeRham.check_smoothness(quartic)
    @test DeRham.zeta_coefficients(quartic_frontend) == DeRham.zeta_coefficients(quartic)

    for p in (2, 3, 5), N in 1:3, m in 1:3
        @test DeRham.reduction_prescale(cubic_frontend, p * (N + m - 1)) ==
              factorial(ZZ(p * (N + m - 1) - 1))
        for e in 1:(N + m - 1)
            @test DeRham.reduction_prescale(cubic_frontend, p * e) ==
                  factorial(ZZ(p * e - 1))
        end
    end

    for N in 1:4, m in 1:4
        @test DeRham.frobenius_series(cubic_frontend, N, m) == DeRham.computeD(N, m)
        for j in 0:(N - 1)
            @test DeRham.frobenius_term_pole_order(cubic_frontend, m, j) == 2 * (m + j)
        end
    end

    for n in 2:5, i in 0:(n - 1)
        @test DeRham.reduction_prescale(cubic_frontend, n - i) == factorial(ZZ(n - 1 - i))
    end

    dimension = DeRham.presentation(cubic_frontend).n
    p = Int(characteristic(parent(cubic)))
    orders = DeRham.pole_orders(DeRham.presentation(cubic_frontend))
    scale_basis = [[nothing, m] for m in orders]
    N_m = [1, 2]
    expected_scale_exponents = [
        dimension - 1 - DeRham._ibm_vp_factorial(p * (m + N_m[m] - 1) - 1, p) for m in orders
    ]
    @test DeRham._frobenius_column_scale_exponents(cubic_frontend, p, N_m, scale_basis) ==
          expected_scale_exponents
end
