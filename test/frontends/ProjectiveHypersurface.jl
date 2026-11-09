@testset "Projective hypersurface frontend" begin
    cubic = DeRham.fermat_hypersurface(3, 3, 7)
    cubic_frontend = DeRham.ProjectiveHypersurface(cubic)
    @test DeRham.polynomial(cubic_frontend) === cubic
    @test DeRham.presentation(cubic_frontend) == DeRham.HypersurfacePresentation(2, 3)
    @test DeRham.weight(cubic_frontend) == 1
    @test DeRham.check_input(cubic_frontend) === nothing
    @test DeRham.zeta_coefficients(cubic_frontend) == DeRham.zeta_coefficients(cubic)

    quartic = DeRham.fermat_hypersurface(4, 4, 3)
    quartic_frontend = DeRham.ProjectiveHypersurface(quartic)
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
            @test DeRham.frobenius_term_pole_order(cubic_frontend, m, j) == 7 * (m + j)
        end
    end

    for n in 2:5, i in 0:(n - 1)
        @test DeRham.reduction_prescale(cubic_frontend, n - i) == factorial(ZZ(n - 1 - i))
    end
end
