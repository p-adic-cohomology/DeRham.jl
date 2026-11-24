@testset "square-root interval products" begin
    for p in (101, 103)
        cases = (p - 1, p, p + 1, 3p + 2, 63, 64, 65)
        for modulus in (p, ZZ(p)^25)
            ring, = residue_ring(ZZ, modulus)
            matrixspace = matrix_space(ring, 2, 2)
            a = matrixspace()
            b = matrixspace()
            for j in 1:2, i in 1:2
                a[i, j] = ring(rand(0:(p - 1)))
                b[i, j] = ring(rand(0:(p - 1)))
            end

            big_modulus = modulus > ZZ(2)^64
            g = big_modulus ? ZZRingElem[ZZ(rand(0:(p - 1))) for _ in 1:2] :
                UInt[data(ring(rand(0:(p - 1)))) for _ in 1:2]
            for interval_length in cases
                start = 11
                stop = start + interval_length - 1
                expected = copy(g)
                expected_temp = matrixspace()
                expected_g_temp = copy(g)
                actual = copy(g)
                actual_temp = matrixspace()
                actual_g_temp = copy(g)

                DeRham.finitediff_prodeval_linear!(
                    a,
                    b,
                    start,
                    stop,
                    expected,
                    expected_temp,
                    expected_g_temp,
                )
                DeRham.sqrt_interval_product!(
                    a,
                    b,
                    start,
                    stop,
                    actual,
                    actual_temp,
                    actual_g_temp,
                )
                @test actual == expected
            end
        end
    end
end

@testset "square-root interval product options" begin
    _, (x, y) = polynomial_ring(GF(7), ["x", "y"])
    frontend = DeRham.ProjectiveHypersurface(x^3 + y^3)
    @test_throws ArgumentError zeta_coefficients(
        frontend;
        interval_product_algorithm = :invalid,
    )
    @test_throws ArgumentError zeta_coefficients(frontend; interval_product_min_length = 0)
    @test_throws ArgumentError zeta_coefficients(frontend; interval_product_min_length = true)
end
