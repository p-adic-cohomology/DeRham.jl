@testset "AKR algorithm selector support" begin
    _, (x, y, z) = polynomial_ring(GF(5), ["x", "y", "z"])
    hypersurface = DeRham.ProjectiveHypersurface(x^3 + y^3 + z^3)
    for algorithm in (:akr, :akr_naive, :akr_levelbylevel)
        @test DeRham._supports_reduction_algorithm(hypersurface, algorithm)
    end
end

@testset "AKR plane curve full-precision equality" begin
    curve = DeRham.fermat_hypersurface(3, 3, 2)
    _, naive = DeRham.zeta_coefficients(curve; givefrobmat = true, algorithm = :akr_naive)
    @test naive == BigInt[2, 0, 1]
    @test DeRham.zeta_coefficients(curve; algorithm = :akr) == naive
    @test DeRham.zeta_coefficients(curve; algorithm = :akr_levelbylevel) == naive
    @test DeRham.zeta_coefficients(curve) == naive
    @test DeRham.zeta_coefficients(curve; algorithm = :akr_naive, termorder = :lex) == naive
    @test DeRham.zeta_coefficients(
        curve;
        algorithm = :akr_levelbylevel,
        termorder = :lex,
    ) == naive
end

@testset "AKR quartic K3 limited-precision matrices" begin
    _, (x, y, z, w) = polynomial_ring(GF(2), ["x", "y", "z", "w"])
    k3 = x^3*y + y^3*z + z^3*w + w^3*x + x*y*z*w
    default_matrix = DeRham.frobenius_matrix_with_precision(k3, 1; algorithm = :default)
    naive_matrix = DeRham.frobenius_matrix_with_precision(k3, 1; algorithm = :akr_naive)
    levelbylevel_matrix =
        DeRham.frobenius_matrix_with_precision(k3, 1; algorithm = :akr_levelbylevel)
    @test naive_matrix == default_matrix
    @test levelbylevel_matrix == default_matrix
end
