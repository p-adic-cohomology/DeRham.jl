@testset "SupportCone algorithm agreement" begin
    fermat_curve = DeRham.fermat_hypersurface(3, 3, 7)
    _, (x, y, z) = polynomial_ring(GF(7), ["x", "y", "z"])
    sparse_quartic = x^4 + y^4 + z^4 + x^2 * y^2

    for curve in (fermat_curve, sparse_quartic)
        supportcone_coefficients = DeRham.zeta_coefficients(
            curve;
            algorithm = :supportcone,
            changef = false,
        )
        depthfirst_coefficients = DeRham.zeta_coefficients(
            curve;
            algorithm = :depthfirst,
            changef = false,
        )
        @test supportcone_coefficients == depthfirst_coefficients
    end

    cubic_surface = DeRham.fermat_hypersurface(4, 3, 7)
    supportcone_matrix = DeRham.frobenius_matrix_with_precision(
        cubic_surface,
        1;
        algorithm = :supportcone,
        S = [0, 1, 2],
        changef = false,
    )
    depthfirst_matrix = DeRham.frobenius_matrix_with_precision(
        cubic_surface,
        1;
        algorithm = :depthfirst,
        S = [0, 1, 2],
        changef = false,
    )
    @test supportcone_matrix !== false
    @test depthfirst_matrix !== false
    @test supportcone_matrix == depthfirst_matrix
end
