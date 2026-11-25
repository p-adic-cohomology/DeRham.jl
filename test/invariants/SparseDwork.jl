@testset "sparse Dwork domains and exact reconstruction" begin
    R, (x, y, z) = polynomial_ring(GF(7), [:x, :y, :z])
    # Direct projective counts over F7, F49 and F343 give traces -2,-38,34.
    f = 2*x^4 + 3*y^4 + 4*z^4 + x^2*y^2
    P = ZZ.([343, 98, 147, 28, 21, 2, 1])
    @test DeRham.sparse_dwork_coefficients(DeRham.ProjectiveHypersurface(f)) == P
    @test DeRham.sparse_dwork_zeta(3*f) == DeRham.zeta_function(P, 7, 1)
    @test_throws ArgumentError DeRham.sparse_dwork_coefficients(x^3+y^3+z^3-3*x*y*z)
    @test_throws ArgumentError DeRham.sparse_dwork_coefficients(x^4+y^4+x*y*z^2)
    @test_throws ArgumentError DeRham.sparse_dwork_coefficients(f+x*y*z^2)
    @test_throws ArgumentError DeRham.sparse_dwork_coefficients(f+1)
    @test_throws ArgumentError DeRham.sparse_dwork_coefficients(zero(R))
    for K in (GF(3), GF(7, 2), QQ)
        _, (a, b, c) = polynomial_ring(K, [:a, :b, :c])
        @test_throws ArgumentError DeRham.sparse_dwork_coefficients(a^3+b^3+c^3)
    end
    @test isnothing(DeRham._sparse_dwork_lift([BigInt[0, 0], BigInt[0, 0]], 7, 1, 1))
end
