@testset "Degree-one functional equation sign" begin
    @test DeRham.sign_fe(2, 1, 3, [2], [3, 1]) == 1
    @test DeRham.sign_fe(2, 1, 3, [2], [6, 1]) == -1
    @test_throws ArgumentError DeRham.sign_fe(2, 1, 3, [1], [3, 1])

    @test_throws ArgumentError DeRham.sign_fe(2, 1, 2, [2], [2, 1])
    @test DeRham.sign_fe(2, 1, 2, [3], [2, 1]) == 1
    @test DeRham.sign_fe(2, 1, 2, [3], [6, 1]) == -1

    @test_throws ArgumentError DeRham.sign_fe(2, 1, 5, [2], [1, 1])
    @test_throws ArgumentError DeRham.sign_fe(2, 1, 5, [2], [10, 1])
end

@testset "K3 over F3 Frobenius weight" begin
    R, = residue_ring(ZZ, ZZ(27))
    FM = matrix_space(R, 2, 2)()
    polygon = DeRham.SlopesPolygon([1, 1])
    coefficients = DeRham._k3f3_zeta_coefficients_from_frobenius(FM, 3, polygon, [2, 2], 0)
    @test coefficients == ZZRingElem[ZZ(9), ZZ(0), ZZ(1)]
end
