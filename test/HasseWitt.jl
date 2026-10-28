@testset "hasse_witt_matrix and a_number" begin
    R, (x, y, z) = polynomial_ring(GF(7), ["x", "y", "z"])
    f = y^2 * z - x^3 - x * z^2 - z^3

    @testset "f form matches the (FM, h) data form" begin
        n = 2
        r_m = fill(0, n)
        r_m[1] = 1
        basis = DeRham.get_basis_of_cohomology(f)
        h = DeRham._hodge_polygon(basis, n).slopelengths[1]
        FM, column_moduli = DeRham.frobenius_matrix_with_precision(
            f,
            r_m;
            basis = basis,
            return_column_moduli = true,
        )
        @test modulus(base_ring(FM)) == BigInt(7)^2
        FM_array = [FM[row, column] for row in axes(FM, 1), column in axes(FM, 2)]

        @test DeRham.hasse_witt_matrix(FM, h; column_moduli) == DeRham.hasse_witt_matrix(f)
        @test DeRham.a_number(FM, h; column_moduli) == DeRham.a_number(f)
        @test DeRham.hasse_witt_matrix(FM_array, h; column_moduli) ==
              DeRham.hasse_witt_matrix(FM, h; column_moduli)
        @test DeRham.a_number(FM_array, h; column_moduli) ==
              DeRham.a_number(FM, h; column_moduli)
    end

    @testset "pinned regression values (the curve is ordinary: a_7 = 3, coprime to 7)" begin
        @test DeRham.hasse_witt_matrix(f) == matrix(GF(7), [3;;])
        @test DeRham.a_number(f) == 0
    end

    @testset "a_number only supports curves" begin
        R5, (a, b, c, d) = polynomial_ring(GF(5), ["a", "b", "c", "d"])
        fk3 = a^4 + b^4 + c^4 + d^4
        @test_throws ArgumentError DeRham.a_number(fk3)
    end
end
