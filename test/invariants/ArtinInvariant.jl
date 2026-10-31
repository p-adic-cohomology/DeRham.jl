@testset "artin_invariant" begin
    @testset "matrix method" begin
        p = 3
        _, (x, y, z, w) = polynomial_ring(GF(p), ["x", "y", "z", "w"])
        fermat = x^4 + y^4 + z^4 + w^4
        F = DeRham.frobenius_matrix_with_precision(fermat, [0, 1, 2])

        @test F !== false
        @test size(F) == (21, 21)
        @test DeRham.artin_invariant(F) == 1

        @test_throws AssertionError DeRham.artin_invariant(
            matrix(base_ring(F), zeros(Int, 20, 21)),
        )

        F_bad_first_column = deepcopy(F)
        F_bad_first_column[1, 1] = one(base_ring(F_bad_first_column))
        @test_throws AssertionError DeRham.artin_invariant(F_bad_first_column)

        R3, _ = residue_ring(ZZ, p)
        @test_throws AssertionError DeRham.artin_invariant(matrix(R3, zeros(Int, 21, 21)))
    end

    @testset "input guards" begin
        _, (a, b, c, d) = polynomial_ring(GF(2), ["a", "b", "c", "d"])
        @test_throws AssertionError DeRham.artin_invariant(a^4 + b^4 + c^4 + d^4)

        _, (u, v, r, s) = polynomial_ring(GF(9), ["u", "v", "r", "s"])
        @test_throws AssertionError DeRham.artin_invariant(u^4 + v^4 + r^4 + s^4)

        _, (u, v, r) = polynomial_ring(GF(3), ["u", "v", "r"])
        @test_throws AssertionError DeRham.artin_invariant(u^4 + v^4 + r^4)

        _, (x, y, z, w) = polynomial_ring(GF(3), ["x", "y", "z", "w"])
        @test_throws AssertionError DeRham.artin_invariant(x^3 + y^3 + z^3 + w^3)
        @test_throws AssertionError DeRham.artin_invariant(x^4 + y^4 + z^4 + w^4 + x)

        _, (a, b, c, d) = polynomial_ring(GF(5), ["a", "b", "c", "d"])
        nonsmooth = (a^2 + b^2)^2 + c^4 + d^4
        @test DeRham.artin_invariant(nonsmooth) === false
    end
end
