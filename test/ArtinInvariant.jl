@testset "artin_invariant" begin
    @testset "polynomial method" begin
        p = 3
        _, (x, y, z, w) = polynomial_ring(GF(p), ["x", "y", "z", "w"])

        sigma_two = x^4 + y^3 * z + y * z^3 + x * y^2 * w + w^4
        @test DeRham.artin_invariant(sigma_two) == 2

        sigma_ten =
            x^4 +
            x^3 * y +
            x * y^3 +
            x^2 * y * z +
            y^3 * z +
            x^3 * w +
            y^2 * z * w +
            z^3 * w +
            x * y * w^2 +
            x * w^3
        @test DeRham.artin_invariant(sigma_ten) == 10
    end

    @testset "matrix method" begin
        p = 3
        _, (x, y, z, w) = polynomial_ring(GF(p), ["x", "y", "z", "w"])
        fermat = x^4 + y^4 + z^4 + w^4
        F, column_moduli = DeRham.frobenius_matrix_with_precision(
            fermat,
            [0, 1, 2];
            return_column_moduli = true,
        )
        @test F !== false
        @test F isa MatElem
        @test size(F) == (21, 21)
        expected_column_moduli = vcat(BigInt(p)^6, fill(BigInt(p)^3, 19), BigInt(p)^2)
        @test length(column_moduli) == size(F, 2)
        @test column_moduli == expected_column_moduli
        @test DeRham.artin_invariant(F; column_moduli) == 1
        @test_throws UndefKeywordError DeRham.artin_invariant(F)
        F_array = [F[row, column] for row in axes(F, 1), column in axes(F, 2)]
        @test DeRham.artin_invariant(F_array; column_moduli) == 1

        low_precision_moduli = copy(column_moduli)
        low_precision_moduli[2] = p
        @test_throws AssertionError DeRham.artin_invariant(
            F;
            column_moduli = low_precision_moduli,
        )

        @test_throws AssertionError DeRham.artin_invariant(
            matrix(parent(F[1, 1]), zeros(Int, 20, 21));
            column_moduli = fill(BigInt(p)^2, 21),
        )

        F_bad_first_column = deepcopy(F)
        F_bad_first_column[1, 1] = one(parent(F_bad_first_column[1, 1]))
        @test_throws AssertionError DeRham.artin_invariant(
            F_bad_first_column;
            column_moduli,
        )

        R3, _ = residue_ring(ZZ, p)
        @test_throws AssertionError DeRham.artin_invariant(
            matrix(R3, zeros(Int, 21, 21));
            column_moduli = fill(BigInt(p), 21),
        )
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
