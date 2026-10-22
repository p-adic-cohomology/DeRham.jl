@testset "divided_frobenius" begin
    @testset "matrix method" begin
        p = 3
        R27, _ = residue_ring(ZZ, p^3)

        F1 = matrix(R27, [3 1; 0 3])
        d1 = DeRham.divided_frobenius(F1, 1)
        @test d1.exponent == 1
        @test d1.p == p
        @test modulus(base_ring(d1.numerator)) == p^2
        @test [lift(ZZ, d1.numerator[r, c]) for r = 1:2, c = 1:2] == [3 1; 0 3]

        F2 = matrix(R27, [0 3; 1 0])
        d2 = DeRham.divided_frobenius(F2, 2)
        @test d2.exponent == 2
        @test d2.p == p
        @test modulus(base_ring(d2.numerator)) == p^3
        @test [lift(ZZ, d2.numerator[r, c]) for r = 1:2, c = 1:2] == [3 0; 0 3]

        R4, _ = residue_ring(ZZ, 2^2)
        @test_throws AssertionError DeRham.divided_frobenius(matrix(R4, [1;;]), 1)

        R9, _ = residue_ring(ZZ, 3^2)
        @test_throws AssertionError DeRham.divided_frobenius(matrix(R9, [1;;]), 2)
        @test_throws AssertionError DeRham.divided_frobenius(matrix(R27, [1;;]), 0)
        @test_throws AssertionError DeRham.divided_frobenius(matrix(R27, [1 2 3; 4 5 6]), 1)
    end

    @testset "quartic method" begin
        p = 3
        _, (x, y, z, w) = polynomial_ring(GF(p), 4)
        f = x^4 + y^4 + z^4 + w^4

        _, (a, b, c, d) = polynomial_ring(GF(5), 4)
        f_nonsmooth = (a^2 + b^2)^2 + c^4 + d^4
        @test DeRham.divided_frobenius(f_nonsmooth, 1) === false
        @test_throws AssertionError DeRham.divided_frobenius(f, 0)

        results = [DeRham.divided_frobenius(f, i) for i = 1:2]
        for (i, result) in enumerate(results)
            @test result.exponent == i
            @test result.p == p
            @test size(result.numerator) == (21, 21)
            @test modulus(base_ring(result.numerator)) == p^(i + 1)
        end

        pole_order_one_column = [lift(ZZ, results[1].numerator[r, 1]) for r = 1:21]
        @test all(mod(entry, p) == 0 for entry in pole_order_one_column)

        pole_order_three_column = [lift(ZZ, results[1].numerator[r, 21]) for r = 1:21]
        @test any(mod(entry, p) != 0 for entry in pole_order_three_column)

        A = 2
        F = DeRham.frobenius_matrix_with_precision(f, [A - 2, A - 1, A])
        Fref = DeRham.frobenius_matrix_with_precision(f, [A + 3, A + 3, A + 3])
        @test size(F) == (21, 21)
        @test size(Fref) == (21, 21)
        @test all(
            mod(lift(ZZ, F[r, c]) - lift(ZZ, Fref[r, c]), p^A) == 0 for r = 1:21, c = 1:21
        )
    end
end

@testset "ogus_nygaard_filtration" begin
    p = 3
    check_congruences = function (F, result)
        F_integer = map_entries(ZZ, lift.(F))
        for (index, basis) in enumerate(result.bases)
            level = index - 1
            image = F_integer^level * basis
            @test all(
                iszero(mod(image[row, column], p^level)) for
                row = 1:nrows(image), column = 1:ncols(image)
            )
        end
    end

    R9, _ = residue_ring(ZZ, p^2)
    F2 = matrix(R9, [0 0; 1 0])
    result2 = DeRham.ogus_nygaard_filtration(F2, 2)
    @test result2.colengths == [0, 1, 1]
    @test result2.bases[end-1] == result2.bases[end]
    check_congruences(F2, result2)
    @test_throws AssertionError DeRham.ogus_nygaard_filtration(F2, 3)

    R27, _ = residue_ring(ZZ, p^3)
    F3 = matrix(R27, [0 0 p^2; p 0 0; 0 1 0])
    result3 = DeRham.ogus_nygaard_filtration(F3, 3)
    @test result3.colengths == [0, 1, 2, 2]
    @test result3.bases[end-1] == result3.bases[end]
    check_congruences(F3, result3)

    R81, _ = residue_ring(ZZ, p^4)
    F4 = matrix(R81, [0 p 0 p; p 0 0 0; 0 1 0 0; 0 0 p^2 0])
    result4 = DeRham.ogus_nygaard_filtration(F4, 4)
    @test result4.colengths == [0, 1, 2, 3, 3]
    @test result4.bases[end-1] == result4.bases[end]
    check_congruences(F4, result4)

    F_bad = matrix(R9, [p 0; 0 p])
    @test_throws ErrorException DeRham.ogus_nygaard_filtration(F_bad, 1)

    _, (x, y, z, w) = polynomial_ring(GF(p), 4)
    f = x^4 + y^4 + z^4 + w^4
    fermat_result = DeRham.ogus_nygaard_filtration(f, 3)
    @test fermat_result.colengths[2] == 1
    @test all(
        fermat_result.colengths[level+1] - fermat_result.colengths[level] in (0, 1) for
        level = 1:(length(fermat_result.colengths)-1)
    )

    high_precision_F = DeRham.frobenius_matrix_with_precision(f, [6, 6, 6])
    check_congruences(high_precision_F, fermat_result)
    @test all(mod(lift(ZZ, high_precision_F[row, 1]), p^2) == 0 for row = 1:21)
    @test all(
        mod(lift(ZZ, high_precision_F[row, column]), p) == 0 for row = 1:21, column = 2:20
    )
    @test any(mod(lift(ZZ, high_precision_F[row, 21]), p) != 0 for row = 1:21)

    for (precision, relative_precision) in ((2, [0, 1, 2]), (3, [1, 2, 3]))
        low_precision_F = DeRham.frobenius_matrix_with_precision(f, relative_precision)
        @test all(
            mod(
                lift(ZZ, low_precision_F[row, column]) -
                lift(ZZ, high_precision_F[row, column]),
                p^precision,
            ) == 0 for row = 1:21, column = 1:21
        )
    end
end
