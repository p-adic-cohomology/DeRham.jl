function DeRham.my_mul!(output::Matrix{Int}, input::Matrix{Int}, scalar::Int)
    for j in axes(output, 2), i in axes(output, 1)
        output[i, j] = input[i, j] * scalar
    end
    return output
end

function DeRham.my_add!(
    output::Matrix{Int},
    left::Matrix{Int},
    right::Matrix{Int},
)
    for j in axes(output, 2), i in axes(output, 1)
        output[i, j] = left[i, j] + right[i, j]
    end
    return output
end

function DeRham.my_matvecmul!(
    output::Vector{Int},
    matrix::Matrix{Int},
    vector::Vector{Int},
)
    for i in axes(matrix, 1)
        output[i] = 0
        for j in axes(matrix, 2)
            output[i] += matrix[i, j] * vector[j]
        end
    end
    return output
end

@testset "square-root interval products" begin
    for p in (101, 103)
        cases = (1, 2, 3, p - 1, p, p + 1, 3p + 2, 63, 64, 65)
        for modulus in (p, ZZ(p)^25)
            ring, = residue_ring(ZZ, modulus)
            matrixspace = matrix_space(ring, 3, 3)
            a = matrixspace()
            b = matrixspace()
            a_values = [1 1 0; 0 1 1; 1 0 1]
            b_values = [1 0 0; 0 2 0; 0 0 3]
            for j in 1:3, i in 1:3
                a[i, j] = ring(a_values[i, j])
                b[i, j] = ring(b_values[i, j])
            end

            big_modulus = modulus > ZZ(2)^64
            g = big_modulus ? ZZRingElem[ZZ(i + 1) for i in 1:3] :
                UInt[data(ring(i + 1)) for i in 1:3]
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

@testset "big-modulus interval points are reduced" begin
    ring = Oscar.Nemo.ZZModRing(ZZ(49), false)
    matrixspace = matrix_space(ring, 3, 3)
    a = matrixspace()
    b = matrixspace()
    a_values = [1 1 0; 0 1 1; 1 0 1]
    b_values = [1 0 0; 0 2 0; 0 0 3]
    for j in 1:3, i in 1:3
        a[i, j] = ring(a_values[i, j])
        b[i, j] = ring(b_values[i, j])
    end
    g = ZZRingElem[ZZ(i + 1) for i in 1:3]
    start = 11
    stop = 138
    expected = copy(g)
    actual = copy(g)

    DeRham.finitediff_prodeval_linear!(
        a,
        b,
        start,
        stop,
        expected,
        matrixspace(),
        copy(g),
    )
    DeRham.sqrt_interval_product!(
        a,
        b,
        start,
        stop,
        actual,
        matrixspace(),
        copy(g),
    )
    @test actual == expected
end

@testset "fast interval product polynomial evaluation" begin
    for modulus in (101, ZZ(101)^25)
        modular_ring, = residue_ring(ZZ, modulus)
        poly_ring, x = polynomial_ring(modular_ring, "x"; cached = false)
        polynomial = poly_ring(3) * x^3 + poly_ring(5) * x + poly_ring(7)
        points = [-3, 0, 1, 104]

        if modulus <= ZZ(2)^64
            point_data = UInt[data(modular_ring(point)) for point in points]
            actual = DeRham._fast_nmod_evaluate(polynomial, point_data)
            expected = UInt[
                data(evaluate(polynomial, modular_ring(point))) for point in points
            ]
        else
            actual = DeRham._fast_fmpz_mod_evaluate(polynomial, points)
            expected = ZZRingElem[
                lift(evaluate(polynomial, modular_ring(point))) for point in points
            ]
        end
        @test actual == expected
    end
end

@testset "square-root interval product dispatch" begin
    ring, = residue_ring(ZZ, 101)
    matrixspace = matrix_space(ring, 2, 2)
    a = matrixspace()
    b = matrixspace()
    for j in 1:2, i in 1:2
        a[i, j] = ring(i + j)
        b[i, j] = ring(2i + j)
    end
    g = UInt[data(ring(i)) for i in 1:2]
    temp = matrixspace()
    g_temp = copy(g)
    sqrt_params = DeRham.ZetaFunctionParams(
        0,
        false,
        :pchunk,
        :invlex,
        false,
        true,
        false,
        false,
        false,
        :sqrt,
        4,
    )
    finite_difference_params = DeRham.ZetaFunctionParams(
        0,
        false,
        :pchunk,
        :invlex,
        false,
        true,
        false,
        false,
        false,
        :finitediff,
        1,
    )

    @test !DeRham._use_sqrt_interval_product(a, b, g, temp, g_temp, sqrt_params, 3)
    @test DeRham._use_sqrt_interval_product(a, b, g, temp, g_temp, sqrt_params, 4)
    @test !DeRham._use_sqrt_interval_product(
        a,
        b,
        g,
        temp,
        g_temp,
        finite_difference_params,
        4,
    )

    expected = copy(g)
    DeRham.finitediff_prodeval_linear!(a, b, 0, 2, expected, matrixspace(), copy(g))
    actual = copy(g)
    DeRham.prodeval_linear!(a, b, 0, 2, actual, matrixspace(), copy(g), sqrt_params)
    @test actual == expected

    unsupported_a = [1 2; 0 1]
    unsupported_b = [1 0; 3 2]
    unsupported_g = [1, 2]
    unsupported_temp = similar(unsupported_a)
    unsupported_g_temp = similar(unsupported_g)
    @test !DeRham._use_sqrt_interval_product(
        unsupported_a,
        unsupported_b,
        unsupported_g,
        unsupported_temp,
        unsupported_g_temp,
        sqrt_params,
        4,
    )

    expected_unsupported = copy(unsupported_g)
    DeRham.finitediff_prodeval_linear!(
        unsupported_a,
        unsupported_b,
        0,
        3,
        expected_unsupported,
        similar(unsupported_a),
        similar(unsupported_g),
    )
    actual_unsupported = copy(unsupported_g)
    DeRham.prodeval_linear!(
        unsupported_a,
        unsupported_b,
        0,
        3,
        actual_unsupported,
        similar(unsupported_a),
        similar(unsupported_g),
        sqrt_params,
    )
    @test actual_unsupported == expected_unsupported
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
    @test_throws ArgumentError zeta_coefficients(
        frontend;
        interval_product_min_length = big(typemax(Int)) + 1,
    )
    @test_throws ArgumentError DeRham.pregen_default_context(
        1,
        3,
        7,
        [-1];
        interval_product_algorithm = :invalid,
    )
    @test_throws ArgumentError DeRham.pregen_default_context(
        1,
        3,
        7,
        [-1];
        interval_product_min_length = 0,
    )
    @test_throws ArgumentError DeRham.pregen_default_context(
        1,
        3,
        7,
        [-1];
        interval_product_min_length = true,
    )

    context = DeRham.pregen_default_context(
        1,
        3,
        7,
        [-1];
        interval_product_algorithm = :sqrt,
        interval_product_min_length = 12,
        lazy = true,
    )
    @test context.interval_product_algorithm === :sqrt
    @test context.interval_product_min_length == 12
    @test DeRham._resolve_interval_product_options(context, nothing, nothing) ==
          (:sqrt, 12)
    @test DeRham._resolve_interval_product_options(context, :finitediff, 64) ==
          (:finitediff, 64)
    @test DeRham._resolve_interval_product_options(nothing, nothing, nothing) ==
          (:finitediff, 64)
end
