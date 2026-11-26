function DeRham.my_mul!(output::Matrix{Int}, input::Matrix{Int}, scalar::Int)
    for j in axes(output, 2), i in axes(output, 1)
        output[i, j] = input[i, j] * scalar
    end
    return output
end

function DeRham.my_add!(output::Matrix{Int}, left::Matrix{Int}, right::Matrix{Int})
    for j in axes(output, 2), i in axes(output, 1)
        output[i, j] = left[i, j] + right[i, j]
    end
    return output
end

function DeRham.my_matvecmul!(output::Vector{Int}, matrix::Matrix{Int}, vector::Vector{Int})
    for i in axes(matrix, 1)
        output[i] = 0
        for j in axes(matrix, 2)
            output[i] += matrix[i, j] * vector[j]
        end
    end
    return output
end

@testset "square-root linear recurrence" begin
    for p in (101, 103)
        cases = (1, 2, 3, p - 1, p, p + 1, 3p + 2, 63, 64, 65)
        for modulus in (p, ZZ(p)^25)
            ring, = residue_ring(ZZ, modulus)
            matrixspace = matrix_space(ring, 3, 3)
            a = matrixspace()
            b = matrixspace()
            a_values = [1 1 0; 0 1 1; 1 0 1]
            b_values = [1 0 0; 0 2 0; 0 0 3]
            for j = 1:3, i = 1:3
                a[i, j] = ring(a_values[i, j])
                b[i, j] = ring(b_values[i, j])
            end

            big_modulus = modulus > ZZ(2)^64
            g =
                big_modulus ? ZZRingElem[ZZ(i + 1) for i = 1:3] :
                UInt[data(ring(i + 1)) for i = 1:3]
            for interval_length in cases
                start = 11
                stop = start + interval_length - 1
                expected = deepcopy(g)
                expected_temp = matrixspace()
                expected_g_temp = deepcopy(g)
                actual = deepcopy(g)
                actual_temp = matrixspace()
                actual_g_temp = deepcopy(g)

                DeRham.finitediff_prodeval_linear!(
                    a,
                    b,
                    start,
                    stop,
                    expected,
                    expected_temp,
                    expected_g_temp,
                )
                DeRham.sqrt_linear_recurrence!(
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

@testset "square-root linear recurrence buffer dimensions" begin
    ring, = residue_ring(ZZ, 101)
    matrixspace = matrix_space(ring, 2, 2)
    a = matrixspace()
    b = matrixspace()
    temp = matrixspace()
    g = UInt[1, 2]
    g_temp = copy(g)

    @test_throws DimensionMismatch DeRham.sqrt_linear_recurrence!(
        a,
        b,
        0,
        0,
        g,
        matrix_space(ring, 1, 1)(),
        g_temp,
    )
    @test_throws DimensionMismatch DeRham.sqrt_linear_recurrence!(
        a,
        b,
        0,
        0,
        g,
        temp,
        g_temp[1:1],
    )
    @test_throws DimensionMismatch DeRham.sqrt_linear_recurrence!(
        a,
        b,
        0,
        0,
        g[1:1],
        temp,
        g_temp,
    )
    other_ring, = residue_ring(ZZ, 103)
    @test_throws ArgumentError DeRham.sqrt_linear_recurrence!(
        a,
        b,
        0,
        0,
        g,
        matrix_space(other_ring, 2, 2)(),
        g_temp,
    )
end

@testset "square-root linear recurrence interval bounds" begin
    ring, = residue_ring(ZZ, 101)
    matrixspace = matrix_space(ring, 1, 1)
    a = matrixspace()
    b = matrixspace()
    b[1, 1] = ring(2)
    g = UInt[3]
    temp = matrixspace()
    g_temp = similar(g)

    DeRham.sqrt_linear_recurrence!(a, b, typemax(Int), typemax(Int), g, temp, g_temp)
    @test g == UInt[6]
    @test_throws OverflowError DeRham.sqrt_linear_recurrence!(
        a,
        b,
        typemin(Int),
        0,
        g,
        temp,
        g_temp,
    )
end

@testset "square-root linear recurrence with empty matrices" begin
    for (modulus, vector) in ((101, UInt[]), (big(101), ZZRingElem[]))
        ring, = residue_ring(ZZ, modulus)
        matrixspace = matrix_space(ring, 0, 0)
        @test DeRham.sqrt_linear_recurrence!(
            matrixspace(),
            matrixspace(),
            0,
            0,
            vector,
            matrixspace(),
            similar(vector),
        ) === vector
    end
end

@testset "square-root recurrence requires reduced UInt inputs" begin
    ring, = residue_ring(ZZ, 101)
    matrixspace = matrix_space(ring, 1, 1)
    input = UInt[101]
    @test_throws ArgumentError DeRham.sqrt_linear_recurrence!(
        matrixspace(),
        matrixspace(),
        0,
        0,
        input,
        matrixspace(),
        similar(input),
    )
end

@testset "square-root linear recurrence requires initialized BigInt entries" begin
    ring = Oscar.Nemo.ZZModRing(ZZ(49), false)
    matrixspace = matrix_space(ring, 1, 1)
    a = matrixspace()
    b = matrixspace()
    temp = matrixspace()
    g = ZZRingElem[ZZ(1)]
    g_temp = ZZRingElem[ZZ(0)]

    @test_throws ArgumentError DeRham.sqrt_linear_recurrence!(
        a,
        b,
        0,
        0,
        g,
        temp,
        similar(g),
    )
    @test_throws ArgumentError DeRham.sqrt_linear_recurrence!(
        a,
        b,
        0,
        0,
        similar(g),
        temp,
        g_temp,
    )
end

@testset "BGS reduces giant-step points modulo the working modulus" begin
    ring = Oscar.Nemo.ZZModRing(ZZ(49), false)
    matrixspace = matrix_space(ring, 3, 3)
    a = matrixspace()
    b = matrixspace()
    a_values = [1 1 0; 0 1 1; 1 0 1]
    b_values = [1 0 0; 0 2 0; 0 0 3]
    for j = 1:3, i = 1:3
        a[i, j] = ring(a_values[i, j])
        b[i, j] = ring(b_values[i, j])
    end
    g = ZZRingElem[ZZ(i + 1) for i = 1:3]
    start = 11
    stop = 138
    expected = deepcopy(g)
    actual = deepcopy(g)

    DeRham.finitediff_prodeval_linear!(
        a,
        b,
        start,
        stop,
        expected,
        matrixspace(),
        deepcopy(g),
    )
    DeRham.sqrt_linear_recurrence!(a, b, start, stop, actual, matrixspace(), deepcopy(g))
    @test actual == expected
end

@testset "BGS multipoint polynomial evaluation" begin
    points = [-3, 0, 1, 104]
    for modulus in (101, ZZ(101)^25)
        modular_ring, = residue_ring(ZZ, modulus)
        poly_ring, x = polynomial_ring(modular_ring, "x"; cached = false)
        polynomial = poly_ring(3) * x^3 + poly_ring(5) * x + poly_ring(7)
        second_polynomial = poly_ring(2) * x^2 + poly_ring(11) * x + poly_ring(13)
        word_modulus = modulus < ZZ(2)^64
        backend =
            word_modulus ? DeRham._NmodBGSEvaluationBackend() :
            DeRham._FmpzModBGSEvaluationBackend()
        coefficient = word_modulus ? data : lift
        buffer = DeRham._bgs_evaluation_buffer(backend, length(points))
        short_buffer = DeRham._bgs_evaluation_buffer(backend, length(points) - 1)
        output = word_modulus ? buffer : last(buffer)
        evaluation_points =
            DeRham._bgs_evaluation_points(backend, modular_ring, points, polynomial)
        try
            @test_throws DimensionMismatch DeRham._evaluate_bgs_polynomial!(
                short_buffer,
                polynomial,
                evaluation_points,
            )
            if !word_modulus
                @test_throws DimensionMismatch DeRham._evaluate_bgs_polynomial!(
                    (Int[], output),
                    polynomial,
                    evaluation_points,
                )
            end
            for poly in (polynomial, second_polynomial, zero(polynomial))
                actual = DeRham._evaluate_bgs_polynomial!(buffer, poly, evaluation_points)
                @test actual === output
                @test actual ==
                      [coefficient(evaluate(poly, modular_ring(point))) for point in points]
            end
        finally
            DeRham._release_bgs_evaluation_points!(backend, evaluation_points)
        end
    end
end

@testset "linear recurrence method selection" begin
    ring, = residue_ring(ZZ, 101)
    matrixspace = matrix_space(ring, 2, 2)
    a = matrixspace()
    b = matrixspace()
    for j = 1:2, i = 1:2
        a[i, j] = ring(i + j)
        b[i, j] = ring(2i + j)
    end
    g = UInt[data(ring(i)) for i = 1:2]
    temp = matrixspace()
    g_temp = copy(g)

    @test !DeRham._use_sqrt_linear_recurrence(a, b, g, temp, g_temp, :sqrt, 4, 3)
    @test DeRham._use_sqrt_linear_recurrence(a, b, g, temp, g_temp, :sqrt, 4, 4)
    @test !DeRham._use_sqrt_linear_recurrence(a, b, g, temp, g_temp, :finitediff, 1, 4)

    wide_matrixspace = matrix_space(ring, 28, 28)
    wide_a = wide_matrixspace()
    wide_b = wide_matrixspace()
    wide_g = UInt[data(ring(i)) for i = 1:28]
    wide_temp = wide_matrixspace()
    wide_g_temp = copy(wide_g)

    @test DeRham._use_sqrt_linear_recurrence(
        wide_a,
        wide_b,
        wide_g,
        wide_temp,
        wide_g_temp,
        :sqrt,
        64,
        64,
    )
    @test DeRham._use_sqrt_linear_recurrence(
        wide_a,
        wide_b,
        wide_g,
        wide_temp,
        wide_g_temp,
        :sqrt,
        64,
        28^2,
    )
    @test !DeRham._use_sqrt_linear_recurrence(a, b, g, temp, g_temp, :sqrt, nothing, 4096)
    @test DeRham._use_sqrt_linear_recurrence(
        wide_a,
        wide_b,
        wide_g,
        wide_temp,
        wide_g_temp,
        :sqrt,
        1,
        1,
    )
    @test !DeRham._use_sqrt_linear_recurrence(
        wide_a,
        wide_b,
        wide_g,
        wide_temp,
        wide_g_temp,
        :sqrt,
        nothing,
        28^2 - 1,
    )
    @test !DeRham._use_sqrt_linear_recurrence(
        wide_a,
        wide_b,
        wide_g,
        wide_temp,
        wide_g_temp,
        :sqrt,
        nothing,
        4096,
    )

    expected = copy(g)
    DeRham.finitediff_prodeval_linear!(a, b, 0, 2, expected, matrixspace(), copy(g))
    actual = copy(g)
    DeRham.evaluate_linear_recurrence!(a, b, 0, 2, actual, matrixspace(), copy(g), :sqrt, 4)
    @test actual == expected

    big_ring, = residue_ring(ZZ, ZZ(101)^25)
    big_matrixspace = matrix_space(big_ring, 2, 2)
    big_a = big_matrixspace()
    big_b = big_matrixspace()
    for j = 1:2, i = 1:2
        big_a[i, j] = big_ring(i + j)
        big_b[i, j] = big_ring(2i + j)
    end
    big_g = ZZRingElem[ZZ(i) for i = 1:2]
    big_expected = deepcopy(big_g)
    DeRham.finitediff_prodeval_linear!(
        big_a,
        big_b,
        0,
        3,
        big_expected,
        big_matrixspace(),
        deepcopy(big_g),
    )
    big_actual = deepcopy(big_g)
    big_sqrt_calls = Ref(0)
    DeRham._with_sqrt_linear_recurrence_test_hook(() -> (big_sqrt_calls[] += 1)) do
        DeRham.evaluate_linear_recurrence!(
            big_a,
            big_b,
            0,
            3,
            big_actual,
            big_matrixspace(),
            deepcopy(big_g),
            :sqrt,
            4,
        )
    end
    @test big_actual == big_expected
    @test big_sqrt_calls[] > 0

    unsupported_a = [1 2; 0 1]
    unsupported_b = [1 0; 3 2]
    unsupported_g = [1, 2]
    unsupported_temp = similar(unsupported_a)
    unsupported_g_temp = similar(unsupported_g)
    @test !DeRham._use_sqrt_linear_recurrence(
        unsupported_a,
        unsupported_b,
        unsupported_g,
        unsupported_temp,
        unsupported_g_temp,
        :sqrt,
        4,
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
    DeRham.evaluate_linear_recurrence!(
        unsupported_a,
        unsupported_b,
        0,
        3,
        actual_unsupported,
        similar(unsupported_a),
        similar(unsupported_g),
        :sqrt,
        4,
    )
    @test actual_unsupported == expected_unsupported
end

@testset "linear recurrence options" begin
    _, (x, y, z) = polynomial_ring(GF(7), ["x", "y", "z"])
    frontend = DeRham.ProjectiveHypersurface(x^3 + y^3 + z^3)
    @test_throws ArgumentError zeta_coefficients(
        frontend;
        linear_recurrence_algorithm = :invalid,
    )
    @test_throws ArgumentError zeta_coefficients(
        frontend;
        linear_recurrence_algorithm = missing,
    )
    @test_throws ArgumentError zeta_coefficients(
        frontend;
        sqrt_linear_recurrence_min_length = 0,
    )
    @test_throws ArgumentError zeta_coefficients(
        frontend;
        sqrt_linear_recurrence_min_length = true,
    )
    @test_throws ArgumentError zeta_coefficients(
        frontend;
        sqrt_linear_recurrence_min_length = big(typemax(Int)) + 1,
    )
    @test_throws ArgumentError DeRham.pregen_default_context(
        1,
        3,
        7,
        [-1];
        linear_recurrence_algorithm = :invalid,
    )
    @test_throws ArgumentError DeRham.pregen_default_context(
        1,
        3,
        7,
        [-1];
        linear_recurrence_algorithm = missing,
    )
    @test_throws ArgumentError DeRham.pregen_default_context(
        1,
        3,
        7,
        [-1];
        sqrt_linear_recurrence_min_length = 0,
    )
    @test_throws ArgumentError DeRham.pregen_default_context(
        1,
        3,
        7,
        [-1];
        sqrt_linear_recurrence_min_length = true,
    )

    context = DeRham.pregen_default_context(
        2,
        3,
        7,
        [-1];
        algorithm = :depthfirst,
        linear_recurrence_algorithm = :sqrt,
        sqrt_linear_recurrence_min_length = 1,
        lazy = true,
    )
    @test context.linear_recurrence_algorithm === :sqrt
    @test context.sqrt_linear_recurrence_min_length == 1
    @test DeRham._resolve_linear_recurrence_options(
        context,
        nothing,
        DeRham._UNSPECIFIED_SQRT_LINEAR_RECURRENCE_MIN_LENGTH,
    ) == (:sqrt, 1)
    @test DeRham._resolve_linear_recurrence_options(context, nothing, nothing) ==
          (:sqrt, nothing)
    @test DeRham._resolve_linear_recurrence_options(context, :finitediff, 64) ==
          (:finitediff, 64)
    @test DeRham._resolve_linear_recurrence_options(nothing, nothing, nothing) ==
          (:finitediff, nothing)

    automatic_context = DeRham.pregen_default_context(
        2,
        3,
        7,
        [-1];
        algorithm = :depthfirst,
        linear_recurrence_algorithm = :sqrt,
        lazy = true,
    )
    @test automatic_context.sqrt_linear_recurrence_min_length === nothing
    @test DeRham._resolve_linear_recurrence_options(
        automatic_context,
        nothing,
        DeRham._UNSPECIFIED_SQRT_LINEAR_RECURRENCE_MIN_LENGTH,
    ) == (:sqrt, nothing)

    sqrt_calls = Ref(0)
    DeRham._with_sqrt_linear_recurrence_test_hook(() -> (sqrt_calls[] += 1)) do
        inherited_coefficients = zeta_coefficients(
            frontend;
            algorithm = :depthfirst,
            changef = false,
            context = context,
        )
        @test inherited_coefficients == BigInt[7, 1, 1]
        @test sqrt_calls[] > 0

        sqrt_calls[] = 0
        automatic_cutoff_coefficients = zeta_coefficients(
            frontend;
            algorithm = :depthfirst,
            changef = false,
            context = context,
            sqrt_linear_recurrence_min_length = nothing,
        )
        @test automatic_cutoff_coefficients == inherited_coefficients
        @test sqrt_calls[] == 0

        overridden_coefficients = zeta_coefficients(
            frontend;
            algorithm = :depthfirst,
            changef = false,
            context = context,
            linear_recurrence_algorithm = :finitediff,
            sqrt_linear_recurrence_min_length = 1,
        )
        @test overridden_coefficients == inherited_coefficients
        @test sqrt_calls[] == 0
    end
    @test_throws ArgumentError zeta_coefficients(
        frontend;
        context = context,
        linear_recurrence_algorithm = :invalid,
    )
end

@testset "pre-generated CPU context residue representation" begin
    context = DeRham.pregen_default_context(
        2,
        3,
        7,
        [-1];
        algorithm = :depthfirst,
        always_use_bigints = true,
        lazy = true,
    )
    @test context.A isa DeRham.ZZModMatrix
    @test context.B isa DeRham.ZZModMatrix
    @test context.temp isa DeRham.ZZModMatrix
    @test context.g isa Vector{ZZRingElem}
    @test context.g_temp isa Vector{ZZRingElem}
end

@testset "square-root settings with fast evaluation disabled" begin
    _, (x, y, z) = polynomial_ring(GF(7), ["x", "y", "z"])
    frontend = DeRham.ProjectiveHypersurface(x^3 + y^3 + z^3)
    for algorithm in (:depthfirst, :pchunk)
        sqrt_calls = Ref(0)
        coefficients =
            DeRham._with_sqrt_linear_recurrence_test_hook(() -> (sqrt_calls[] += 1)) do
                zeta_coefficients(
                    frontend;
                    algorithm,
                    changef = false,
                    fastevaluation = false,
                    linear_recurrence_algorithm = :sqrt,
                    sqrt_linear_recurrence_min_length = 1,
                )
            end
        @test coefficients == BigInt[7, 1, 1]
        @test sqrt_calls[] == 0
    end
end

@testset "square-root linear recurrence in DepthFirst and PChunk" begin
    _, (x, y, z) = polynomial_ring(GF(7), ["x", "y", "z"])
    frontend = DeRham.ProjectiveHypersurface(x^3 + y^3 + z^3)
    default_coefficients =
        zeta_coefficients(frontend; algorithm = :depthfirst, changef = false)

    for algorithm in (:depthfirst, :pchunk), always_use_bigints in (false, true)
        sqrt_calls = Ref(0)
        sqrt_coefficients =
            DeRham._with_sqrt_linear_recurrence_test_hook(() -> (sqrt_calls[] += 1)) do
                zeta_coefficients(
                    frontend;
                    algorithm,
                    changef = false,
                    always_use_bigints,
                    linear_recurrence_algorithm = :sqrt,
                    sqrt_linear_recurrence_min_length = 1,
                )
            end

        @test sqrt_coefficients == default_coefficients
        @test sqrt_calls[] > 0
    end
end

@testset "square-root linear recurrence with a separate prime" begin
    _, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"])
    sqrt_calls = Ref(0)
    coefficients =
        DeRham._with_sqrt_linear_recurrence_test_hook(() -> (sqrt_calls[] += 1)) do
            zeta_coefficients(
                x^3 + y^3 + z^3,
                7;
                algorithm = :depthfirst,
                changef = false,
                fastevaluation = true,
                linear_recurrence_algorithm = :sqrt,
                sqrt_linear_recurrence_min_length = 1,
            )
        end
    @test coefficients == BigInt[7, 1, 1]
    @test sqrt_calls[] > 0
end
