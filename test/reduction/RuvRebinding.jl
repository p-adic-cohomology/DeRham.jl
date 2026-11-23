@testset "preconfigured R_uv rebinding" begin
    coefficient_ring, = residue_ring(ZZ, 5)
    matspace = matrix_space(coefficient_ring, 1, 1)
    first_matrix = matspace()
    second_matrix = matspace()
    first_matrix[1, 1] = 1
    second_matrix[1, 1] = 2

    pep = DeRham.PregenLazyPEP{typeof(first_matrix)}(_ -> [first_matrix])
    point = [1]
    original_matrix = only(pep[point])
    compute_calls = Ref(0)

    DeRham.bind_ruv_function!(pep, _ -> begin
        compute_calls[] += 1
        [second_matrix]
    end)

    @test compute_calls[] == 0
    @test only(pep[point]) === original_matrix
    @test original_matrix[1, 1] == 2
    @test compute_calls[] == 1

    backing_compute_calls = Ref(0)
    backing = DeRham.PregenLazyPEP{typeof(first_matrix)}(point -> begin
        matrix = matspace()
        matrix[1, 1] = point[1]
        [matrix]
    end)
    create = matrices -> begin
        matrix = matspace()
        copyto!(matrix, only(matrices))
        [matrix]
    end
    convert = (destination, source) -> copyto!(only(destination), only(source))
    cached_pep = DeRham.CachePEP{typeof(first_matrix),typeof(first_matrix)}(
        backing,
        create,
        convert,
        1,
    )
    first_point = [1]
    second_point = [2]
    cached_first = only(cached_pep[first_point])
    cached_second = only(cached_pep[second_point])
    @test cached_pep.tempV[] == first_point

    DeRham.bind_ruv_function!(cached_pep, point -> begin
        backing_compute_calls[] += 1
        matrix = matspace()
        matrix[1, 1] = point[1] + 2
        [matrix]
    end)

    @test backing_compute_calls[] == 0
    @test only(cached_pep[second_point]) === cached_second
    @test cached_second[1, 1] == 4
    @test only(cached_pep[first_point]) === cached_first
    @test cached_first[1, 1] == 3
    @test backing_compute_calls[] == 2
end

@testset "setup rebinds a reused context across polynomials" begin
    _, (x, y, z) = polynomial_ring(GF(5), ["x", "y", "z"])
    S = [0, 1, 2]
    first_polynomial = x^3 + y^3 + z^3
    second_polynomial = x^3 + y^3 + z^3 + 3 * x * y * z
    make_context =
        () -> DeRham.pregen_default_context(
            2,
            3,
            5,
            S;
            algorithm = :depthfirst,
            fastevaluation = false,
        )
    run =
        (f, context) -> DeRham.zeta_coefficients(
            f;
            S,
            algorithm = :depthfirst,
            context,
            changef = false,
            fastevaluation = false,
        )

    reused_context = make_context()
    reused_first = run(first_polynomial, reused_context)
    reused_second = run(second_polynomial, reused_context)
    fresh_first = run(first_polynomial, make_context())
    fresh_second = run(second_polynomial, make_context())

    @test reused_first !== false
    @test reused_second !== false
    @test reused_first == fresh_first
    @test reused_second == fresh_second
    @test !isempty(reused_context.Ruvs.Vs)
end
