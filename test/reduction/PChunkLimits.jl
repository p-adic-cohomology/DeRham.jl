@testset "pchunk representation limits" begin
    modulus = BigInt(2)^64
    ring, = residue_ring(ZZ, modulus)
    _, variables = polynomial_ring(ring, ["x0", "x1", "x2"])
    f = sum(variable^3 for variable in variables)
    fe = DeRham.ProjectiveHypersurface(f)
    P = DeRham.cohomology_quotient(fe)
    configuration = DeRham.reduction_setup(fe, [2], :pchunk)
    params = DeRham.ZetaFunctionParams(
        0,
        false,
        :pchunk,
        :invlex,
        false,
        true,
        false,
        false,
        false,
    )
    cache = DeRham.controlled_reduction_cache(2, 3, [2], params, P)
    g_length = binomial(DeRham._reduction_g_degree(P) + P.n, P.n)
    matrix_space = Oscar.matrix_space(ring, g_length, g_length)

    @test_throws ArgumentError DeRham.reduction_execution_setup(
        configuration,
        fe,
        f,
        P,
        matrix_space(),
        5,
        cache,
        params,
        nothing,
    )

    small_ring, = residue_ring(ZZ, 5)
    _, small_variables = polynomial_ring(small_ring, ["x0", "x1", "x2"])
    small_f = sum(variable^3 for variable in small_variables)
    small_fe = DeRham.ProjectiveHypersurface(small_f)
    small_P = DeRham.cohomology_quotient(small_fe)
    small_configuration = DeRham.reduction_setup(small_fe, [2], :pchunk)
    bigint_params = DeRham.ZetaFunctionParams(
        0,
        false,
        :pchunk,
        :invlex,
        false,
        true,
        true,
        false,
        false,
    )
    small_cache = DeRham.controlled_reduction_cache(2, 3, [2], bigint_params, small_P)
    small_g_length = binomial(DeRham._reduction_g_degree(small_P) + small_P.n, small_P.n)
    small_matrix_space = Oscar.matrix_space(small_ring, small_g_length, small_g_length)

    @test_throws ArgumentError DeRham.reduction_execution_setup(
        small_configuration,
        small_fe,
        small_f,
        small_P,
        small_matrix_space(),
        5,
        small_cache,
        bigint_params,
        nothing,
    )
end
