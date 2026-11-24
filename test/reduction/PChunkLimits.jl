function _pchunk_execution_setup_at_modulus(modulus; always_use_bigints = false)
    ring, = residue_ring(ZZ, BigInt(modulus))
    _, variables = polynomial_ring(ring, ["x0", "x1", "x2"])
    f = sum(variable^3 for variable in variables)
    fe = DeRham.ProjectiveHypersurface(f)
    P = DeRham.cohomology_quotient(fe)
    configuration =
        DeRham.reduction_setup(fe, [2], :pchunk; always_use_bigints = always_use_bigints)
    params = DeRham.ZetaFunctionParams(
        0,
        false,
        :pchunk,
        :invlex,
        false,
        true,
        always_use_bigints,
        false,
        false,
    )
    cache = DeRham.controlled_reduction_cache(2, 3, [2], params, P)
    g_length = binomial(DeRham._reduction_g_degree(P) + P.n, P.n)
    matrix_space = Oscar.matrix_space(ring, g_length, g_length)
    return DeRham.reduction_execution_setup(
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
end

@testset "pchunk representation limits" begin
    oversized_chunk = BigInt(typemax(Int)) + 1
    policy = DeRham.PChunkPolicy(3, [2], oversized_chunk)
    @test policy.maximum_chunk == typemax(Int)

    largest_allowed_modulus = BigInt(2)^64 - 1
    setup = _pchunk_execution_setup_at_modulus(largest_allowed_modulus)
    @test setup.backend == :cpu_uint

    rejected_setup = try
        _pchunk_execution_setup_at_modulus(BigInt(2)^64)
        nothing
    catch err
        err
    end
    @test rejected_setup isa ArgumentError
    if rejected_setup isa ArgumentError
        @test occursin(
            "pchunk currently supports only small-modulus CPU execution",
            sprint(showerror, rejected_setup),
        )
    end

    @test_throws ArgumentError _pchunk_execution_setup_at_modulus(
        5;
        always_use_bigints = true,
    )
end
