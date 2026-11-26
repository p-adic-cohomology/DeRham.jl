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

@testset "PChunk frontend at the CPU modulus boundary" begin
    curve = DeRham.fermat_hypersurface(3, 3, 2)
    basis = DeRham.get_basis_of_cohomology(curve)
    hodge, relative, series, _ = DeRham.calculate_precision_information(curve, basis)
    @test DeRham.zeta_coefficients(
        curve;
        algorithm = :pchunk,
        changef = false,
        precision_info = (hodge, relative, series, 64),
    ) == BigInt[2, 0, 1]
end

@testset "pchunk CPU representations" begin
    oversized_chunk = BigInt(typemax(Int)) + 1
    policy = DeRham.PChunkPolicy(3, [2], oversized_chunk)
    @test policy.maximum_chunk == typemax(Int)

    largest_allowed_modulus = BigInt(2)^64 - 1
    setup = _pchunk_execution_setup_at_modulus(largest_allowed_modulus)
    @test setup.backend == :cpu_uint
    @test setup.representation == :uint

    big_modulus_setup = _pchunk_execution_setup_at_modulus(BigInt(2)^64)
    @test big_modulus_setup.backend == :cpu_bigint
    @test big_modulus_setup.representation == :bigint
    big_modulus_context = big_modulus_setup.context_factory()
    @test big_modulus_context.A isa DeRham.ZZModMatrix
    @test big_modulus_context.g isa Vector{DeRham.ZZRingElem}

    forced_bigint_setup = _pchunk_execution_setup_at_modulus(5; always_use_bigints = true)
    @test forced_bigint_setup.backend == :cpu_bigint
    @test forced_bigint_setup.representation == :bigint
    forced_bigint_context = forced_bigint_setup.context_factory()
    @test forced_bigint_context.g isa Vector{DeRham.ZZRingElem}
    @test DeRham._matrix_engine_vector(forced_bigint_context.g) isa
          Vector{DeRham.ZZRingElem}
end
