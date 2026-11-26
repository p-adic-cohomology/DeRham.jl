@testset "matrix reduction engine" begin
    function test_engine(fastevaluation; kwargs...)
        coefficient_ring, = residue_ring(ZZ, 5)
        matspace = matrix_space(coefficient_ring, 1, 1)
        params = DeRham.ZetaFunctionParams(
            0,
            false,
            :depthfirst,
            :invlex,
            false,
            fastevaluation,
            false,
            false,
            false,
            :sqrt,
            1,
        )
        compute = _ -> begin
            constant = matspace()
            linear = matspace()
            constant[1, 1] = 1
            linear[1, 1] = 1
            [constant, linear]
        end
        pep = DeRham.LazyPEP{typeof(matspace())}(compute)
        context = DeRham.default_context(matspace, pep, params)
        directions = DeRham.ReductionDirectionSet(
            (u, k, v) -> length(u) == 1 && sum(v) == 1 && all(u .>= k .* v),
        )
        fill_vector! = (vector, key) -> fill!(vector, 1)
        preconfigured = DeRham.PreconfiguredReductionEngine(context)
        engine = DeRham.MatrixReductionEngine(
            DeRham.ReductionEngineGeometry(0, 1),
            directions,
            preconfigured,
            fill_vector!;
            fastevaluation,
            kwargs...,
        )
        return engine
    end

    @test_throws ArgumentError test_engine(true; linear_recurrence_algorithm = :invalid)
    @test_throws ArgumentError test_engine(true; sqrt_linear_recurrence_min_length = 0)

    for fastevaluation in (false, true)
        engine = test_engine(fastevaluation)
        @test engine.linear_recurrence_algorithm === :finitediff
        @test engine.sqrt_linear_recurrence_min_length === nothing
        first = DeRham.materialize(engine, DeRham.ReductionChainKey([3], 3))
        second = DeRham.materialize(engine, DeRham.ReductionChainKey([3], 3))
        state = @inferred DeRham._matrix_chain_state(engine, first)
        @test fieldtype(typeof(state), :vector) === typeof(engine.context.g)
        chains = [first, second]
        @test first !== second
        @test DeRham.reduce_batched!(engine, chains, [[1], [1]], 2) === chains
        survivors = DeRham.combine!(engine, 1)
        @test length(survivors) == 1
        u, vector = DeRham.finish!(engine, only(survivors))
        @test u == (1,)
        value = only(Array(vector))
        @test (value isa Integer ? value == 2 : Int(DeRham.lift(value)) == 2)
        @test isempty(engine.live)
    end
end
