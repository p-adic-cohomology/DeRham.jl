@testset "pchunk direction selection" begin
    R, variables = polynomial_ring(GF(5), ["x$i" for i = 0:4])
    fe = DeRham.ProjectiveHypersurface(sum(variable^3 for variable in variables))
    S = [0, 1, 2]
    directions = DeRham.ReductionDirectionSet(DeRham.reduction_directions(fe, S))
    key = DeRham.ReductionChainKey([11, 8, 13, 9, 6], 7)

    direction = DeRham._pchunk_choose_direction(key, 3, S)
    @test direction == [1, 1, 1, 0, 0]
    chunk = DeRham._depthfirst_chunk_length(collect(key.u), direction, 3)
    @test chunk == 3
    @test directions(key.u, chunk, direction)
    @test any(!iszero(rem(exponent, 2)) for exponent in key.u)

    invalid_key = DeRham.ReductionChainKey([0, 0, 0, 0, 0], 7)
    @test_throws ArgumentError DeRham._pchunk_choose_direction(invalid_key, 3, S)
end

@testset "pchunk reduction policy" begin
    geometry = DeRham.ReductionEngineGeometry(3, 1)
    S = [0, 1, 2]
    directions = DeRham.ReductionDirectionSet(
        (u, k, v) ->
            length(u) == length(v) &&
            sum(v) == 3 &&
            all(v .>= 0) &&
            all(u .>= k .* v) &&
            all(i -> u[i+1] == 0 || v[i+1] != 0, S),
    )
    engine = DeRham.AccountingReductionEngine(geometry, directions)
    policy = DeRham.PChunkPolicy(3, S, 3)
    inputs = [(key = ([11, 8, 13, 9, 6], 7), terms = [])]
    terminal_chains = DeRham.ReductionChainKey[]
    on_terminal = chain -> push!(terminal_chains, chain)

    @test isnothing(
        DeRham.run_reduction_policy!(
            policy,
            engine,
            [(7, inputs)],
            directions,
            on_terminal,
        ),
    )

    result = DeRham.accounting_result(engine)
    reductions = filter(event -> event.operation == :reduce, engine.events)
    @test length(terminal_chains) == 1
    @test only(terminal_chains) == DeRham.ReductionChainKey([5, 2, 7, 9, 6], 1)
    @test result.live_chains == 0
    @test result.steps == 6
    @test result.chunks == 2
    @test [event.chunk_length for event in reductions] == [3, 3]
    @test all(event -> event.direction == (1, 1, 1, 0, 0), reductions)
    @test all(
        event -> directions(event.key_before.u, event.chunk_length, event.direction),
        reductions,
    )
end
