@testset "accounting reduction engine" begin
    geometry = DeRham.ReductionEngineGeometry(5, 2)
    directions = DeRham.ReductionDirectionSet(
        (u, k, v) ->
            length(u) == length(v) && sum(v) == 3 && all(v .>= 0) && all(u .>= k .* v),
    )
    engine = DeRham.AccountingReductionEngine(geometry, directions)
    key = DeRham.ReductionChainKey([5, 4, 3], 5)
    chain = DeRham.materialize(engine, key)

    @test_throws ArgumentError DeRham.materialize(
        engine,
        DeRham.ReductionChainKey([1, 1, 1], 1),
    )
    @test_throws ArgumentError DeRham.finish!(engine, chain)
    @test_throws ArgumentError DeRham.reduce!(engine, chain, [1, 1, 0], 1)
    @test_throws ArgumentError DeRham.reduce!(engine, chain, [1, 1, 1], 4)
    @test DeRham.reduce!(engine, chain, [1, 1, 1], 1) === chain
    @test chain.key == DeRham.ReductionChainKey([4, 3, 2], 4)
    @test DeRham.reduce!(engine, chain, [1, 1, 1], 2) === chain
    @test chain.key == DeRham.ReductionChainKey([2, 1, 0], 2)

    @test DeRham.finish!(engine, chain) == chain.key

    result = DeRham.accounting_result(engine)
    @test result.cost == 3 * engine.M + 2 * engine.T
    @test result.live_chains == 0
    @test result.peak_live_chains == 1
    @test result.steps == 3
    @test result.chunks == 2
    @test result.duplicates_combined == 0
    @test [event.operation for event in engine.events] == [:materialize, :reduce, :reduce, :finish]

    duplicates = DeRham.AccountingReductionEngine(geometry, directions)
    first = DeRham.materialize(duplicates, DeRham.ReductionChainKey([1, 2, 3], 2))
    second = DeRham.materialize(duplicates, DeRham.ReductionChainKey([1, 2, 3], 2))
    unique = DeRham.materialize(duplicates, DeRham.ReductionChainKey([3, 2, 1], 2))
    survivors = DeRham.combine!(duplicates, 2)
    @test length(survivors) == 2
    @test all(handle.identifier != second.identifier for handle in survivors)
    @test any(handle.identifier == first.identifier for handle in survivors)
    @test any(handle.identifier == unique.identifier for handle in survivors)
    duplicate_result = DeRham.accounting_result(duplicates)
    @test duplicate_result.cost == duplicates.D
    @test duplicate_result.live_chains == 2
    @test duplicate_result.peak_live_chains == 3
    @test duplicate_result.duplicates_combined == 1

    batch = DeRham.AccountingReductionEngine(geometry, directions)
    batch_chains = [
        DeRham.materialize(batch, DeRham.ReductionChainKey([4, 4, 4], 4)),
        DeRham.materialize(batch, DeRham.ReductionChainKey([3, 4, 4], 4)),
    ]
    @test DeRham.reduce_batched!(batch, batch_chains, [[1, 1, 1], [1, 1, 1]], 1) ===
          batch_chains
    @test DeRham.accounting_result(batch).cost == 2 * batch.M + 2 * batch.T
    @test DeRham.accounting_result(batch).steps == 2
    @test DeRham.accounting_result(batch).chunks == 2
    @test last(batch.events).operation == :reduce_batched
end
