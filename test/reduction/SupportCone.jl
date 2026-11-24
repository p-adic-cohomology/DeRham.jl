function supportcone_test_directions(degree, S)
    return DeRham.ReductionDirectionSet(
        (u, k, v) ->
            sum(v) == degree &&
            all(index -> u[index] >= k * v[index], eachindex(u)) &&
            all(index -> u[index+1] == 0 || v[index+1] != 0, S),
    )
end

@testset "SupportCone planner" begin
    fermat_support = [(4, 0, 0), (0, 4, 0), (0, 0, 4)]
    fermat = DeRham.SupportConePlanner(fermat_support)
    directions = supportcone_test_directions(4, Int[])
    key = DeRham.ReductionChainKey((8, 8, 8), 6)
    plan = DeRham._supportcone_plan(fermat, key, 1, directions)

    @test fermat.support == [(0, 0, 4), (0, 4, 0), (4, 0, 0)]
    @test length(fermat.subcones) == 1
    @test plan !== nothing
    @test plan.subcone == (1, 2, 3)
    @test [(chunk.direction, chunk.steps) for chunk in plan.chunks] ==
          [((0, 0, 4), 1), ((0, 4, 0), 2), ((4, 0, 0), 2)]
    @test sum(chunk.steps for chunk in plan.chunks) == key.order - 1
    repeated_plan = DeRham._supportcone_plan(
        fermat,
        DeRham.ReductionChainKey(collect(key.u), key.order),
        1,
        directions,
    )
    @test (
        plan.subcone,
        [(chunk.direction, chunk.steps) for chunk in plan.chunks],
    ) == (
        repeated_plan.subcone,
        [(chunk.direction, chunk.steps) for chunk in repeated_plan.chunks],
    )

    current_u = collect(key.u)
    for chunk in plan.chunks
        @test directions(Tuple(current_u), chunk.steps, chunk.direction)
        current_u .-= chunk.steps .* chunk.direction
    end
    @test current_u == [0, 0, 4]

    overlapping_support = [(4, 0, 0), (0, 4, 0), (0, 0, 4), (2, 2, 0)]
    overlapping = DeRham.SupportConePlanner(overlapping_support)
    overlapping_key = DeRham.ReductionChainKey((8, 8, 4), 5)
    overlapping_plan =
        DeRham._supportcone_plan(overlapping, overlapping_key, 1, directions)

    @test length(overlapping.subcones) == 3
    @test DeRham._supportcone_floors(overlapping, overlapping.subcones[1], overlapping_key.u) ==
          BigInt[1, 0, 4]
    @test DeRham._supportcone_floors(overlapping, overlapping.subcones[2], overlapping_key.u) ==
          BigInt[1, 2, 2]
    @test overlapping_plan.subcone == (1, 2, 3)
    @test [(chunk.direction, chunk.steps) for chunk in overlapping_plan.chunks] ==
          [((0, 0, 4), 1), ((2, 2, 0), 3)]
    @test sum(chunk.steps for chunk in overlapping_plan.chunks) ==
          overlapping_key.order - 1
    overlapping_u = collect(overlapping_key.u)
    for chunk in overlapping_plan.chunks
        @test directions(Tuple(overlapping_u), chunk.steps, chunk.direction)
        overlapping_u .-= chunk.steps .* chunk.direction
    end
    @test overlapping_u == [2, 2, 0]

    incomplete = DeRham.SupportConePlanner([(4, 0, 0), (0, 4, 0)])
    @test isempty(incomplete.subcones)
    @test DeRham._supportcone_plan(incomplete, overlapping_key, 1, directions) === nothing

    restricted_directions = supportcone_test_directions(4, [0])
    restricted_key = DeRham.ReductionChainKey((4, 4, 4), 4)
    @test DeRham._supportcone_plan(fermat, restricted_key, 1, restricted_directions) === nothing
end

@testset "SupportCone policy" begin
    directions = supportcone_test_directions(2, Int[])
    geometry = DeRham.ReductionEngineGeometry(2, 1)
    planner = DeRham.SupportConePlanner([(2, 0), (0, 2)])
    policy = DeRham.SupportConePolicy(planner, 2, Int[])
    engine = DeRham.AccountingReductionEngine(geometry, directions)
    inputs = [
        (4, [(key = ([4, 4], 4), terms = [])]),
        (3, [(key = ([4, 2], 3), terms = [])]),
    ]
    terminal_keys = DeRham.ReductionChainKey[]
    on_terminal = key -> push!(terminal_keys, key)

    @test isnothing(
        DeRham.run_reduction_policy!(policy, engine, inputs, directions, on_terminal),
    )
    result = DeRham.accounting_result(engine)
    @test terminal_keys == [DeRham.ReductionChainKey((2, 0), 1)]
    @test result.live_chains == 0
    @test result.peak_live_chains == 2
    @test result.steps == 3
    @test result.chunks == 3
    @test result.duplicates_combined == 1
    @test result.fallback_keys == 0
    @test result.fallback_steps == 0

    events = engine.events
    duplicate_index = findfirst(
        event -> event.operation == :combine && event.duplicates == 1,
        events,
    )
    @test events[duplicate_index - 1].operation == :reduce
    @test events[duplicate_index - 1].key_after.order == 3
    @test events[duplicate_index + 1].operation == :reduce
    @test events[duplicate_index + 1].key_before.order == 3
    @test all(
        event -> directions(event.key_before.u, event.chunk_length, event.direction),
        filter(event -> event.operation == :reduce, events),
    )
    @test all(
        event -> event.operation != :finish ||
                 event.key.order == geometry.terminal_order,
        events,
    )

    fallback_engine = DeRham.AccountingReductionEngine(geometry, directions)
    fallback_policy = DeRham.SupportConePolicy(
        DeRham.SupportConePlanner([(2, 0)]),
        2,
        Int[],
    )
    fallback_inputs = [(4, [(key = ([4, 4], 4), terms = [])])]
    fallback_terminals = DeRham.ReductionChainKey[]
    @test isnothing(
        DeRham.run_reduction_policy!(
            fallback_policy,
            fallback_engine,
            fallback_inputs,
            directions,
            key -> push!(fallback_terminals, key),
        ),
    )
    fallback_result = DeRham.accounting_result(fallback_engine)
    @test fallback_terminals == [DeRham.ReductionChainKey((0, 2), 1)]
    @test fallback_result.live_chains == 0
    @test fallback_result.steps == 3
    @test fallback_result.fallback_keys ==
          count(event -> event.operation == :fallback, fallback_engine.events)
    @test fallback_result.fallback_keys == 2
    @test fallback_result.fallback_steps == 3
end
