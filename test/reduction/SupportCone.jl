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
