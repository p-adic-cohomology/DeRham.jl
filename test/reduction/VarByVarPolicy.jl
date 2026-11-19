@testset "varbyvar reduction policy" begin
    geometry = DeRham.ReductionEngineGeometry(5, 2)
    directions = DeRham.ReductionDirectionSet(
        (u, k, v) ->
            sum(v) == 3 &&
            all(i -> v[i] >= 0 && u[i] >= k * v[i], eachindex(v)) &&
            (u[3] == 0 || v[3] > 0),
    )
    engine = DeRham.AccountingReductionEngine(geometry, directions)
    policy = DeRham.VarByVarPolicy(3)
    inputs = [
        (key = ([4, 3, 2], 5), terms = []),
        (key = ([4, 3, 2], 5), terms = []),
    ]
    terminal_chains = DeRham.ReductionChainKey[]
    on_terminal = chain -> begin
        @test count(event -> event.operation == :finish, engine.events) ==
              length(terminal_chains) + 1
        @test DeRham.accounting_result(engine).live_chains == 0
        push!(terminal_chains, chain)
    end

    @test isnothing(
        DeRham.run_reduction_policy!(policy, engine, [(5, inputs)], directions, on_terminal),
    )

    result = DeRham.accounting_result(engine)
    reductions = filter(event -> event.operation == :reduce, engine.events)
    combinations = filter(event -> event.operation == :combine, engine.events)
    first_combine = findfirst(event -> event.operation == :combine, engine.events)
    @test length(terminal_chains) == 1
    @test only(terminal_chains).u == (0, 0, 0)
    @test only(terminal_chains).order == geometry.terminal_order
    @test result.live_chains == 0
    @test result.peak_live_chains == 2
    @test result.steps == 4
    @test result.chunks == 4
    @test result.duplicates_combined == 1
    @test length(reductions) == 4
    @test all(event -> event.chunk_length == 1, reductions)
    @test reductions[1].key_before.order == reductions[2].key_before.order == 5
    @test reductions[3].key_before.order == 4
    @test reductions[4].key_before.order == 3
    @test [event.direction for event in reductions] ==
          [(0, 1, 2), (0, 1, 2), (1, 2, 0), (3, 0, 0)]
    @test [event.order for event in combinations] == [4, 3, 2]
    @test combinations[1].duplicates == 1
    @test count(
        event -> event.operation == :reduce,
        engine.events[1:first_combine],
    ) == 2
    @test all(
        event -> directions(event.key_before.u, event.chunk_length, event.direction),
        reductions,
    )
end
