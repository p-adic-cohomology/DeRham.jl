struct OpaqueDepthFirstHandle <: DeRham.AbstractReductionChainHandle
    identifier::Int
end

mutable struct OpaqueDepthFirstEngine <: DeRham.AbstractReductionEngine
    geometry::DeRham.ReductionEngineGeometry
    next_identifier::Int
    live::Dict{Int,Tuple{DeRham.ReductionChainKey,Int}}
end

function DeRham.reduction_engine_geometry(engine::OpaqueDepthFirstEngine)
    return engine.geometry
end

function DeRham.materialize(engine::OpaqueDepthFirstEngine, key::DeRham.ReductionChainKey)
    throw(ArgumentError("an input record is required"))
end

function DeRham.materialize(
    engine::OpaqueDepthFirstEngine,
    key::DeRham.ReductionChainKey,
    input,
)
    identifier = engine.next_identifier
    engine.next_identifier += 1
    engine.live[identifier] = (key, input.value)
    return OpaqueDepthFirstHandle(identifier)
end

function DeRham.reduction_chain_key(
    engine::OpaqueDepthFirstEngine,
    chain::OpaqueDepthFirstHandle,
)
    return engine.live[chain.identifier][1]
end

function DeRham.reduce!(
    engine::OpaqueDepthFirstEngine,
    chain::OpaqueDepthFirstHandle,
    direction,
    chunk_length::Integer,
)
    key, value = engine.live[chain.identifier]
    next_key = DeRham.ReductionChainKey(
        ntuple(i -> key.u[i] - chunk_length * direction[i], length(key.u)),
        key.order - chunk_length,
    )
    engine.live[chain.identifier] = (next_key, value)
    return chain
end

function DeRham.finish!(engine::OpaqueDepthFirstEngine, chain::OpaqueDepthFirstHandle)
    key, value = pop!(engine.live, chain.identifier)
    return (key.u, value)
end

@testset "depthfirst reduction policy" begin
    geometry = DeRham.ReductionEngineGeometry(5, 2)
    directions = DeRham.ReductionDirectionSet(
        (u, k, v) ->
            length(u) == length(v) && sum(v) == 3 && all(v .>= 0) && all(u .>= k .* v),
    )
    engine = DeRham.AccountingReductionEngine(geometry, directions)
    policy = DeRham.DepthFirstPolicy(3, Int[])
    inputs = [(key = ([4, 4, 4], 4), terms = []), (key = ([3, 4, 4], 4), terms = [])]

    terminal_chains = DeRham.ReductionChainKey[]
    on_terminal =
        chain -> begin
            finish_count = count(event -> event.operation == :finish, engine.events)
            @test finish_count == length(terminal_chains) + 1
            @test DeRham.accounting_result(engine).live_chains == 0
            push!(terminal_chains, chain)
        end
    @test isnothing(
        DeRham.run_reduction_policy!(
            policy,
            engine,
            [(4, inputs)],
            directions,
            on_terminal,
        ),
    )

    result = DeRham.accounting_result(engine)
    @test length(terminal_chains) == 2
    @test all(chain -> chain.order == geometry.terminal_order, terminal_chains)
    @test result.live_chains == 0
    @test result.peak_live_chains == 1
    @test result.steps == 4
    @test result.chunks == 4
    @test all(
        event -> directions(event.key_before.u, event.chunk_length, event.direction),
        filter(event -> event.operation == :reduce, engine.events),
    )

    opaque_engine = OpaqueDepthFirstEngine(
        DeRham.ReductionEngineGeometry(5, 4),
        1,
        Dict{Int,Tuple{DeRham.ReductionChainKey,Int}}(),
    )
    duplicate_key_inputs =
        [(key = ([4, 4, 4], 4), value = 2), (key = ([4, 4, 4], 4), value = 3)]
    opaque_chains = Int[]
    opaque_consumer = chain -> push!(opaque_chains, last(chain))
    @test isnothing(
        DeRham.run_reduction_policy!(
            policy,
            opaque_engine,
            [(4, duplicate_key_inputs)],
            directions,
            opaque_consumer,
        ),
    )
    @test sort(opaque_chains) == [2, 3]
end
