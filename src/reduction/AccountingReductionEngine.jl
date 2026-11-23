"""
    AccountingReductionResult

The symbolic cost and move counts collected by `AccountingReductionEngine`.
"""
struct AccountingReductionResult{C}
    cost::C
    live_chains::Int
    peak_live_chains::Int
    steps::Int
    chunks::Int
    duplicates_combined::Int
end

"""
    AccountingReductionEngine(geometry, directions)

An engine for the abstract reduction game. It tracks chain keys and operations
without storing or computing vectors. `directions` is the frontend's predicate
for membership in V_S(u, k).
"""
mutable struct AccountingReductionEngine{F,R,P} <: AbstractReductionEngine
    geometry::ReductionEngineGeometry
    directions::ReductionDirectionSet{F}
    cost_ring::R
    M::P
    T::P
    D::P
    cost::P
    live::Dict{Int,ReductionChainKey}
    events::Vector{NamedTuple}
    next_identifier::Int
    peak_live_chains::Int
    steps::Int
    chunks::Int
    duplicates_combined::Int
end

function AccountingReductionEngine(
    geometry::ReductionEngineGeometry,
    directions::ReductionDirectionSet{F},
) where {F}
    cost_ring, (M, T, D) = polynomial_ring(ZZ, ["M", "T", "D"])
    return AccountingReductionEngine{F,typeof(cost_ring),typeof(M)}(
        geometry,
        directions,
        cost_ring,
        M,
        T,
        D,
        cost_ring(0),
        Dict{Int,ReductionChainKey}(),
        NamedTuple[],
        1,
        0,
        0,
        0,
        0,
    )
end

AccountingReductionEngine(geometry::ReductionEngineGeometry, directions) =
    AccountingReductionEngine(geometry, ReductionDirectionSet(directions))

reduction_engine_geometry(engine::AccountingReductionEngine) = engine.geometry

function _accounting_key(engine::AccountingReductionEngine, chain::ReductionChainHandle)
    key = get(engine.live, chain.identifier, nothing)
    key === nothing && throw(ArgumentError("the chain is not live in this engine"))
    key == chain.key || throw(ArgumentError("the chain handle is stale"))
    return key
end

reduction_chain_key(engine::AccountingReductionEngine, chain::ReductionChainHandle) =
    _accounting_key(engine, chain)

function _accounting_event!(engine::AccountingReductionEngine, event::NamedTuple)
    push!(engine.events, event)
    return nothing
end

function materialize(engine::AccountingReductionEngine, key::ReductionChainKey)
    key.order >= engine.geometry.terminal_order ||
        throw(ArgumentError("cannot materialize a chain below the terminal order"))
    identifier = engine.next_identifier
    engine.next_identifier += 1
    engine.live[identifier] = key
    engine.peak_live_chains = max(engine.peak_live_chains, length(engine.live))
    _accounting_event!(
        engine,
        (; operation = :materialize, identifier, key, live_chains = length(engine.live)),
    )
    return ReductionChainHandle(key, identifier)
end

function reduce!(
    engine::AccountingReductionEngine,
    chain::ReductionChainHandle,
    v,
    k::Integer,
)
    key = _accounting_key(engine, chain)
    k >= 1 || throw(ArgumentError("a reduction chunk must have positive length"))
    key.order - k >= engine.geometry.terminal_order ||
        throw(ArgumentError("a reduction chunk cannot pass the terminal order"))
    engine.directions(key.u, k, v) ||
        throw(ArgumentError("the reduction direction is outside V_S(u, k)"))

    next_u = ntuple(i -> key.u[i] - k * v[i], length(key.u))
    next_key = ReductionChainKey(next_u, key.order - k)
    engine.live[chain.identifier] = next_key
    chain.key = next_key
    engine.cost += engine.M * k + engine.T
    engine.steps += k
    engine.chunks += 1
    _accounting_event!(
        engine,
        (
            operation = :reduce,
            identifier = chain.identifier,
            key_before = key,
            key_after = next_key,
            direction = Tuple(v),
            chunk_length = Int(k),
        ),
    )
    return chain
end

function reduce_batched!(engine::AccountingReductionEngine, chains, directions, k::Integer)
    k >= 1 || throw(ArgumentError("a reduction chunk must have positive length"))
    length(chains) == length(directions) ||
        throw(ArgumentError("each batched chain must have one direction"))
    for (chain, direction) in zip(chains, directions)
        reduce!(engine, chain, direction, k)
    end
    _accounting_event!(
        engine,
        (; operation = :reduce_batched, chunk_length = Int(k), count = length(chains)),
    )
    return chains
end

function combine!(engine::AccountingReductionEngine, order::Integer)
    identifiers =
        sort!([identifier for (identifier, key) in engine.live if key.order == order])
    surviving = Int[]
    seen = Dict{Tuple,Int}()
    duplicates = 0
    for identifier in identifiers
        key = engine.live[identifier]
        if haskey(seen, key.u)
            delete!(engine.live, identifier)
            duplicates += 1
        else
            seen[key.u] = identifier
            push!(surviving, identifier)
        end
    end
    engine.cost += engine.D * duplicates
    engine.duplicates_combined += duplicates
    _accounting_event!(
        engine,
        (
            operation = :combine,
            order = Int(order),
            duplicates = duplicates,
            live_chains = length(engine.live),
        ),
    )
    return [
        ReductionChainHandle(engine.live[identifier], identifier) for
        identifier in surviving
    ]
end

function finish!(engine::AccountingReductionEngine, chain::ReductionChainHandle)
    key = _accounting_key(engine, chain)
    key.order == engine.geometry.terminal_order ||
        throw(ArgumentError("only terminal chains can be finished"))
    delete!(engine.live, chain.identifier)
    _accounting_event!(
        engine,
        (;
            operation = :finish,
            identifier = chain.identifier,
            key,
            live_chains = length(engine.live),
        ),
    )
    return key
end

function accounting_result(engine::AccountingReductionEngine)
    return AccountingReductionResult(
        engine.cost,
        length(engine.live),
        engine.peak_live_chains,
        engine.steps,
        engine.chunks,
        engine.duplicates_combined,
    )
end

public AccountingReductionResult, AccountingReductionEngine, accounting_result
