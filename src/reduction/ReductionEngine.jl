"""
    AbstractReductionEngine

An execution engine for controlled reduction. An engine owns materialized
chains, their vector representation and placement, and the operations that
move and combine them. The engine receives only the reduction degree and
terminal order from the cohomology quotient.
"""
abstract type AbstractReductionEngine end

"""
    AbstractReductionPolicy

A strategy that chooses which chain to reduce, its direction and chunk length,
and when chains at one order are combined. A policy operates through the
`AbstractReductionEngine` interface and does not inspect vector data or
representation.
"""
abstract type AbstractReductionPolicy end

"""
    ReductionEngineGeometry(reduction_degree, terminal_order)

The quotient data an engine needs: the degree of its chain vectors and the
order at which reduction stops.
"""
struct ReductionEngineGeometry
    reduction_degree::Int
    terminal_order::Int
end

ReductionEngineGeometry(P::CohomologyQuotient) =
    ReductionEngineGeometry(_reduction_g_degree(P), terminal_pole_order(P))

"""
    ReductionChainKey(u, order)

An immutable identity for one reduction chain. `u` is the exponent vector after
the frontend's shift and degree adjustment; `order` identifies the chain's
current reduction order.
"""
struct ReductionChainKey{N}
    u::NTuple{N,Int}
    order::Int
end

function ReductionChainKey(u::AbstractVector{<:Integer}, order::Integer)
    return ReductionChainKey(Tuple(Int.(u)), Int(order))
end

"""
    AbstractReductionChainHandle

An opaque reference to a live chain owned by a reduction engine.
"""
abstract type AbstractReductionChainHandle end

"""
    ReductionChainHandle(key, identifier)

A generic chain handle containing its current key and an engine-owned
identifier. `reduce!` updates the key as the chain moves. Engines may use their
own subtypes of `AbstractReductionChainHandle` when a different handle
representation is useful.
"""
mutable struct ReductionChainHandle{K<:ReductionChainKey,I} <: AbstractReductionChainHandle
    key::K
    identifier::I
end

"""
    ReductionDirectionSet(contains)

A predicate for membership in the frontend's direction collection V_S(u, k).
Call it as `directions(u, k, v)` to test a candidate direction.
"""
struct ReductionDirectionSet{F}
    contains::F
end

(directions::ReductionDirectionSet)(u, k, v) = directions.contains(u, k, v)

"""
    reduction_engine_geometry(engine)

Return the `ReductionEngineGeometry` consumed by the engine. It contains only
the reduction degree and terminal order, not the cohomology quotient.
"""
function reduction_engine_geometry end

"""
    reduction_chain_key(engine, chain)

Return the current key of a live chain without inspecting its engine-owned
handle representation.
"""
function reduction_chain_key end

"""
    materialize(engine, key[, input])

Allocate a chain for `key`, fill its vector with the frontend input in the
prescaled normalization, and return an opaque live-chain handle. When supplied,
`input` is an opaque frontend record associated with the chain. The engine
chooses the vector representation and placement.
"""
function materialize end

materialize(engine::AbstractReductionEngine, key::ReductionChainKey, input) =
    materialize(engine, key)

"""
    reduce!(engine, chain, v, k)

Apply a reduction chunk of length `k` in direction `v` to `chain`, preserving
the prescaled normalization used by the frontend, and update the handle's key.
"""
function reduce! end

"""
    reduce_batched!(engine, chains, directions, k)

Reduce each chain in `chains` by a chunk of common length `k`, using the
corresponding direction in `directions`. Version one may implement this as a
loop over `reduce!`.
"""
function reduce_batched! end

function reduce_batched!(engine::AbstractReductionEngine, chains, directions, k::Integer)
    k >= 1 || throw(ArgumentError("a reduction chunk must have positive length"))
    length(chains) == length(directions) ||
        throw(ArgumentError("each batched chain must have one direction"))
    for (chain, direction) in zip(chains, directions)
        reduce!(engine, chain, direction, k)
    end
    return chains
end

"""
    combine!(engine, order)

For the given order, combine live chains with equal exponent vectors. Vector
engines add their chain vectors; accounting engines record the duplicate
removals. Return the surviving handles at that order. The engine owns the
resulting live-chain bookkeeping.
"""
function combine! end

"""
    finish!(engine, chain)

Read back a terminal chain as its exponent vector and reduced vector, then free
the live chain from the engine. Return the engine's terminal-chain result; a
matrix engine returns its exponent vector and reduced vector, while an
accounting engine returns only the key. The frontend converts matrix results to
a polynomial.
"""
function finish! end

"""
    run_reduction_policy!(policy, engine, keys_by_order, directions, on_terminal)

Run `policy` against `engine` using frontend input records grouped by order and
the frontend's `V_S` direction predicate. The policy passes each input record
to the engine when materializing its chain, then chooses directions and chunk
lengths through the engine interface. It does not inspect input records,
matrices, vectors, representations, or device placement. Each finished chain
is passed to `on_terminal` immediately after the engine finishes it, and
implementations do not retain terminal chains.
"""
function run_reduction_policy! end

public AbstractReductionEngine,
    AbstractReductionPolicy,
    ReductionEngineGeometry,
    ReductionChainKey,
    AbstractReductionChainHandle,
    ReductionChainHandle,
    ReductionDirectionSet,
    reduction_engine_geometry,
    reduction_chain_key,
    materialize,
    reduce!,
    reduce_batched!,
    combine!,
    finish!,
    run_reduction_policy!
