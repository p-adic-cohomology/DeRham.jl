mutable struct MatrixReductionChainState{V}
    key::ReductionChainKey
    vector::V
end

"""
    MatrixReductionEngine(geometry, directions, context, fill_vector!; kwargs...)

An execution engine backed by the existing partially evaluated R_uv matrices
and reduction kernels. The engine owns a separate vector for every live chain;
its context supplies shared matrix and scratch buffers for sequential moves.
`fill_vector!(vector, key)` fills each new chain in the frontend's prescaled
normalization.

Keyword arguments:
`evaluation_kernel = eval_to_linear!` and `fastevaluation = false`.
"""
mutable struct MatrixReductionEngine{D,C,F,K,V} <: AbstractReductionEngine
    geometry::ReductionEngineGeometry
    directions::D
    context::C
    fill_vector!::F
    evaluation_kernel::K
    fastevaluation::Bool
    live::Dict{Int,MatrixReductionChainState{V}}
    next_identifier::Int
end

function MatrixReductionEngine(
    geometry::ReductionEngineGeometry,
    directions::ReductionDirectionSet,
    context::ControlledReductionContext,
    fill_vector!;
    evaluation_kernel = eval_to_linear!,
    fastevaluation = false,
)
    vector_type = typeof(_matrix_engine_vector(context.g))
    live = Dict{Int,MatrixReductionChainState{vector_type}}()
    return MatrixReductionEngine(
        geometry,
        directions,
        context,
        fill_vector!,
        evaluation_kernel,
        fastevaluation,
        live,
        1,
    )
end

function MatrixReductionEngine(
    geometry::ReductionEngineGeometry,
    directions::ReductionDirectionSet,
    preconfigured::PreconfiguredReductionEngine,
    fill_vector!;
    kwargs...,
)
    return MatrixReductionEngine(
        geometry,
        directions,
        preconfigured.context,
        fill_vector!;
        kwargs...,
    )
end

reduction_engine_geometry(engine::MatrixReductionEngine) = engine.geometry

_matrix_engine_vector(g::Vector{ZZRingElem}) = [ZZ(0) for _ in eachindex(g)]
_matrix_engine_vector(g) = similar(g)

function _matrix_chain_state(engine::MatrixReductionEngine, chain::ReductionChainHandle)
    state = get(engine.live, chain.identifier, nothing)
    state === nothing && throw(ArgumentError("the chain is not live in this engine"))
    state.key == chain.key || throw(ArgumentError("the chain handle is stale"))
    return state
end

function _materialize_matrix_chain!(engine, key, input)
    key.order >= engine.geometry.terminal_order ||
        throw(ArgumentError("cannot materialize a chain below the terminal order"))
    vector = _matrix_engine_vector(engine.context.g)
    if input !== nothing && applicable(engine.fill_vector!, vector, key, input)
        engine.fill_vector!(vector, key, input)
    elseif applicable(engine.fill_vector!, vector, key)
        engine.fill_vector!(vector, key)
    else
        throw(ArgumentError("the matrix engine requires an input record to fill a chain"))
    end
    identifier = engine.next_identifier
    engine.next_identifier += 1
    engine.live[identifier] = MatrixReductionChainState(key, vector)
    return ReductionChainHandle(key, identifier)
end

materialize(engine::MatrixReductionEngine, key::ReductionChainKey) =
    _materialize_matrix_chain!(engine, key, nothing)

materialize(engine::MatrixReductionEngine, key::ReductionChainKey, input) =
    _materialize_matrix_chain!(engine, key, input)

reduction_chain_key(engine::MatrixReductionEngine, chain::ReductionChainHandle) =
    _matrix_chain_state(engine, chain).key

function _combine_matrix_vectors!(destination::CuModArray, source)
    my_add!(destination, destination, source)
    return destination
end

function _combine_matrix_vectors!(destination::KaratsubaArray, source)
    my_add!(destination, destination, source)
    return destination
end

function _combine_matrix_vectors!(destination, source)
    destination .= destination .+ source
    return destination
end

function reduce!(engine::MatrixReductionEngine, chain::ReductionChainHandle, v, k::Integer)
    state = _matrix_chain_state(engine, chain)
    key = state.key
    direction = collect(v)
    k >= 1 || throw(ArgumentError("a reduction chunk must have positive length"))
    key.order - k >= engine.geometry.terminal_order ||
        throw(ArgumentError("a reduction chunk cannot pass the terminal order"))
    engine.directions(key.u, k, direction) ||
        throw(ArgumentError("the reduction direction is outside V_S(u, k)"))

    U = collect(key.u) .- k .* direction
    matrices = engine.context.Ruvs[direction]
    engine.evaluation_kernel(
        engine.context.B,
        engine.context.A,
        engine.context.temp,
        matrices,
        U,
        direction,
    )

    if engine.fastevaluation
        finitediff_prodeval_linear!(
            engine.context.B,
            engine.context.A,
            0,
            k - 1,
            state.vector,
            engine.context.temp,
            engine.context.g_temp,
        )
    else
        for i = 1:k
            my_mul!(engine.context.temp, engine.context.B, k - i)
            my_add!(engine.context.temp, engine.context.temp, engine.context.A)
            my_matvecmul!(engine.context.g_temp, engine.context.temp, state.vector)
            my_copy!(state.vector, engine.context.g_temp)
        end
    end

    next_key = ReductionChainKey(
        ntuple(i -> key.u[i] - k * direction[i], length(key.u)),
        key.order - k,
    )
    state.key = next_key
    chain.key = next_key
    return chain
end

function combine!(engine::MatrixReductionEngine, order::Integer)
    identifiers = sort!([
        identifier for (identifier, state) in engine.live if state.key.order == order
    ])
    surviving = Int[]
    seen = Dict{Tuple,Int}()
    for identifier in identifiers
        state = engine.live[identifier]
        if haskey(seen, state.key.u)
            survivor = engine.live[seen[state.key.u]]
            _combine_matrix_vectors!(survivor.vector, state.vector)
            delete!(engine.live, identifier)
        else
            seen[state.key.u] = identifier
            push!(surviving, identifier)
        end
    end
    return [
        ReductionChainHandle(engine.live[identifier].key, identifier) for
        identifier in surviving
    ]
end

function finish!(engine::MatrixReductionEngine, chain::ReductionChainHandle)
    state = _matrix_chain_state(engine, chain)
    state.key.order == engine.geometry.terminal_order ||
        throw(ArgumentError("only terminal chains can be finished"))
    delete!(engine.live, chain.identifier)
    return (state.key.u, state.vector)
end

public MatrixReductionEngine
