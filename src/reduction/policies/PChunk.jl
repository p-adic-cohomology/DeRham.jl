function _pchunk_choose_direction(key::ReductionChainKey, degree, S)
    sum(key.u) >= degree ||
        throw(ArgumentError("pchunk has no direction for the chain key"))
    return rev_chooseV(collect(key.u), degree, S)
end

struct PChunkPolicy <: AbstractReductionPolicy
    degree::Int
    S::Vector{Int}
    maximum_chunk::Int
end

function PChunkPolicy(degree::Integer, S, maximum_chunk::Integer)
    maximum_chunk >= 1 || throw(ArgumentError("pchunk size must be positive"))
    return PChunkPolicy(Int(degree), Int.(collect(S)), Int(maximum_chunk))
end

function _pchunk_chunk_length(u, direction, maximum)
    chunk = maximum
    for i in eachindex(u)
        if direction[i] > 0
            chunk = min(chunk, fld(u[i], direction[i]))
        end
    end
    return chunk
end

function run_reduction_policy!(
    policy::PChunkPolicy,
    engine::AbstractReductionEngine,
    keys_by_order,
    directions::ReductionDirectionSet,
    on_terminal,
)
    terminal_order = reduction_engine_geometry(engine).terminal_order

    for (_, inputs) in keys_by_order
        for input in inputs
            key = ReductionChainKey(input.key[1], input.key[2])
            chain = materialize(engine, key, input)
            key = reduction_chain_key(engine, chain)
            key.order >= terminal_order ||
                throw(ArgumentError("pchunk chain is below the terminal order"))

            while key.order > terminal_order
                maximum_chunk = min(key.order - terminal_order, policy.maximum_chunk)
                direction = _pchunk_choose_direction(key, policy.degree, policy.S)
                chunk = _pchunk_chunk_length(key.u, direction, maximum_chunk)
                chunk >= 1 || throw(ArgumentError("pchunk selected an empty chunk"))
                directions(key.u, chunk, direction) ||
                    throw(ArgumentError("pchunk selected an invalid reduction direction"))
                reduce!(engine, chain, direction, chunk)
                key = reduction_chain_key(engine, chain)
            end

            key.order == terminal_order ||
                throw(ArgumentError("pchunk chain is below the terminal order"))
            on_terminal(finish!(engine, chain))
        end
    end

    return nothing
end

function reducetransform_pchunk(
    FT,
    N_m,
    S,
    f,
    P,
    pseudoInverseMat,
    p,
    cache,
    params,
    context,
    fe,
    configuration,
)
    result, _ = _reducetransform_matrix_policy(
        FT,
        S,
        f,
        P,
        pseudoInverseMat,
        p,
        cache,
        params,
        context,
        fe,
        configuration,
        PChunkPolicy(P.d, S, p),
    )
    return result
end

public PChunkPolicy
