struct DepthFirstPolicy <: AbstractReductionPolicy
    degree::Int
    S::Vector{Int}
end

function _depthfirst_chunk_length(u, direction, maximum)
    chunk = maximum
    for i = 1:length(u)
        if direction[i] > 0
            chunk = min(chunk, fld(u[i], direction[i]))
        end
    end
    return chunk
end

function run_reduction_policy!(
    policy::DepthFirstPolicy,
    engine::AbstractReductionEngine,
    keys_by_order,
    directions::ReductionDirectionSet,
    on_terminal,
)
    terminal_order = reduction_engine_geometry(engine).terminal_order

    for (_, inputs) in keys_by_order
        for input in inputs
            chain =
                materialize(engine, ReductionChainKey(input.key[1], input.key[2]), input)
            key = reduction_chain_key(engine, chain)
            while key.order > terminal_order
                direction = rev_chooseV(collect(key.u), policy.degree, policy.S)
                maximum_chunk = key.order - terminal_order
                chunk = _depthfirst_chunk_length(key.u, direction, maximum_chunk)
                chunk >= 1 || throw(ArgumentError("depthfirst selected an empty chunk"))
                directions(key.u, chunk, direction) || throw(
                    ArgumentError("depthfirst selected an invalid reduction direction"),
                )
                reduce!(engine, chain, direction, chunk)
                key = reduction_chain_key(engine, chain)
            end
            on_terminal(finish!(engine, chain))
        end
    end

    return nothing
end

public DepthFirstPolicy
