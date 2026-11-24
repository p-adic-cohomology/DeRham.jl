
"""
    _varbyvar_evaluation_kernel(params, d, n, modulus)

Select the existing evaluation implementation for varbyvar. The generic
evaluator supports the GPU array wrappers too; the named GPU kernels fuse the
coefficient evaluation and remain selected where they are the tuned path.
"""
function _varbyvar_evaluation_kernel(params, d, n, modulus)
    if params.use_gpu && !(ZZ(2)^25 < modulus < ZZ(2)^106)
        return eval_to_linear_gpu!
    elseif params.use_gpu && d == 3 && (n == 4 || n == 5)
        return eval_to_linear_gpu_karatsuba!
    end
    return eval_to_linear!
end

function reducetransform_varbyvar(
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
    result, execution_setup = _reducetransform_matrix_policy(
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
        VarByVarPolicy(P.d),
    )
    if params.verbose > 0
        println(
            "Created $(length(allpoints(execution_setup.Ruv))) of $(length(cache[P.d])) possible V",
        )
    end
    if params.verbose > 1
        println("V that were created: \n$(allpoints(execution_setup.Ruv))")
    end
    return result
end

function varbyvar_chooseV!(v, I, d)
    fill!(v, 0)
    remaining = d
    remaining <= 0 && return v

    for index in reverse(eachindex(I))
        if I[index] > 0
            amount = min(I[index], remaining)
            v[index] = amount
            remaining -= amount
            remaining == 0 && return v
        end
    end

    throw(ArgumentError("varbyvar could not choose a reduction direction"))
end

function varbyvar_chooseV(I, d)
    return varbyvar_chooseV!(zeros(Int, length(I)), I, d)
end

"""
    VarByVarPolicy(degree)

Experimental controlled-reduction policy that chooses reduction directions
one variable at a time. `degree` is the defining polynomial degree used to
select each direction.

Select this policy through `algorithm=:varbyvar` in `zeta_coefficients` or a
related coefficient-based entry point.
"""
struct VarByVarPolicy <: AbstractReductionPolicy
    degree::Int
end

function _varbyvar_chunk_length(u, direction, maximum)
    chunk = maximum
    for i in eachindex(u)
        if direction[i] > 0
            chunk = min(chunk, fld(u[i], direction[i]))
        end
    end
    return chunk
end

function run_reduction_policy!(
    policy::VarByVarPolicy,
    engine::AbstractReductionEngine,
    keys_by_order,
    directions::ReductionDirectionSet,
    on_terminal,
)
    terminal_order = reduction_engine_geometry(engine).terminal_order
    chains = nothing

    for (_, inputs) in keys_by_order
        for input in inputs
            key = ReductionChainKey(input.key[1], input.key[2])
            chain = materialize(engine, key, input)
            if chains === nothing
                chains = [chain]
            else
                push!(chains, chain)
            end
        end
    end

    chains === nothing && return nothing

    first_key = reduction_chain_key(engine, first(chains))
    seen_keys = Set{Tuple{typeof(first_key.u),Int}}()
    duplicate_orders = Set{Int}()
    direction = zeros(Int, length(first_key.u))
    chain_orders = Int[]
    sizehint!(chain_orders, length(chains))
    next_chains = similar(chains, 0)
    sizehint!(next_chains, length(chains))

    while !isempty(chains)
        for chain in chains
            key = reduction_chain_key(engine, chain)
            key.order == terminal_order && continue
            key.order > terminal_order ||
                throw(ArgumentError("varbyvar chain is below the terminal order"))

            varbyvar_chooseV!(direction, key.u, policy.degree)
            maximum_chunk = key.order - terminal_order
            chunk = _varbyvar_chunk_length(key.u, direction, maximum_chunk)
            chunk >= 1 || throw(ArgumentError("varbyvar selected an empty chunk"))
            directions(key.u, chunk, direction) ||
                throw(ArgumentError("varbyvar selected an invalid reduction direction"))
            reduce!(engine, chain, direction, chunk)
        end

        empty!(seen_keys)
        empty!(duplicate_orders)
        empty!(chain_orders)
        for chain in chains
            key = reduction_chain_key(engine, chain)
            push!(chain_orders, key.order)
            key_signature = (key.u, key.order)
            if key_signature in seen_keys
                push!(duplicate_orders, key.order)
            else
                push!(seen_keys, key_signature)
            end
        end

        empty!(next_chains)
        for (chain, order) in zip(chains, chain_orders)
            order in duplicate_orders && continue
            if order == terminal_order
                on_terminal(finish!(engine, chain))
            else
                push!(next_chains, chain)
            end
        end

        for order in duplicate_orders
            survivors = combine!(engine, order)
            if order == terminal_order
                for chain in survivors
                    on_terminal(finish!(engine, chain))
                end
            else
                append!(next_chains, survivors)
            end
        end
        chains, next_chains = next_chains, chains
    end

    return nothing
end

public VarByVarPolicy
