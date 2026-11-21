
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

struct VarByVarEngineState{E}
    engine::E
end

function _varbyvar_engine_factory(execution_setup, fe, P, S, cache, params, directions)
    return () -> begin
        fill_vector! =
            (vector, key, input) ->
                reduction_input_vector!(vector, fe, input, S, cache, params)
        preconfigured = if execution_setup.preconfigured_engine === nothing
            PreconfiguredReductionEngine(execution_setup.context_factory())
        else
            execution_setup.preconfigured_engine
        end
        engine = MatrixReductionEngine(
            ReductionEngineGeometry(P),
            directions,
            preconfigured,
            fill_vector!;
            evaluation_kernel = execution_setup.kernel,
            fastevaluation = params.fastevaluation,
        )
        return VarByVarEngineState(engine)
    end
end

function _reducepoly_varbyvar!(
    runtime::VarByVarEngineState,
    pol,
    S,
    P,
    fe,
    PR,
    cache,
    params,
    policy,
    directions,
)
    keys_by_order = reduction_input_keys(fe, pol, S)
    terminal_polynomial = Ref(PR())
    on_terminal =
        chain -> begin
            terminal_polynomial[] = reduction_add_terminal_chain(
                fe,
                terminal_polynomial[],
                chain,
                PR,
                params,
            )
        end
    run_reduction_policy!(policy, runtime.engine, keys_by_order, directions, on_terminal)
    return _finish_reduction_terminal_polynomial(fe, terminal_polynomial[], PR, S)
end

function _reducetransform_varbyvar_engine(
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
)
    execution_setup = reduction_execution_setup(
        configuration,
        fe,
        f,
        P,
        pseudoInverseMat,
        p,
        cache,
        params,
        context,
    )
    directions = ReductionDirectionSet(reduction_directions(fe, S))
    policy = VarByVarPolicy(P.d)
    engine_factory =
        _varbyvar_engine_factory(execution_setup, fe, P, S, cache, params, directions)
    result = similar(FT)

    if execution_setup.threads
        runtime_tlv = OhMyThreads.TaskLocalValue{VarByVarEngineState}(engine_factory)
        Threads.@threads for i in eachindex(FT)
            local runtime = runtime_tlv[]
            if 0 < params.verbose
                println("Reducing vector $i in thread $(Threads.threadid())")
                @time result[i] = _reducepoly_varbyvar!(
                    runtime,
                    FT[i],
                    S,
                    P,
                    fe,
                    parent(f),
                    cache,
                    params,
                    policy,
                    directions,
                )
            else
                result[i] = _reducepoly_varbyvar!(
                    runtime,
                    FT[i],
                    S,
                    P,
                    fe,
                    parent(f),
                    cache,
                    params,
                    policy,
                    directions,
                )
            end
        end
    else
        runtime = engine_factory()
        for i in eachindex(FT)
            if 0 < params.verbose
                println("Reducing vector $i")
                @time result[i] = _reducepoly_varbyvar!(
                    runtime,
                    FT[i],
                    S,
                    P,
                    fe,
                    parent(f),
                    cache,
                    params,
                    policy,
                    directions,
                )
            else
                result[i] = _reducepoly_varbyvar!(
                    runtime,
                    FT[i],
                    S,
                    P,
                    fe,
                    parent(f),
                    cache,
                    params,
                    policy,
                    directions,
                )
            end
        end
    end

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
    return _reducetransform_varbyvar_engine(
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
    )
end

function varbyvar_chooseV(I, d)
    v = zeros(Int, length(I))
    n = length(I)
    J = copy(I)
    i = 0
    while i < d
        for j in eachindex(J)
            if J[n-j+1] > 0
                J[n-j+1] = J[n-j+1] - 1
                v[n-j+1] = v[n-j+1] + 1
                i = i + 1
                break
            elseif j == n
                throw(ArgumentError("varbyvar could not choose a reduction direction"))
            end
        end
    end
    return v
end

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
    chains = AbstractReductionChainHandle[]

    for (_, inputs) in keys_by_order
        for input in inputs
            key = ReductionChainKey(input.key[1], input.key[2])
            push!(chains, materialize(engine, key, input))
        end
    end

    isempty(chains) && return nothing

    first_key = reduction_chain_key(engine, first(chains))
    seen_keys = Set{Tuple{typeof(first_key.u),Int}}()
    duplicate_orders = Set{Int}()
    chain_orders = Int[]
    sizehint!(chain_orders, length(chains))
    next_chains = AbstractReductionChainHandle[]
    sizehint!(next_chains, length(chains))

    while !isempty(chains)
        for chain in chains
            key = reduction_chain_key(engine, chain)
            key.order == terminal_order && continue
            key.order > terminal_order ||
                throw(ArgumentError("varbyvar chain is below the terminal order"))

            direction = varbyvar_chooseV(collect(key.u), policy.degree)
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
