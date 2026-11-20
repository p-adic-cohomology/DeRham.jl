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

struct PChunkEngineState{E}
    engine::E
end

function _pchunk_engine_factory(execution_setup, fe, P, S, cache, params, directions)
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
        return PChunkEngineState(engine)
    end
end

function _reducepoly_pchunk!(
    runtime::PChunkEngineState,
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
    policy = PChunkPolicy(P.d, S, p)
    engine_factory =
        _pchunk_engine_factory(execution_setup, fe, P, S, cache, params, directions)
    result = similar(FT)

    if execution_setup.threads
        runtime_tlv = OhMyThreads.TaskLocalValue{PChunkEngineState}(engine_factory)
        Threads.@threads for i in eachindex(FT)
            local runtime = runtime_tlv[]
            if 0 < params.verbose
                println("Reducing vector $i in thread $(Threads.threadid())")
                @time result[i] = _reducepoly_pchunk!(
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
                result[i] = _reducepoly_pchunk!(
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
                @time result[i] = _reducepoly_pchunk!(
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
                result[i] = _reducepoly_pchunk!(
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

    return result
end

public PChunkPolicy
