struct DepthFirstEngineState{E}
    engine::E
end

function _depthfirst_engine_factory(execution_setup, fe, P, S, cache, params, directions)
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
        return DepthFirstEngineState(engine)
    end
end

function _reducepoly_depthfirst!(
    runtime::DepthFirstEngineState,
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

function reducetransform_depthfirst(
    FT,
    N_m,
    S,
    f,
    P,
    fe,
    pseudoInverseMat,
    p,
    cache,
    params,
    context,
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
    policy = DepthFirstPolicy(P.d, S)
    engine_factory =
        _depthfirst_engine_factory(execution_setup, fe, P, S, cache, params, directions)
    result = similar(FT)

    if execution_setup.threads
        runtime_tlv = OhMyThreads.TaskLocalValue{DepthFirstEngineState}(engine_factory)
        Threads.@threads for i = 1:length(FT)
            local runtime = runtime_tlv[]
            if 0 < params.verbose
                println("Reducing vector $i in thread $(Threads.threadid())")
                @time result[i] = _reducepoly_depthfirst!(
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
                result[i] = _reducepoly_depthfirst!(
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
                @time result[i] = _reducepoly_depthfirst!(
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
                result[i] = _reducepoly_depthfirst!(
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

    (0 < params.verbose) && begin
        println(
            "Created $(length(allpoints(execution_setup.Ruv))) of $(length(cache[P.d])) possible V",
        )
    end
    (1 < params.verbose) && begin
        println("V that were created: \n$(allpoints(execution_setup.Ruv))")
    end

    return result
end
