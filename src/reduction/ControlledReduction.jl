
function reducetransform(
    FT,
    N_m,
    S,
    f,
    P,
    fe::ZetaFrontend,
    pseudoInverseMat,
    p,
    params,
    cache,
    context;
    setup = nothing,
)
    setup === nothing && (
        setup = reduction_setup(
            fe,
            S,
            params.algorithm;
            use_gpu = params.use_gpu,
            use_threads = params.use_threads,
            always_use_bigints = params.always_use_bigints,
            fastevaluation = params.fastevaluation,
            context = context,
        )
    )
    if setup.policy == :pchunk
        reducetransform_pchunk(
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
            setup,
        )
    elseif setup.policy == :depthfirst
        reducetransform_depthfirst(
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
            setup,
        )
    elseif setup.policy == :varbyvar
        reducetransform_varbyvar(
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
            setup,
        )
    elseif setup.policy == :akr
        reducetransform_akr(FT, N_m, S, f, pseudoInverseMat, p, cache, params)
    else
        println(setup.policy)
        throw(ArgumentError("Unsupported Algorithm: "))
    end
end

function reducetransform(
    FT,
    N_m,
    S,
    f,
    P::HypersurfaceQuotient,
    pseudoInverseMat,
    p,
    params,
    cache,
    context;
    setup = nothing,
)
    return reducetransform(
        FT,
        N_m,
        S,
        f,
        P,
        ProjectiveHypersurface(f),
        pseudoInverseMat,
        p,
        params,
        cache,
        context;
        setup = setup,
    )
end

function reducetransform(
    FT,
    N_m,
    S,
    f,
    P::DoubleCoverQuotient,
    pseudoInverseMat,
    p,
    params,
    cache,
    context;
    setup = nothing,
)
    return reducetransform(
        FT,
        N_m,
        S,
        f,
        P,
        DoubleCover(f),
        pseudoInverseMat,
        p,
        params,
        cache,
        context;
        setup = setup,
    )
end

function reducetransform(
    FT,
    N_m,
    S,
    f,
    pseudoInverseMat,
    p,
    params,
    cache,
    context;
    setup = nothing,
)
    fe = ProjectiveHypersurface(f)
    P = cohomology_quotient(fe)
    return reducetransform(
        FT,
        N_m,
        S,
        f,
        P,
        fe,
        pseudoInverseMat,
        p,
        params,
        cache,
        context;
        setup = setup,
    )
end

function reducetransform(FT, N_m, S, f, pseudoInverseMat, p, params, cache; setup = nothing)
    return reducetransform(
        FT,
        N_m,
        S,
        f,
        pseudoInverseMat,
        p,
        params,
        cache,
        nothing;
        setup = setup,
    )
end

struct MatrixReductionRuntime{E}
    engine::E
end

function _matrix_reduction_runtime_factory(execution_setup, fe, P, S, cache, params, directions)
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
        return MatrixReductionRuntime(engine)
    end
end

function _reducepoly_matrix_policy!(
    runtime::MatrixReductionRuntime,
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
    on_terminal = chain -> begin
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

function _reducetransform_matrix_policy(
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
    policy,
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
    runtime_factory =
        _matrix_reduction_runtime_factory(execution_setup, fe, P, S, cache, params, directions)
    result = similar(FT)

    reduce_one! = function (runtime, i, threaded)
        if params.verbose > 0
            if threaded
                println("Reducing vector $i in thread $(Threads.threadid())")
            else
                println("Reducing vector $i")
            end
            @time result[i] = _reducepoly_matrix_policy!(
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
            result[i] = _reducepoly_matrix_policy!(
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

    if execution_setup.threads
        runtime_tlv = OhMyThreads.TaskLocalValue{MatrixReductionRuntime}(runtime_factory)
        Threads.@threads for i in eachindex(FT)
            reduce_one!(runtime_tlv[], i, true)
        end
    else
        runtime = runtime_factory()
        for i in eachindex(FT)
            reduce_one!(runtime, i, false)
        end
    end

    return result, execution_setup
end
