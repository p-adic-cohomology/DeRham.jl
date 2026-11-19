
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

function reducechain_varbyvar(u, g, m, S, f, P, p, context, cache, params, kernel = nothing)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    terminal = terminal_pole_order(P)
    PR = parent(f)
    R = coefficient_ring(parent(f))

    J = copy(u)

    modm = Integer(modulus(base_ring(parent(f))))
    if kernel === nothing
        kernel = _varbyvar_evaluation_kernel(params, d, n, modm)
    end

    gMat = g
    mins = copy(J)
    tempv = copy(J)

    (4 < params.verbose) && println("Starting: J = $J")

    l = 1
    for i in eachindex(J)
        if J[n+1-i+1] > 0
            l = n+1-i+1
            break
        end
    end

    highpole = true
    if m == terminal
        highpole == false
    end
    while highpole && J[l] > 0
        V = varbyvar_chooseV(J, d)

        (4 < params.verbose) && print("Chose V = $V; ")
        (6 < params.verbose) && begin end

        K = 0
        while true
            @. tempv = mins - V
            isLessThanZero = false
            for j in tempv
                if j < 0
                    isLessThanZero = true
                    break
                end
            end
            if isLessThanZero == true
                break
            end
            if m - K == terminal
                break
            end
            @. mins = tempv
            K = K+1
        end

        (4 < params.verbose) && println("Getting Ruv matrices; ")
        matrices = context.Ruvs[V]
        (4 < params.verbose) && println("Computing A and B; ")
        kernel(context.B, context.A, context.temp, matrices, mins, V)
        (4 < params.verbose) && println("Starting the reduction steps; ")
        gMat = finitediff_prodeval_linear!(
            context.B,
            context.A,
            0,
            K-1,
            gMat,
            context.temp,
            context.g_temp,
        )
        @. J = J - K*V
        m = m - K
        if m == terminal
            highpole = false
        end

        (4 < params.verbose) && print("After $(lpad(K,4,' ')) steps,")
        (4 < params.verbose) && println("J = $J")

    end

    return ((J, gMat), m)
end

function reducepoly_varbyvar(pol, S, f, P, fe, p, context, cache, params, kernel = nothing)
    terminal = terminal_pole_order(P)
    PR = parent(f)

    vectype = typeof(context.g)
    allcostadata = Vector{Tuple{Tuple{Vector{Int},vectype},Int}}()
    for (_, inputs) in reduction_input_keys(fe, pol, S)
        for input in inputs
            u, order = input.key
            g = reduction_input_vector!(
                _zero_reduction_vector(context.g),
                fe,
                input,
                S,
                cache,
                params,
            )
            push!(allcostadata, ((u, g), order))
        end
    end

    notallred = true
    # TODO: prime the jitter here?

    while notallred
        for i in eachindex(allcostadata)
            allcostadata[i] = reducechain_varbyvar(
                allcostadata[i][1]...,
                allcostadata[i][2],
                S,
                f,
                P,
                p,
                context,
                cache,
                params,
                kernel,
            )
        end
        if params.use_gpu == true
            remove_duplicates_gpu!(allcostadata)
        else
            remove_duplicates!(allcostadata)
        end
        for i = 1:length(allcostadata)
            if allcostadata[i][2] > terminal
                break
            elseif i == length(allcostadata)
                notallred = false
            end
        end
    end
    terminal_chains = [entry[1] for entry in allcostadata]
    return reduction_terminal_polynomial(fe, terminal_chains, PR, S, params)
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
            directions(key.u, chunk, direction) || throw(
                ArgumentError("varbyvar selected an invalid reduction direction"),
            )
            reduce!(engine, chain, direction, chunk)
        end

        orders = sort!(unique([reduction_chain_key(engine, chain).order for chain in chains]))
        next_chains = AbstractReductionChainHandle[]
        for order in orders
            survivors = combine!(engine, order)
            if order == terminal_order
                for chain in survivors
                    on_terminal(finish!(engine, chain))
                end
            else
                append!(next_chains, survivors)
            end
        end
        chains = next_chains
    end

    return nothing
end

public VarByVarPolicy
