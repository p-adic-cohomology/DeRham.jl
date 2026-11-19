
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
    setup = reduction_execution_setup(
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
    Ruv = setup.Ruv
    result = similar(FT)

    if !setup.context_provided
        #TODO: right now, it usually seems better to do lazy computations,
        #  since not all of the Ruv are used. However, I know that for some
        #  classes of examples, they are all pretty much always used. For
        #  such examples, it's better to use an EagerPEP and do threads.
        if setup.threads
            context_tlv =
                OhMyThreads.TaskLocalValue{setup.context_type}(setup.context_factory)
            Threads.@threads for i = 1:length(FT)
                local context = context_tlv[]

                pol = FT[i]
                if (0 < params.verbose)
                    println("Reducing vector $i")
                    @time reduction = reducepoly_varbyvar(
                        pol,
                        S,
                        f,
                        P,
                        fe,
                        p,
                        context,
                        cache,
                        params,
                        setup.kernel,
                    )
                else
                    reduction = reducepoly_varbyvar(
                        pol,
                        S,
                        f,
                        P,
                        fe,
                        p,
                        context,
                        cache,
                        params,
                        setup.kernel,
                    )
                end
                result[i] = reduction

            end
        else

            context = setup.context_factory()
            for i = 1:length(FT) #pol in FT


                pol = FT[i]
                if (0 < params.verbose)
                    println("Reducing vector $i")
                    @time reduction = reducepoly_varbyvar(
                        pol,
                        S,
                        f,
                        P,
                        fe,
                        p,
                        context,
                        cache,
                        params,
                        setup.kernel,
                    )
                else
                    reduction = reducepoly_varbyvar(
                        pol,
                        S,
                        f,
                        P,
                        fe,
                        p,
                        context,
                        cache,
                        params,
                        setup.kernel,
                    )
                end
                result[i] = reduction

            end
        end
    else
        context = setup.context_factory()
        for i = 1:length(FT) #pol in FT
            pol = FT[i]
            if (0 < params.verbose)
                println("Reducing vector $i")
                @time reduction = reducepoly_varbyvar(
                    pol,
                    S,
                    f,
                    P,
                    fe,
                    p,
                    context,
                    cache,
                    params,
                    setup.kernel,
                )
            else
                reduction = reducepoly_varbyvar(
                    pol,
                    S,
                    f,
                    P,
                    fe,
                    p,
                    context,
                    cache,
                    params,
                    setup.kernel,
                )
            end
            result[i] = reduction
        end
    end

    (0 < params.verbose) && begin
        println("Created $(length(allpoints(Ruv))) of $(length(cache[P.d])) possible V")
    end
    (1 < params.verbose) && begin
        println("V that were created: \n$(allpoints(Ruv))")
    end

    return result
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
                throw(error)
            end
        end
    end
    return v
end
