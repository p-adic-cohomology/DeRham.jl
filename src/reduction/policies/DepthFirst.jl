
# what in here could possibly be causing a lock conflict?
"""
    reducechain_depthfirst

Iteratively execute reduction chunks in the navie strategy until
the vector g is reduced to the terminal pole order of P

fields
------
u -  vector of ints
g - vector
m - pole order
S - vector of ints
f - polynomial
P - cohomology quotient
p - prime number
context - ControlledReductionContext
cache - the GradedExpCache used for this controlled reduction
params - the ControlledReductionParamaters
"""
function reducechain_depthfirst(u, g, m, S, f, P, p, context, cache, params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    reduction_degree = _reduction_g_degree(P)
    terminal = terminal_pole_order(P)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    (9 < params.verbose) && println("u = $u")

    J = copy(u)

    gMat = context.g
    mins = similar(J)
    tempv = similar(J)
    (4 < params.verbose) && println("Starting: J = $J")
    (5 < params.verbose) && begin
        g_poly = vector_to_polynomial(g, n, reduction_degree, PR, params.termorder)
        if params.always_use_bigints || params.use_gpu
            println("Starting: g = $((gMat)) = $g_poly")
        else
            println("Starting: g = $(Int.(gMat)) = $g_poly")
        end
    end

    firsttime = true


    while m > terminal
        V = rev_chooseV(J, d, S)

        (4 < params.verbose) && print("Chose V = $V; ")
        (6 < params.verbose) && begin end
        @. mins = J
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
        matrices = context.Ruvs[V]

        (6 < params.verbose && V == [0, 0, 0, 3] && firsttime) && begin
            println(matrices[5][:, 25])
            firsttime=false
            error()
        end

        eval_to_linear!(context.B, context.A, context.temp, matrices, mins, V)

        i = 1
        if params.fastevaluation == false
            params.verbose == 11 && println("gMat before is $gMat")
            while i <= K
                gMat = (context.A+context.B*(K-i))*gMat
                i = i+1
                if params.verbose == 11
                    A = context.A
                    B = context.B
                    g = vector_to_polynomial(
                        gMat,
                        n,
                        reduction_degree,
                        PR,
                        params.termorder,
                    )
                    println("gMat after $i is $gMat = $g")
                end
            end
        else
            gMat = prodeval_linear!(
                context.B,
                context.A,
                0,
                K-1,
                gMat,
                context.temp,
                context.g_temp,
                params,
            )
        end
        @. J = J - K*V
        m = m - K
        (4 < params.verbose) && print("After $(lpad(K,4,' ')) steps,")
        (4 < params.verbose) && println("J = $J")
        if (5 < params.verbose)
            CUDA.@allowscalar g =
                vector_to_polynomial(gMat, n, reduction_degree, PR, params.termorder)
            if params.always_use_bigints || params.use_gpu
                println("g = $((gMat)) = $g")
            elseif params.fastevaluation
                println("g = $(Int.(gMat)) = $g")
            else
                println("g = $(gMat) = $g")
            end
        end

    end
    return (J, gMat)
end

function reducepoly_depthfirst(pol, S, f, P, fe, p, context, cache, params)
    PR = parent(f)
    terminal_chains = (
        reducechain_depthfirst(
            input.key[1],
            reduction_input_vector!(context.g, fe, input, S, cache, params),
            input.key[2],
            S,
            f,
            P,
            p,
            context,
            cache,
            params,
        ) for (_, inputs) in reduction_input_keys(fe, pol, S) for input in inputs
    )
    return reduction_terminal_polynomial(fe, terminal_chains, PR, S, params)
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
)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    g = _reduction_g_degree(P)
    g_length = binomial(g+n, n)

    MS1 = matrix_space(coefficient_ring(parent(f)), g_length, g_length)
    m = Integer(modulus(base_ring(MS1)))

    ruv_function = reduction_ruv_function(fe, S, f, pseudoInverseMat, cache, params)
    computeRuv = if 3 < params.verbose
        V -> begin
            println("Computing Ruv for V = $V for the first time.")
            @time ruv_function(V)
        end
    else
        ruv_function
    end

    result = similar(FT)

    if context == nothing
        #TODO: right now, it usually seems better to do lazy computations,
        #  since not all of the Ruv are used. However, I know that for some
        #  classes of examples, they are all pretty much always used. For
        #  such examples, it's better to use an EagerPEP and do threads.
        lazy_Ruv = length(S) < d || d < n

        if (0 < params.verbose)
            println("Creating the Ruv PEP object...")
            @time Ruv = select_Ruv_PEP(n, d, S, params, computeRuv, lazy_Ruv, MS1, cache)
        else
            Ruv = select_Ruv_PEP(n, d, S, params, computeRuv, lazy_Ruv, MS1, cache)
        end


        if params.use_threads

            context_tlv = OhMyThreads.TaskLocalValue{default_context_type(MS1, params)}(
                () -> default_context(MS1, Ruv, params),
            )

            Threads.@threads for i = 1:length(FT)
                local context = context_tlv[]

                pol = FT[i]
                if (0 < params.verbose)
                    println("Reducing vector $i in thread $(Threads.threadid())")
                    @time reduction =
                        reducepoly_depthfirst(pol, S, f, P, fe, p, context, cache, params)
                else
                    reduction =
                        reducepoly_depthfirst(pol, S, f, P, fe, p, context, cache, params)
                end
                result[i] = reduction

            end
        else
            context = default_context(MS1, Ruv, params)
            for i = 1:length(FT) #pol in FT

                pol = FT[i]
                if (0 < params.verbose)
                    println("Reducing vector $i")
                    @time reduction =
                        reducepoly_depthfirst(pol, S, f, P, fe, p, context, cache, params)
                else
                    reduction =
                        reducepoly_depthfirst(pol, S, f, P, fe, p, context, cache, params)
                end
                result[i] = reduction

            end
        end
    else
        if (ZZ(2)^25 < m < ZZ(2)^106)
            compute_gpu = V -> KaratsubaMat.(computeRuv(V))
        else
            compute_gpu = V -> cuMod.(computeRuv(V))
        end

        compute_float = V -> float_entries.(computeRuv(V))

        if 3 < n && d == 3 && S == [n] && params.use_gpu
            context.Ruvs.backing.compute = compute_float
        elseif params.use_gpu
            context.Ruvs.compute = compute_gpu
        else
            context.Ruvs.compute = computeRuv
        end

        for i = 1:length(FT) #pol in FT


            pol = FT[i]
            if (0 < params.verbose)
                println("Reducing vector $i in thread $(Threads.threadid())")
                @time reduction =
                    reducepoly_depthfirst(pol, S, f, P, fe, p, context, cache, params)
            else
                reduction =
                    reducepoly_depthfirst(pol, S, f, P, fe, p, context, cache, params)
            end
            result[i] = reduction
        end
    end

    (0 < params.verbose) && begin
        println("Created $(length(allpoints(Ruv))) of $(length(cache[d])) possible V")
    end
    (1 < params.verbose) && begin
        println("V that were created: \n$(allpoints(Ruv))")
    end

    return result
end
