
function _pchunk_choose_direction(key::ReductionChainKey, degree, S)
    exponents = collect(key.u)
    sum(exponents) >= degree ||
        throw(ArgumentError("pchunk has no direction for the chain key"))
    return rev_chooseV(exponents, degree, S)
end

"""
    reducechain_pchunk(key, g, S, f, P, pseudoInverseMat, p, Ruv, cache, A, B, temp, g_temp, params)

takes a single monomial in Frobenius and reduces one chunk toward
`terminal_pole_order(P)`.


if the reduction hits the end, returns u as the "true" value, otherwise returns it in Costa's format
(i.e. entries will be multiplies of p in Costa's format)

fields
------
key - chain key containing the exponent vector and current pole order
g - vector
S - vector of ints
f - polynomial
P - cohomology quotient
pseudoInverseMat - output of pseudo_inverse_controlled_lifted
p - prime number
Ruv - output of computeRuvS
cache - the GradedExpCache used for this controlled reduction
A - matrix
B - matrix
temp - preallocated storage matrix
g_temp - preallocated storage vector
params - the ControlledReductionParamaters
"""
function reducechain_pchunk(
    key::ReductionChainKey,
    g,
    S,
    f,
    P,
    pseudoInverseMat,
    p,
    Ruv,
    cache,
    A,
    B,
    temp,
    g_temp,
    params,
    directions,
)
    verbose = params.verbose
    d = total_degree(f)
    terminal = terminal_pole_order(P)

    (4 < verbose) && println("Expanded U: $(key.u)")

    gMat = g
    (4 < verbose) && println("Before reduction chunk, U is $(key.u)")
    maximum_chunk = min(key.order - terminal, p)
    maximum_chunk <= 0 && return (key, gMat)

    direction = _pchunk_choose_direction(key, d, S)
    chunk = _depthfirst_chunk_length(collect(key.u), direction, maximum_chunk)
    chunk >= 1 || throw(ArgumentError("pchunk selected an empty chunk"))
    directions(key.u, chunk, direction) ||
        throw(ArgumentError("pchunk selected an invalid reduction direction"))
    (4 < verbose) && println("Reducing with V = $direction for $chunk steps")

    matrices = Ruv[direction]
    next_u = collect(key.u) .- chunk .* direction
    B, A = eval_to_linear!(B, A, temp, matrices, next_u, direction)

    if params.fastevaluation
        gMat = finitediff_prodeval_linear!(B, A, 0, chunk - 1, gMat, temp, g_temp)
    else
        for i = 1:chunk
            my_mul!(temp, B, chunk - i)
            my_add!(temp, temp, A)
            my_matvecmul!(g_temp, temp, gMat)
            my_copy!(gMat, g_temp)
        end
    end

    return (ReductionChainKey(next_u, key.order - chunk), gMat)
end

"""
    reducepoly_pchunk(pol,S,f,pseudoInverseMat,p,Ruvs,termorder)

Implements Costa's algorithm for controlled reduction,
reducing each chain in chunks of at most `p` steps.

fields
------
pol - polynomial
S - vector of ints
f - polynomial
pseudoInverseMat - output of pseudo_inverse_controlled_lifted
p - prime number
Ruv - output of computeRuvS
cache - the GradedExpCache used for this controlled reduction
A - matrix
B - matrix
temp - preallocated storage matrix
params - the ControlledReductionParamaters
"""
function reducepoly_pchunk(
    pol,
    S,
    f,
    P,
    pseudoInverseMat,
    p,
    Ruv,
    cache,
    A,
    B,
    temp,
    params,
    directions,
)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    reduction_degree = _reduction_g_degree(P)
    g_length = binomial(_reduction_g_degree(P)+n, n)
    terminal = terminal_pole_order(P)

    allcostadata = []
    key_indices = Dict{ReductionChainKey,Int}()
    g_temp = zeros(UInt, g_length)

    for polpol in pol
        for term_with_order in wholeterms(polpol)
            order = term_with_order[2]
            order < terminal && continue
            g = zeros(UInt, g_length)
            term_costadata =
                costadata_of_initial_term!(term_with_order, g, n, d, P, p, S, cache, params)
            key = ReductionChainKey(tweak(term_costadata[1], reduction_degree), order)
            existing_index = get(key_indices, key, nothing)
            if existing_index === nothing
                push!(allcostadata, (key, term_costadata[2]))
                key_indices[key] = length(allcostadata)
            else
                allcostadata[existing_index] =
                    (key, allcostadata[existing_index][2] .+ term_costadata[2])
            end
        end
    end

    for i in eachindex(allcostadata)
        key, gMat = allcostadata[i]
        while key.order > terminal
            key, gMat = reducechain_pchunk(
                key,
                gMat,
                S,
                f,
                P,
                pseudoInverseMat,
                p,
                Ruv,
                cache,
                A,
                B,
                temp,
                g_temp,
                params,
                directions,
            )
        end
        allcostadata[i] = (collect(key.u), gMat)
    end

    return poly_of_end_costadatas(allcostadata, PR, P, p, d, n, S, params)
end

"""
    reducetransform_pchunk(FT,N_m,S,f,pseudoInverseMat,p,cache,params)

trying to emulate Costa's controlled reduction, changes the order that polynomials are reduced, starts from highest pole order and accumulates the lower order poles as reduction proceeds

fields
------
FT - applyFrobeniusToBasis
N_m - the precision
S - vector of ints
f - polynomial
pseudoInverseMat - output of pseudo_inverse_controlled_lifted
p - prime number
cache - the GradedExpCache used for this controlled reduction
params - the ControlledReductionParamaters
"""
function reducetransform_pchunk(FT, N_m, S, f, P, pseudoInverseMat, p, cache, params, fe)

    d = total_degree(f)
    n = nvars(parent(f)) - 1
    g_length = binomial(_reduction_g_degree(P)+n, n)
    MS1 = matrix_space(coefficient_ring(parent(f)), g_length, g_length)
    A = MS1()
    B = MS1()
    temp = MS1()
    ruv_function = reduction_ruv_function(fe, S, f, pseudoInverseMat, cache, params)
    computeRuv = if 3 < params.verbose
        V -> begin
            println("Computing Ruv for V = $V for the first time.")
            @time ruv_function(V)
        end
    else
        ruv_function
    end
    Ruv = LazyPEP{typeof(MS1())}(computeRuv)
    directions = ReductionDirectionSet(reduction_directions(fe, S))

    result = []
    i = 1
    for pol in FT
        (0 < params.verbose) && println("Reducing vector $i")
        i += 1
        if (0 < params.verbose)
            @time reduction = reducepoly_pchunk(
                pol,
                S,
                f,
                P,
                pseudoInverseMat,
                p,
                Ruv,
                cache,
                A,
                B,
                temp,
                params,
                directions,
            )
        else
            reduction = reducepoly_pchunk(
                pol,
                S,
                f,
                P,
                pseudoInverseMat,
                p,
                Ruv,
                cache,
                A,
                B,
                temp,
                params,
                directions,
            )
        end

        push!(result, reduction)
    end

    return result
end
