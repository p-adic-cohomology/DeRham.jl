"""
    AKRNaiveDenseReduction()

Experimental dense reduction that applies each AKR pole-reduction chain
independently.
"""
struct AKRNaiveDenseReduction <: AbstractDenseReduction end

"""
    AKRLevelByLevelDenseReduction()

Experimental dense reduction that combines terms at each pole order before
applying the next AKR reduction step.
"""
struct AKRLevelByLevelDenseReduction <: AbstractDenseReduction end

function _reducechain_akr_naive(g, m, f, p, picache, params; target_order = nothing)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    target_order === nothing && (target_order = n)
    target_order >= n ||
        throw(ArgumentError("target order is below the terminal pole order"))

    while m > target_order
        gtemp = picache[m]*g
        len = Int(length(gtemp)/(n+1))
        g = derivative(
            vector_to_polynomial(gtemp[1:len], n, d*(m-1)-n, PR, params.termorder),
            1,
        )
        for i = 1:n
            gpolytemp = derivative(
                vector_to_polynomial(
                    gtemp[(i*len+1):((i+1)*len)],
                    n,
                    d*(m-1)-n,
                    PR,
                    params.termorder,
                ),
                i+1,
            )
            g = g + gpolytemp
        end
        g = polynomial_to_vector(g, n+1, params.termorder)
        if size(g)[1] == 0
            g = fill(R(0), binomial(d*(m-1)-1, n))
        end
        m = m - 1
    end

    return g
end

function _reducepoly_akr_naive(pol, S, f, p, picache, params; target_order = nothing)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    target_order === nothing && (target_order = n)
    target_order >= n ||
        throw(ArgumentError("target order is below the terminal pole order"))
    result_by_order = Dict{Int,typeof(PR())}()
    input_orders = sort!(unique(Int(term[2]) for term in pol))

    for order in input_orders
        for term in termsoforder(pol, order)
            if order > target_order
                g = polynomial_to_vector(term[1], n + 1, params.termorder)
                reduced = _reducechain_akr_naive(
                    g,
                    order,
                    f,
                    p,
                    picache,
                    params;
                    target_order = target_order,
                )
                reduced_poly = vector_to_polynomial(
                    reduced,
                    n,
                    d * target_order - n - 1,
                    PR,
                    params.termorder,
                )
                result = get!(result_by_order, target_order) do
                    PR()
                end
                result_by_order[target_order] = result + reduced_poly
            else
                result = get!(result_by_order, order) do
                    PR()
                end
                result_by_order[order] = result + term[1]
            end
        end
    end

    isempty(result_by_order) && return [[PR(), target_order]]
    return [
        [result_by_order[order], order] for order in sort!(collect(keys(result_by_order)))
    ]
end

function _build_akr_picache(
    highm,
    n,
    d,
    S,
    f,
    pseudoInverseMat,
    cache,
    params,
    p;
    target_order = n,
)
    target_order >= n ||
        throw(ArgumentError("target order is below the terminal pole order"))
    M = factor(modulus(base_ring(parent(f))))[p]
    picache = Vector{typeof(pseudoInverseMat)}(undef, highm)
    picache[n] = pseudoInverseMat
    for i = max(n+1, target_order+1):highm
        l = d*i - n - 1
        generate_degree_forward(cache, l)
        generate_degree_forward(cache, l-(d-1))
        generate_degree_forward(cache, l-d)
        pi_new = pseudo_inverse_controlled_lifted(f, S, l, M, params, cache)
        MS = matrix_space(base_ring(parent(f)), nrows(pi_new), ncols(pi_new))
        pseudo_inverse_mat = MS()
        for i = 1:nrows(pi_new)
            for j = 1:ncols(pi_new)
                pseudo_inverse_mat[i, j] = ZZ(pi_new[i, j])
            end
        end
        picache[i] = pseudo_inverse_mat
    end
    return picache
end

function _highest_akr_pole_order(FT)
    return FT[length(FT)][length(FT[length(FT)])][2]
end

function reducetransform_akr_naive(
    FT,
    N_m,
    S,
    f,
    pseudoInverseMat,
    p,
    cache,
    params;
    target_order = nvars(parent(f)) - 1,
)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    highm = _highest_akr_pole_order(FT)
    picache = _build_akr_picache(
        highm,
        n,
        d,
        S,
        f,
        pseudoInverseMat,
        cache,
        params,
        p;
        target_order = target_order,
    )
    reduction = AKRNaiveDenseReduction()
    context = (; S, f, p, picache, params, terminal_order = n)

    result = similar(FT)
    for i = 1:length(FT)
        pol = FT[i]
        if 0 < params.verbose
            println("Reducing vector $i")
            @time reduced = dense_reduce(reduction, pol, highm, target_order, context)
        else
            reduced = dense_reduce(reduction, pol, highm, target_order, context)
        end
        result[i] = reduced
    end

    return result
end

"""
    _reduce_akr_level(g, m, f, picache, params)

Apply one AKR reduction step to the combined coefficient vector at pole order
`m`, returning the coefficient vector for pole order `m - 1`.
"""
function _reduce_akr_level(g, m, f, picache, params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    gtemp = picache[m]*g
    len = Int(length(gtemp)/(n+1))
    reduced = derivative(
        vector_to_polynomial(gtemp[1:len], n, d*(m-1)-n, PR, params.termorder),
        1,
    )
    for i = 1:n
        component = vector_to_polynomial(
            gtemp[(i*len+1):((i+1)*len)],
            n,
            d*(m-1)-n,
            PR,
            params.termorder,
        )
        reduced += derivative(component, i+1)
    end
    reduced_vector = polynomial_to_vector(reduced, n+1, params.termorder)
    if size(reduced_vector)[1] == 0
        reduced_vector = fill(R(0), binomial(d*(m-1)-1, n))
    end
    return reduced_vector
end

"""
    _reducepoly_akr_levelbylevel(pol, S, f, p, picache, params)

Reduce all terms at each pole order together. The derivative still uses the
existing polynomial-vector conversions used by the naive AKR implementation.
"""
function _reducepoly_akr_levelbylevel(pol, S, f, p, picache, params; target_order = nothing)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    target_order === nothing && (target_order = n)
    target_order >= n ||
        throw(ArgumentError("target order is below the terminal pole order"))
    result_by_order = Dict{Int,typeof(PR())}()
    input_orders = sort!(unique(Int(term[2]) for term in pol))
    isempty(input_orders) && return [[PR(), target_order]]
    for order in input_orders
        order <= target_order || continue
        for term in termsoforder(pol, order)
            result = get!(result_by_order, order) do
                PR()
            end
            result_by_order[order] = result + term[1]
        end
    end

    highm = maximum(input_orders)
    if highm <= target_order
        return [
            [result_by_order[order], order] for
            order in sort!(collect(keys(result_by_order)))
        ]
    end
    degree = d*highm - n - 1
    g = fill(R(0), binomial(degree+n, n))

    for m = highm:-1:(target_order+1)
        for term in termsoforder(pol, m)
            g .+= polynomial_to_vector(term[1], n+1, params.termorder)
        end
        g = _reduce_akr_level(g, m, f, picache, params)
    end

    result = vector_to_polynomial(g, n, d * target_order - n - 1, PR, params.termorder)
    target_result = get!(result_by_order, target_order) do
        PR()
    end
    result_by_order[target_order] = target_result + result
    return [
        [result_by_order[order], order] for order in sort!(collect(keys(result_by_order)))
    ]
end

function reducetransform_akr_levelbylevel(
    FT,
    N_m,
    S,
    f,
    pseudoInverseMat,
    p,
    cache,
    params;
    target_order = nvars(parent(f)) - 1,
)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    highm = _highest_akr_pole_order(FT)
    picache = _build_akr_picache(
        highm,
        n,
        d,
        S,
        f,
        pseudoInverseMat,
        cache,
        params,
        p;
        target_order = target_order,
    )
    reduction = AKRLevelByLevelDenseReduction()
    context = (; S, f, p, picache, params, terminal_order = n)

    result = similar(FT)
    for i = 1:length(FT)
        pol = FT[i]
        if 0 < params.verbose
            println("Reducing vector $i")
            @time reduced = dense_reduce(reduction, pol, highm, target_order, context)
        else
            reduced = dense_reduce(reduction, pol, highm, target_order, context)
        end
        result[i] = reduced
    end

    return result
end

function _validate_akr_dense_reduction(pol, m_src, m_tgt, context)
    m_tgt >= context.terminal_order ||
        throw(ArgumentError("target order is below the terminal pole order"))
    m_src >= m_tgt || throw(ArgumentError("source order must be at least the target order"))
    all(term -> Int(term[2]) <= m_src, pol) ||
        throw(ArgumentError("input contains a term above the source pole order"))
    return nothing
end

function dense_reduce(
    reduction::AKRNaiveDenseReduction,
    pol,
    m_src::Integer,
    m_tgt::Integer,
    context::NamedTuple,
)
    _validate_akr_dense_reduction(pol, m_src, m_tgt, context)
    return _reducepoly_akr_naive(
        pol,
        context.S,
        context.f,
        context.p,
        context.picache,
        context.params;
        target_order = Int(m_tgt),
    )
end

function dense_reduce(
    reduction::AKRLevelByLevelDenseReduction,
    pol,
    m_src::Integer,
    m_tgt::Integer,
    context::NamedTuple,
)
    _validate_akr_dense_reduction(pol, m_src, m_tgt, context)
    return _reducepoly_akr_levelbylevel(
        pol,
        context.S,
        context.f,
        context.p,
        context.picache,
        context.params;
        target_order = Int(m_tgt),
    )
end
