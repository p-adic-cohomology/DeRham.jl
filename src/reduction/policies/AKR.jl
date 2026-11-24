function _reducechain_akr_naive(g, m, f, p, picache, params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))

    while m > n
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

function _reducepoly_akr_naive(pol, S, f, p, picache, params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    result = PR()
    for term in pol
        terms = termsoforder(pol, term[2])
        for t in terms
            g = polynomial_to_vector(t[1], n+1, params.termorder)
            reduced = _reducechain_akr_naive(g, t[2], f, p, picache, params)
            reduced_poly = vector_to_polynomial(reduced, n, d*n-n-1, PR, params.termorder)
            result += reduced_poly
        end
    end

    return [[result, n]]
end

function _build_akr_picache(highm, n, d, S, f, pseudoInverseMat, cache, params, p)
    M = factor(modulus(base_ring(parent(f))))[p]
    picache = Vector{typeof(pseudoInverseMat)}(undef, highm)
    picache[n] = pseudoInverseMat
    for i = (n+1):highm
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

function reducetransform_akr_naive(FT, N_m, S, f, pseudoInverseMat, p, cache, params)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    highm = _highest_akr_pole_order(FT)
    picache = _build_akr_picache(highm, n, d, S, f, pseudoInverseMat, cache, params, p)

    result = similar(FT)
    for i = 1:length(FT)
        pol = FT[i]
        if 0 < params.verbose
            println("Reducing vector $i")
            @time reduction = _reducepoly_akr_naive(pol, S, f, p, picache, params)
        else
            reduction = _reducepoly_akr_naive(pol, S, f, p, picache, params)
        end
        result[i] = reduction
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
function _reducepoly_akr_levelbylevel(pol, S, f, p, picache, params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    isempty(pol) && return [[PR(), n]]
    highm = maximum(Int(term[2]) for term in pol)
    degree = d*highm - n - 1
    g = fill(R(0), binomial(degree+n, n))

    for m = highm:-1:(n+1)
        for term in termsoforder(pol, m)
            g .+= polynomial_to_vector(term[1], n+1, params.termorder)
        end
        g = _reduce_akr_level(g, m, f, picache, params)
    end

    result = vector_to_polynomial(g, n, d*n-n-1, PR, params.termorder)
    return [[result, n]]
end

function reducetransform_akr_levelbylevel(FT, N_m, S, f, pseudoInverseMat, p, cache, params)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    highm = _highest_akr_pole_order(FT)
    picache = _build_akr_picache(highm, n, d, S, f, pseudoInverseMat, cache, params, p)

    result = similar(FT)
    for i = 1:length(FT)
        pol = FT[i]
        if 0 < params.verbose
            println("Reducing vector $i")
            @time reduction = _reducepoly_akr_levelbylevel(pol, S, f, p, picache, params)
        else
            reduction = _reducepoly_akr_levelbylevel(pol, S, f, p, picache, params)
        end
        result[i] = reduction
    end

    return result
end
