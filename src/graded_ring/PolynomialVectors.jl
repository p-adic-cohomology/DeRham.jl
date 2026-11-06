
"""
    polynomial_to_vector!(v, f, n, order=:lex, cache=nothing; vars_reversed=false)

Convert (homogeneous) polynomial to vector form with specified order. Default is lexicographic.
Update v with this vector. Assumes v has the correct length.

v - vector
f - an oscar polynomial
n - the number of variables (minus 1???)
order - a symbol which denotes the term order
cache - a GradedExpCache which gives us the variables
vars_reversed - reverses the order of basis vectors at various places
>>>if you don't know what this is, ignore it.
"""
function polynomial_to_vector!(
    v,
    f,
    n,
    order = :lex,
    cache = nothing;
    vars_reversed = false,
)

    #TODO: remove vars_reversed and set it based on the cache's property

    #TODO: if f turns out to be zero, we don't know what the degree should be.
    #
    #How best to fix this?
    R = base_ring(parent(f))
    d = total_degree(f)

    if cache != nothing
        mon_vec = cache[d]
    else
        mon_vec = gen_exp_vec(n, d, order, vars_reversed = vars_reversed)
    end
    res = v
    for i in eachindex(mon_vec)
        if vars_reversed
            reverse!(mon_vec[i])
            res[i] = coeff(f, mon_vec[i])
            reverse!(mon_vec[i])
        else
            res[i] = coeff(f, mon_vec[i])
        end
    end

    res
end

"""
    polynomial_to_vector(f, n, order=:lex, cache=nothing; vars_reversed=false)

Convert (homogeneous) polynomial to vector form with specified order. Default is lexicographic.

f - an oscar polynomial
n - the number of variables (minus 1???)
order - a symbol which denotes the term order
cache - a GradedExpCache which gives us the variables
vars_reversed - should we reverse the variables?
"""
function polynomial_to_vector(f, n, order = :lex, cache = nothing; vars_reversed = false)

    #TODO: remove vars_reversed and set it based on the cache's property

    #TODO: if f turns out to be zero, we don't know what the degree should be.
    #
    #How best to fix this?
    R = base_ring(parent(f))
    d = total_degree(f)

    if cache != nothing
        mon_vec = cache[d]
    else
        mon_vec = gen_exp_vec(n, d, order, vars_reversed = vars_reversed)
    end
    res = fill(R(0), length(mon_vec))
    for i in eachindex(mon_vec)
        if vars_reversed
            reverse!(mon_vec[i])
            res[i] = coeff(f, mon_vec[i])
            reverse!(mon_vec[i])
        else
            res[i] = coeff(f, mon_vec[i])
        end
    end

    res
end

#
"""
vector_to_polynomial(vect, n, d, PR, order=:lex)

Convert vector to polynomial with specified order. Default is lexicographic.

vect - vector of coefficients
n - number of variables-1
d - homogeneous degree
PR - polynomial ring to be the parent of the return value
order - a symbol which denotes the term order
vars_reversed - reverses the order of basis vectors at various places
>>>if you don't know what this is, ignore it.
"""
function vector_to_polynomial(vect, n, d, PR, order = :lex, vars_reversed = false)

    C = MPolyBuildCtx(PR)
    R = base_ring(PR)
    exp_vecs = gen_exp_vec(n+1, d, order)
    @assert length(vect) == length(exp_vecs) "vector has incorrect length for the specified degree"
    for i in eachindex(vect)
        v = vect[i]
        if v isa AbstractFloat
            v = Integer(vect[i])
        end

        if vars_reversed
            push_term!(C, R(v), reverse(exp_vecs[i]))
        else
            push_term!(C, R(v), exp_vecs[i])
        end
    end
    res = finish(C)

    res
end

# Computes the relations
function compute_relations(monomials, partials)
    result = []
    for i in axes(monomials, 1)
        for j in axes(partials, 1)
            push!(result, monomials[i]*partials[j])
        end
    end
    result
end

# Converts vector of homogeneous polynomials to a matrix of their coefficents
function convert_p_to_m(polys, expvecs; vars_reversed = false)
    R = coefficient_ring(parent(polys[1]))
    MS = matrix_space(R, length(polys), length(expvecs))
    result = MS()
    for i in axes(polys, 1)
        for j in axes(expvecs, 1)
            exp_vec = expvecs[j]

            if vars_reversed
                reverse!(exp_vec)
            end

            result[i, j] = coeff(polys[i], exp_vec)

            if vars_reversed
                reverse!(exp_vec)
            end
        end
    end
    return result
end

# Converts Matrix of coefficents to vector of polynomials, each row is one polynomial
function convert_m_to_p(mat, expvec, R, PR)
    result = []
    B = MPolyBuildCtx(PR)
    for i in axes(mat, 1)
        for j in axes(expvec, 1)
            push_term!(B, mat[i, j], expvec[j])
        end
        push!(result, finish(B))
    end
    result
end

"""
    liftCoefficients(R, PR, f)

Lifts the coefficeints of f to the ring R.

Works by lifting coefficients to ZZ and then
converting elements of ZZ to elements of R.
Thus, this method is mostly useful for
lifting something mod p^m to p^n,
for m < n.

INPUTS:
* "f" -- the polynomial to be lifted
* "R" -- the ring for the coefficients to end up in
* "PR" -- the polynomial ring (over R) for the result to end up in
"""
function liftCoefficients(R, PR, f, positiveLift = true)
    t = terms(f)
    sum = 0
    for i in t
        ev = exponent_vector(i, 1)
        c = coeff(i, 1)
        B = MPolyBuildCtx(PR)
        charBaseField = characteristic(parent(f))
        if positiveLift && (lift(ZZ, c) > div(charBaseField, 2))
            push_term!(B, R(lift(ZZ, c)-charBaseField), ev)
        else
            push_term!(B, R(lift(ZZ, c)), ev)
        end
        sum = sum + finish(B)
    end
    return sum
end
