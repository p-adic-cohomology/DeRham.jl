function _graded_monomial_order(order)
    order === :lex && return (false, false)
    order === :neglex && return (false, true)
    order === :invlex && return (true, false)
    throw(ArgumentError("Unsupported order '$order'"))
end

_graded_variable(i, n, reverse_variables) = reverse_variables ? n + 1 - i : i

"""
Find the indexed exponent vector in `gen_exp_vec` order. Lexicographic
indexing skips whole blocks of monomials at each variable.
"""
function _graded_exponent_at!(exponent, degree, index, ordering)
    reverse_variables, descending = ordering
    N = length(exponent)
    if descending
        index = binomial(degree + N - 1, N - 1) + 1 - index
    end
    remaining = degree
    for i = 1:(N-1)
        trailing = N - i
        count = binomial(remaining + trailing, trailing)
        target = count - index
        lo, hi = 0, remaining
        while lo < hi
            mid = (lo + hi) ÷ 2
            if binomial(mid + trailing, trailing) > target
                hi = mid
            else
                lo = mid + 1
            end
        end
        exponent[_graded_variable(i, N, reverse_variables)] = remaining - lo
        index -= count - binomial(lo + trailing, trailing)
        remaining = lo
    end
    exponent[_graded_variable(N, N, reverse_variables)] = remaining
    return exponent
end

function _graded_monomial_index(exponent, degree, ordering)
    reverse_variables, descending = ordering
    N = length(exponent)
    remaining = degree
    index = 1
    for i = 1:(N-1)
        trailing = N - i
        next = remaining - exponent[_graded_variable(i, N, reverse_variables)]
        index +=
            binomial(remaining + trailing, trailing) - binomial(next + trailing, trailing)
        remaining = next
    end
    return descending ? binomial(degree + N - 1, N - 1) + 1 - index : index
end

"""
   gen_exp_vec(n, d, order; vars_reversed=false)

Returns all nonnegative integer lists of length n who entires sum to d

These are the exponent vectors for all the homogeneous monomials of
degree d, in n variables.

TODO: give this function @memoize. perhaps for some big examples the
storage required to store the result isn't worth it. However, for
tests where we're running many similar examples of medium size,
I think it'll be nicer.

INPUTS:
* "n" -- integer
* "d" -- integer
* "order" -- string, monomial ordering, defaulted to lexicographic ordering. Also supports neglex
"""
function gen_exp_vec(n, d, order::Symbol = :lex; vars_reversed = false)
    reverse_variables, descending = _graded_monomial_order(order)
    result = Vector{Vector{Int64}}(undef, binomial(n+d-1, d))
    for i = 1:binomial(n+d-1, d)
        result[i] = zeros(Int64, n)
    end
    if !descending
        for i = 1:n
            variable = _graded_variable(i, n, reverse_variables)
            dtemp = copy(d)
            k = 0
            while k <= (length(result) - 1)
                if i > 1
                    dtemp = copy(d)
                    for j = 1:n
                        dtemp = dtemp - result[length(result)-k][j]
                    end
                end
                if i == n && dtemp > 0
                    result[length(result)-k][variable] = dtemp
                    k = k + 1
                    continue
                end
                if dtemp == 0
                    k = k + 1
                    continue
                end
                if dtemp == 1 && i > 1
                    for j = i:n
                        result[length(result)-(k+j-i)][_graded_variable(
                            j,
                            n,
                            reverse_variables,
                        )] = 1
                    end
                    k = k + n - i + 1
                    continue
                end
                dtemp2 = copy(d - dtemp)
                if i == 1 || dtemp2 == 0
                    while dtemp >= 0
                        for j = 1:binomial(n-i+d-dtemp-1, d-dtemp)
                            result[length(result)-(j+k-1)][variable] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-1, d-dtemp)
                        dtemp = dtemp - 1
                    end
                else
                    while dtemp >= 0
                        for j = 1:binomial(n-i+d-dtemp-dtemp2-1, d-dtemp-dtemp2)
                            result[length(result)-(j+k-1)][variable] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-dtemp2-1, d-dtemp-dtemp2)
                        dtemp = dtemp - 1
                    end
                end
            end
        end
    else
        for i = 1:n
            dtemp = copy(d)
            k = 1
            while k <= length(result)
                if i > 1
                    dtemp = copy(d)
                    for j = 1:n
                        dtemp = dtemp - result[k][j]
                    end
                end
                if i == n && dtemp > 0
                    result[k][i] = dtemp
                    k = k + 1
                    continue
                end
                if dtemp == 0
                    k = k + 1
                    continue
                end
                if dtemp == 1 && i > 1
                    for j = i:n
                        result[k+j-i][j] = 1
                    end
                    k = k + n - i + 1
                    continue
                end
                dtemp2 = copy(d - dtemp)
                if i == 1 || dtemp2 == 0
                    while dtemp >= 0
                        for j = 1:binomial(n-i+d-dtemp-1, d-dtemp)
                            result[j+k-1][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-1, d-dtemp)
                        dtemp = dtemp - 1
                    end
                else
                    while dtemp >= 0
                        for j = 1:binomial(n-i+d-dtemp-dtemp2-1, d-dtemp-dtemp2)
                            result[j+k-1][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-dtemp2-1, d-dtemp-dtemp2)
                        dtemp = dtemp - 1
                    end
                end
            end
        end
    end

    if vars_reversed
        for i = 1:binomial(n+d-1, d)
            reverse!(result[i])
        end
    end
    return result
end

"""
    gen_mon

fields
------
exp_vecs - output of gen_mon
PR - polynomial ring
"""
function gen_mon(exp_vecs, PR)
    result = zeros(PR, length(exp_vecs))
    R = base_ring(PR)
    B = MPolyBuildCtx(PR)
    for i in axes(exp_vecs, 1)
        push_term!(B, one(R), exp_vecs[i])
        result[i] = finish(B)
    end
    result
end

"""
    compute_monomials

Computes the monomials in n variables, of degree d, in the
variable order order, as Oscar polynonmials.

fields
------
n - int
d - int
PR - polynomial ring
order - the term ordering that should be used in vector representations
cache - the GradedExpCache used for this controlled reduction
vars_reversed - reverses the order of basis vectors at various places
>>>if you don't know what this is, ignore it.
"""
function compute_monomials(n, d, PR, order = :lex, cache = nothing; vars_reversed = false)
    if n < 0 || d < 0
        return []
    end
    if cache != nothing
        exp_vecs = cache[d]
    else
        exp_vecs = gen_exp_vec(n, d, order)
    end

    if vars_reversed
        reverse!.(exp_vecs)
        res = gen_mon(exp_vecs, PR)
        reverse!.(exp_vecs)
        res
    else
        gen_mon(exp_vecs, PR)
    end
end
