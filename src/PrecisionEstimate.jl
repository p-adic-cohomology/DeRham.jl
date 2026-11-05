# module PrecisionEstimate

#p = 41
#r = 7 # p-adic precision
#n = 2


"""
    frobenius_precision()

Determines the necessary p-adic precision r for the Frobenius matrix
"""
function frobenius_precision(k, q)
    if k == 2 && q == 7
        r = 2
    elseif k == 21 && q == 7
        r = 10
    end
    return r
end

function _ibm_floor_log(p::Int, a::Int)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    a >= 1 || throw(ArgumentError("a must be positive"))

    value = a
    exponent = 0
    while value >= p
        value = div(value, p)
        exponent += 1
    end
    exponent
end

"""
    ibm_derham_bound(p,m,n)

Returns the exact integer baseline
`f₀(m) = sum(floor(log_p(max(1, m-i))) for i in 1:n)` from
Theorem 3.4.6 of Abbott–Kedlaya–Roe. The integer logarithms are computed
without floating-point arithmetic. The AKR estimate assumes that `p` is prime.
The result is computed in `Int`.
"""
function ibm_derham_bound(p::Int, m::Int, n::Int)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    m >= 1 || throw(ArgumentError("m must be positive"))
    n >= 1 || throw(ArgumentError("n must be positive"))

    total = 0
    log_value = _ibm_floor_log(p, max(1, m - 1))
    threshold = 1
    for _ = 1:log_value
        threshold *= p
    end
    for i = 1:n
        value = max(1, m - i)
        while value < threshold
            threshold = div(threshold, p)
            log_value -= 1
        end
        total += log_value
    end
    total
end

function _ibm_vp_factorial(n::Int, p::Int)
    n >= 0 || throw(ArgumentError("factorial valuation requires n ≥ 0"))
    p >= 2 || throw(ArgumentError("p must be at least 2"))

    quotient = n
    result = 0
    while quotient > 0
        quotient = div(quotient, p)
        result += quotient
    end
    result
end

function _ibm_binomial_valuation(n::Int, k::Int, p::Int)
    0 <= k <= n || throw(ArgumentError("invalid binomial arguments"))
    _ibm_vp_factorial(n, p) - _ibm_vp_factorial(k, p) - _ibm_vp_factorial(n - k, p)
end

"""
    ibm_lattice_correction(p, n)

Returns the conservative correction from the AKR H-lattice bound to the
coordinate lattice V bound, `v_p((n-1)!)` when `p < n`, and zero otherwise.
The AKR correction assumes that `p` is prime.
"""
function ibm_lattice_correction(p::Int, n::Int)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    n >= 1 || throw(ArgumentError("n must be positive"))
    p >= n && return 0
    _ibm_vp_factorial(n - 1, p)
end

"""
    ibm_constant(p, m, ibm_bound, n)

Computes the exact supremum of `A((m+l)p) - l - v_p(binomial(m+l-1,l))`
over all `l ≥ 1`, where `A` is the stored upper bound and is `f₀` outside the
stored range. Here `p` is prime. The tail ends only after the all-l certificate
`(p(m+L)-1)^n ≤ p^(L+S)` holds, with `L ≥ 2`, `m+L ≥ 2n+1`, and `L+S ≥ 0`;
the certificate bounds every term for `l ≥ L` by the maximum `S` already seen.
"""
function ibm_constant(p::Int, m::Int, ibm_bound::AbstractVector{Int}, n::Int)
    _ibm_constant(p, m, ibm_bound, n, false)
end

"""
    _ibm_constant(p, m, ibm_bound, n, stored_capped_by_f0)

Shared implementation of `ibm_constant`. When `stored_capped_by_f0` is `false`,
the stored terms from the certificate threshold `L` through the end of the
stored range are also included, so the result is the exact supremum for an
arbitrary stored vector. When it is `true`, the caller guarantees
`ibm_bound[k] ≤ f₀(k)` for every stored `k`. Each stored term with `l ≥ L` is
then at most the corresponding `f₀` term, which the certificate already bounds
by the returned maximum, so that scan cannot change the result and is omitted.
`_ibm_crank_bound_h` satisfies the guarantee because its bounds start at `f₀`
and are only ever replaced by minima.
The prefix maximum is updated at every `L`, while the expensive tail certificate
is checked at exponentially spaced checkpoints; intervening candidate terms are
still included in the maximum.
"""
function _ibm_constant(
    p::Int,
    m::Int,
    ibm_bound::AbstractVector{Int},
    n::Int,
    stored_capped_by_f0::Bool,
)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    m >= 1 || throw(ArgumentError("m must be positive"))
    n >= 1 || throw(ArgumentError("n must be positive"))

    stored_length = length(ibm_bound)

    function bound_at(index::Int)
        if index <= stored_length
            return ibm_bound[index]
        end
        ibm_derham_bound(p, index, n)
    end

    function nn(l::Int)
        index = (m + l) * p
        g = _ibm_binomial_valuation(m + l - 1, l, p)
        bound_at(index) - l - g
    end

    best = nn(1)
    L = 2
    next_certificate_check = 2
    while true
        if BigInt(m) + L >= 2 * BigInt(n) + 1 &&
           BigInt(L) + best >= 0 &&
           L >= next_certificate_check &&
           (BigInt(p) * (BigInt(m) + L) - 1)^n <= BigInt(p)^(BigInt(L) + best)
            if !stored_capped_by_f0
                last_stored_l = div(stored_length, p) - m
                for stored_l = L:last_stored_l
                    best = max(best, nn(stored_l))
                end
            end
            return best
        end
        if BigInt(m) + L >= 2 * BigInt(n) + 1 &&
           BigInt(L) + best >= 0 &&
           L >= next_certificate_check
            next_certificate_check = 2 * L
        end
        best = max(best, nn(L))
        L += 1
    end
end

function _ibm_crank_bound_h(p::Int, n::Int, K::Int)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    n >= 1 || throw(ArgumentError("n must be positive"))
    K >= 1 || throw(ArgumentError("K must be positive"))

    bounds = [ibm_derham_bound(p, k, n) for k = 1:K]
    m_max = div(K, p)

    for m = 1:m_max
        pm = p * m
        N = _ibm_constant(p, m, bounds, n, true)
        base = n - 1 + bounds[m]
        updated = min(bounds[pm], max(base, N))
        bounds[pm] = updated

        first_backfill = pm - p + 1
        for k = first_backfill:pm
            bounds[k] = min(bounds[k], updated)
        end
    end

    bounds
end

"""
    ibm_crank_bound(p, n, K)

Returns a single-pass certified upper bound on the coordinate lattice V for
pole orders `1:K`, where `K` is the exact stored range. The one-pass H-lattice
bound is pointwise at most `f₀`; the return values add
`ibm_lattice_correction(p,n)`. No equality with the fully converged AKR bound
or minimality is asserted. The AKR bound assumes that `p` is prime. During the
pass, requests outside `1:K` use exact `f₀` values, rather than padded storage.
"""
function ibm_crank_bound(p::Int, n::Int, K::Int)
    h_bound = _ibm_crank_bound_h(p, n, K)
    correction = ibm_lattice_correction(p, n)
    [value + correction for value in h_bound]
end

"""
    calculate_relative_precision(HP, slope, hodge_numbers, weight, p)

Calculates the vector of relative precisions r_m

INPUTS:
* "polygon" -- A SlopesPolygon struct describing a polygon whose values describe the
divisibility of the roots. Note: it is traditional to use the hodge polygon here,
but if you know the Newton Polygon (e.g. you have a K3 surface and you already have the
Artin-Mazur height) then the Newton Polygon will do just as well.
* "weight" -- integer, the motivic weight of the cohomology group, equals to the dimension of the hypersurface
* "p" -- integer, prime number
"""
function calculate_relative_precision(polygon, weight, p)

    #   * "HP" -- list, the i-th item corresponds to the height of above i in the Hodge polygon
    #   * "slope" -- list, the i-th item corresponds to the slope of the i-th segment in the Hodge polygon
    #   * "hodge_numbers" -- list, the list of Hodge numbers
    HP = polygon.values
    slope = Int.(polygon.slopesbefore)
    hodge_numbers = polygon.slopelengths

    r_vector = [0 for i = 1:length(hodge_numbers)]
    Pdeg = length(HP) - 1  # degree of the L-polynomial P_n(X,T)
    max_digits = 0

    for i = 1:(div(Pdeg, 2)+1)
        r = ceil(Int, log(2*Pdeg/i)/log(p) + i*weight*0.5) - HP[i+1]

        if r >= max_digits
            max_digits = r
            for j = 0:slope[i+1]
                r_vector[j+1] = r
            end

            for j = (slope[i+1]+1):(length(hodge_numbers)-1)
                r = r-1
                r_vector[j+1] = r
            end
        end
    end

    return reverse(r_vector)
end

function calculate_series_precision(p, n, r_m)

    # TODO: prove that this is correct
    #
    # We may assume in (1.8) of Costa's thesis that the worst case is i=0

    bounds = ibm_crank_bound(p, n, p*n)

    #TODO: why does ibm_derham_bound sometimes give answers that are better than the crank?

    N = [(r_m[m] == 0 ? 0 : r_m[m] - m + 1 + bounds[p*m]) for m = 1:n]

    #println("Series precision: $N")
    N
end

function series_precision(p, n, d, r_m; use_default_table = true)
    if use_default_table
        N = series_precision_default(p, n, d)
        if N !== nothing
            return N
        end
    end
    calculate_series_precision(p, n, r_m)
end

"""
    algorithm_precision()

Determines M such that the algorithm works in Z/p^M Z
"""
function algorithm_precision(p, n, d, r_m, N_m)
    s_m = [i+x-1 for (i, x) in enumerate(N_m)]
    s_m_valuation = [valuation(ZZ(factorial(big(p*s-1))), ZZ(p)) for s in s_m]

    maximum([r_m[m] + s_m_valuation[m] - m + 1 for m = 1:length(r_m)])
end
