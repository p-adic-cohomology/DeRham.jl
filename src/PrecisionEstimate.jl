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

"""
    ibm_derham_bound(p,m,n)

Returns the exact integer baseline
`f₀(m) = sum(floor(log_p(max(1, m-i))) for i in 1:n)` from
Theorem 3.4.6 of Abbott–Kedlaya–Roe. The integer logarithms are computed
without floating-point arithmetic. The AKR estimate assumes that `p` is prime.
The return value is checked to fit in `Int`.
"""
function _ibm_checked_int(value::Integer, description::AbstractString)
    typemin(Int) <= value <= typemax(Int) ||
        throw(OverflowError("$description does not fit in Int"))
    Int(value)
end

function _ibm_floor_log(p::Integer, a::Integer)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    a >= 1 || throw(ArgumentError("a must be positive"))

    base = BigInt(p)
    value = BigInt(a)
    exponent = 0
    while value >= base
        value = div(value, base)
        exponent += 1
    end
    exponent
end

function ibm_derham_bound(p::Integer, m::Integer, n::Integer)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    m >= 1 || throw(ArgumentError("m must be positive"))
    n >= 1 || throw(ArgumentError("n must be positive"))

    n_int = _ibm_checked_int(n, "n")
    p_big = BigInt(p)
    m_big = BigInt(m)
    total = BigInt(0)
    log_value = _ibm_floor_log(p_big, max(BigInt(1), m_big - 1))
    threshold = p_big^log_value
    for i in 1:n_int
        value = max(BigInt(1), m_big - i)
        while value < threshold
            threshold = div(threshold, p_big)
            log_value -= 1
        end
        total += log_value
    end
    _ibm_checked_int(total, "ibm_derham_bound")
end

function _ibm_vp_factorial(n::Integer, p::Integer)
    n >= 0 || throw(ArgumentError("factorial valuation requires n ≥ 0"))
    p >= 2 || throw(ArgumentError("p must be at least 2"))

    quotient = BigInt(n)
    prime = BigInt(p)
    result = BigInt(0)
    while quotient > 0
        quotient = div(quotient, prime)
        result += quotient
    end
    result
end

function _ibm_binomial_valuation(n::Integer, k::Integer, p::Integer)
    0 <= k <= n || throw(ArgumentError("invalid binomial arguments"))
    _ibm_vp_factorial(n, p) - _ibm_vp_factorial(k, p) - _ibm_vp_factorial(n - k, p)
end

"""
    ibm_lattice_correction(p, n)

Returns the conservative correction from the AKR H-lattice bound to the
coordinate lattice V bound, `v_p((n-1)!)` when `p < n`, and zero otherwise.
The AKR correction assumes that `p` is prime.
"""
function ibm_lattice_correction(p::Integer, n::Integer)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    n >= 1 || throw(ArgumentError("n must be positive"))
    p >= n && return 0
    _ibm_checked_int(_ibm_vp_factorial(n - 1, p), "ibm_lattice_correction")
end

"""
    ibm_constant(p, m, ibm_bound, n)

Computes the exact supremum of `A((m+l)p) - l - v_p(binomial(m+l-1,l))`
over all `l ≥ 1`, where `A` is the stored upper bound and is `f₀` outside the
stored range. Here `p` is prime. The tail ends only after the all-l certificate
`(p(m+L)-1)^n ≤ p^(L+S)` holds, with `L ≥ 2`, `m+L ≥ 2n+1`, and `L+S ≥ 0`;
the certificate bounds every term for `l ≥ L` by the maximum `S` already seen.
"""
function ibm_constant(
    p::Integer,
    m::Integer,
    ibm_bound::AbstractVector{<:Integer},
    n::Integer,
)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    m >= 1 || throw(ArgumentError("m must be positive"))
    n >= 1 || throw(ArgumentError("n must be positive"))

    p_big = BigInt(p)
    m_big = BigInt(m)
    n_big = BigInt(n)
    n_int = _ibm_checked_int(n_big, "n")
    stored_length = length(ibm_bound)

    function bound_at(index::BigInt)
        if index <= stored_length
            return BigInt(ibm_bound[_ibm_checked_int(index, "bound index")])
        end
        BigInt(ibm_derham_bound(p_big, index, n_int))
    end

    function nn(l::BigInt)
        index = (m_big + l) * p_big
        g = _ibm_binomial_valuation(m_big + l - 1, l, p_big)
        bound_at(index) - l - g
    end

    best = nn(BigInt(1))
    L = BigInt(2)
    while true
        exponent = L + best
        if m_big + L >= 2 * n_big + 1 &&
           exponent >= 0 &&
           (p_big * (m_big + L) - 1)^n_int <= p_big^exponent
            return _ibm_checked_int(best, "ibm_constant")
        end
        best = max(best, nn(L))
        L += 1
    end
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
function _ibm_crank_bound_h(p::Integer, n::Integer, K::Integer)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    n >= 1 || throw(ArgumentError("n must be positive"))
    K >= 1 || throw(ArgumentError("K must be positive"))

    p_big = BigInt(p)
    n_int = _ibm_checked_int(n, "n")
    K_int = _ibm_checked_int(K, "K")
    bounds = [ibm_derham_bound(p_big, k, n_int) for k in 1:K_int]
    m_max = _ibm_checked_int(div(BigInt(K_int), p_big), "maximum m")

    for m in 1:m_max
        pm = _ibm_checked_int(p_big * m, "p*m index")
        N = BigInt(ibm_constant(p_big, m, bounds, n_int))
        base = BigInt(n_int - 1) + bounds[m]
        updated = _ibm_checked_int(min(BigInt(bounds[pm]), max(base, N)), "updated bound")
        bounds[pm] = updated

        first_backfill = _ibm_checked_int(BigInt(pm) - p_big + 1, "back-fill index")
        for k in first_backfill:pm
            bounds[k] = min(bounds[k], updated)
        end
    end

    bounds
end

function ibm_crank_bound(p::Integer, n::Integer, K::Integer)
    h_bound = _ibm_crank_bound_h(p, n, K)
    correction = ibm_lattice_correction(p, n)
    [_ibm_checked_int(BigInt(value) + correction, "coordinate-lattice bound") for value in h_bound]
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
