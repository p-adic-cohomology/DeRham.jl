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

"""
    ibm_derham_bound(p,m,n)

Returns the exact integer baseline
`f₀(m) = sum(floor(log_p(max(1, m-i))) for i in 1:n)` from
Theorem 3.4.6 of Abbott–Kedlaya–Roe. The integer logarithms are computed
without floating-point arithmetic. The AKR estimate assumes that `p` is prime.
The return value is checked to fit in `Int`.
"""
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
    for i = 1:n_int
        value = max(BigInt(1), m_big - i)
        while value < threshold
            threshold = div(threshold, p_big)
            log_value -= 1
        end
        total += log_value
    end
    _ibm_checked_int(total, "ibm_derham_bound")
end

function _ibm_floor_log_prefix_sum(p::Integer, upper::Integer)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    upper <= 0 && return BigInt(0)

    base = BigInt(p)
    value = BigInt(upper)
    power = base
    exponent = BigInt(1)
    total = BigInt(0)
    while power <= value
        next_power = power * base
        last = min(value, next_power - 1)
        total += exponent * (last - power + 1)
        power = next_power
        exponent += 1
    end
    total
end

function _ibm_f0_from_prefix_sum(p::Integer, m::Integer, n::Integer)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    m >= 1 || throw(ArgumentError("m must be positive"))
    n >= 1 || throw(ArgumentError("n must be positive"))

    m_big = BigInt(m)
    n_big = BigInt(n)
    upper = m_big - 1
    lower = max(BigInt(0), m_big - n_big - 1)
    _ibm_floor_log_prefix_sum(p, upper) - _ibm_floor_log_prefix_sum(p, lower)
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
    p::Integer,
    m::Integer,
    ibm_bound::AbstractVector{<:Integer},
    n::Integer,
    stored_capped_by_f0::Bool,
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
        _ibm_f0_from_prefix_sum(p_big, index, n_big)
    end

    function nn(l::BigInt)
        index = (m_big + l) * p_big
        g = _ibm_binomial_valuation(m_big + l - 1, l, p_big)
        bound_at(index) - l - g
    end

    best = nn(BigInt(1))
    L = BigInt(2)
    next_certificate_check = BigInt(2)
    while true
        exponent = L + best
        if m_big + L >= 2 * n_big + 1 &&
           exponent >= 0 &&
           L >= next_certificate_check &&
           (p_big * (m_big + L) - 1)^n_int <= p_big^exponent
            if !stored_capped_by_f0
                last_stored_l = div(BigInt(stored_length), p_big) - m_big
                for stored_l = L:last_stored_l
                    best = max(best, nn(stored_l))
                end
            end
            return _ibm_checked_int(best, "ibm_constant")
        end
        if m_big + L >= 2 * n_big + 1 && exponent >= 0 && L >= next_certificate_check
            next_certificate_check = 2 * L
        end
        best = max(best, nn(L))
        L += 1
    end
end

function _ibm_crank_bound_h(p::Integer, n::Integer, K::Integer; fixed_point::Bool = false)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    n >= 1 || throw(ArgumentError("n must be positive"))
    K >= 1 || throw(ArgumentError("K must be positive"))

    p_big = BigInt(p)
    n_int = _ibm_checked_int(n, "n")
    K_int = _ibm_checked_int(K, "K")
    bounds = [
        _ibm_checked_int(_ibm_f0_from_prefix_sum(p_big, k, n_int), "ibm_derham_bound")
        for k = 1:K_int
    ]
    m_max = _ibm_checked_int(div(BigInt(K_int), p_big), "maximum m")

    function update_bound!(m, changed_indices)
        pm = _ibm_checked_int(p_big * m, "p*m index")
        N = BigInt(_ibm_constant(p_big, m, bounds, n_int, true))
        base = BigInt(n_int - 1) + bounds[m]
        updated = _ibm_checked_int(min(BigInt(bounds[pm]), max(base, N)), "updated bound")
        if updated < bounds[pm]
            bounds[pm] = updated
            if changed_indices !== nothing
                push!(changed_indices, pm)
            end
        end

        first_backfill = _ibm_checked_int(BigInt(pm) - p_big + 1, "back-fill index")
        for k = first_backfill:pm
            backfilled = min(bounds[k], updated)
            if backfilled < bounds[k]
                bounds[k] = backfilled
                if changed_indices !== nothing
                    push!(changed_indices, k)
                end
            end
        end
    end

    if fixed_point
        pending = collect(1:m_max)
        queued = trues(m_max)
        unqueued = BitSet()
        head = 1
        while head <= length(pending)
            m = pending[head]
            head += 1
            queued[m] = false
            push!(unqueued, m)
            changed_indices = Int[]
            update_bound!(m, changed_indices)
            for index in changed_indices
                if index <= m_max && !queued[index]
                    delete!(unqueued, index)
                    push!(pending, index)
                    queued[index] = true
                end
                if rem(BigInt(index), p_big) == 0
                    last_dependent = min(
                        m_max,
                        _ibm_checked_int(div(BigInt(index), p_big), "dependent index"),
                    )
                    while !isempty(unqueued)
                        dependent = first(unqueued)
                        dependent <= last_dependent || break
                        delete!(unqueued, dependent)
                        push!(pending, dependent)
                        queued[dependent] = true
                    end
                end
            end
        end
    else
        for m = 1:m_max
            update_bound!(m, nothing)
        end
    end

    bounds
end

"""
    ibm_crank_bound(p, n, K)

Returns a certified upper bound on the coordinate lattice V for pole orders
`1:K`, where `K` is the exact stored range. By default, the bound uses one
ordered pass; set `fixed_point=true` to repeat passes on this same finite range
until no entry changes. The H-lattice values are pointwise at most `f₀`; the
return values add `ibm_lattice_correction(p,n)`. The AKR bound assumes that `p`
is prime. During each pass, requests outside `1:K` use exact `f₀` values,
rather than padded storage.
"""
function ibm_crank_bound(p::Integer, n::Integer, K::Integer; fixed_point = false)
    h_bound = _ibm_crank_bound_h(p, n, K; fixed_point = fixed_point)
    correction = ibm_lattice_correction(p, n)
    [
        _ibm_checked_int(BigInt(value) + correction, "coordinate-lattice bound") for
        value in h_bound
    ]
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

"""
    _validate_series_precision_request(p, n, r_m)

Validates the integer requests used by the certified selector in “Certified
series precision: proof of the fixed algorithm.”
"""
function _validate_series_precision_request(
    p::Integer,
    n::Integer,
    r_m::AbstractVector{<:Integer},
)
    p >= 2 || throw(ArgumentError("p must be at least 2"))
    n >= 1 || throw(ArgumentError("n must be positive"))
    length(r_m) == n || throw(ArgumentError("r_m must have length n"))
    all(r -> r >= 0, r_m) || throw(ArgumentError("relative precisions must be nonnegative"))
    if p < n && any(iszero, r_m)
        throw(ArgumentError("zero relative precision requests require p ≥ n"))
    end
    nothing
end

"""
    _ibm_certified_tail_cutoff(p, n, m, r, correction)

Finds the least all-tail cutoff from the integer certificate in “Certified
series precision: proof of the fixed algorithm.”
"""
function _ibm_certified_tail_cutoff(
    p::Integer,
    n::Integer,
    m::Integer,
    r::Integer,
    correction::Integer,
)
    n_int = _ibm_checked_int(n, "n")
    p_big = BigInt(p)
    n_big = BigInt(n)
    m_big = BigInt(m)
    r_big = BigInt(r)
    correction_big = BigInt(correction)

    lower = max(BigInt(0), 2 * n_big + 1 - m_big, r_big + correction_big - m_big + 1)
    upper = max(16 * n_big^2, 2 * (r_big + correction_big + 1), m_big) - m_big

    certified =
        cutoff -> begin
            exponent = m_big - 1 + cutoff - r_big - correction_big
            pole_order = m_big + cutoff
            exponent >= 0 && (p_big * pole_order - 1)^n_int <= p_big^exponent
        end

    lower <= upper || error("certified tail cutoff upper bound is inconsistent")
    while lower < upper
        midpoint = div(lower + upper, 2)
        if certified(midpoint)
            upper = midpoint
        else
            lower = midpoint + 1
        end
    end
    lower
end

"""
    calculate_series_precision(p, n, r_m)

Returns retained-term counts whose omitted Frobenius tails meet every requested
relative precision. The all-tail certificate and its proof are given in
“Certified series precision: proof of the fixed algorithm.”
"""
function calculate_series_precision(p::Integer, n::Integer, r_m::AbstractVector{<:Integer})
    _validate_series_precision_request(p, n, r_m)
    n_int = _ibm_checked_int(n, "n")
    correction = ibm_lattice_correction(p, n)
    cutoffs = [_ibm_certified_tail_cutoff(p, n, m, r_m[m], correction) for m = 1:n_int]

    maximum_pole_order = BigInt(1)
    p_big = BigInt(p)
    for m = 1:n_int
        cutoff = cutoffs[m]
        if cutoff > 0
            maximum_pole_order = max(maximum_pole_order, p_big * (BigInt(m) + cutoff - 1))
        end
    end
    bound_length = _ibm_checked_int(maximum_pole_order, "series bound range")
    bounds = ibm_crank_bound(p, n, bound_length; fixed_point = true)

    selected = Vector{Int}(undef, n_int)
    for m = 1:n_int
        last_failure = BigInt(-1)
        cutoff = cutoffs[m]
        for index = BigInt(0):(cutoff-1)
            pole_index = _ibm_checked_int(p_big * (BigInt(m) + index), "bound index")
            margin = BigInt(m - 1) + index - bounds[pole_index]
            if margin < r_m[m]
                last_failure = index
            end
        end
        certified_count = max(last_failure + 1, BigInt(n_int - m + 1))
        selected[m] = _ibm_checked_int(certified_count, "series precision")
    end
    selected
end

const _series_precision_table_cache_lock = ReentrantLock()
const _series_precision_table_cache = Dict{Any,Vector{Int}}()
const _series_precision_table_cache_limit = 64

"""
    _cached_table_series_precision(p, n, d, r_m, table)

Caches the certified selector result for repeated table-backed requests, while
retaining the table as a lower bound. The selector certificate is given in
“Certified series precision: proof of the fixed algorithm.”
"""
function _cached_table_series_precision(p, n, d, r_m, table)
    key = (BigInt(p), BigInt(n), BigInt(d), Tuple(BigInt(request) for request in r_m))
    cached = lock(_series_precision_table_cache_lock) do
        get(_series_precision_table_cache, key, nothing)
    end
    cached === nothing || return copy(cached)

    calculated = calculate_series_precision(p, n, r_m)
    result = max.(table, calculated)
    lock(_series_precision_table_cache_lock) do
        cached = get(_series_precision_table_cache, key, nothing)
        if cached === nothing
            if length(_series_precision_table_cache) >= _series_precision_table_cache_limit
                empty!(_series_precision_table_cache)
            end
            _series_precision_table_cache[key] = result
            cached = result
        end
        copy(cached)
    end
end

"""
    series_precision(p, n, d, r_m; use_default_table=true)

Uses the certified selector from “Certified series precision: proof of the fixed
algorithm,” while retaining any existing table value as a lower bound.
"""
function series_precision(p, n, d, r_m; use_default_table = true)
    use_default_table || return calculate_series_precision(p, n, r_m)
    table = series_precision_default(p, n, d)
    table === nothing && return calculate_series_precision(p, n, r_m)
    _cached_table_series_precision(p, n, d, r_m, table)
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
