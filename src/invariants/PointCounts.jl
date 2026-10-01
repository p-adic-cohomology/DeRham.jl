"""
Counts points in the most naive way possible

Returns a tuple containing
a list of the points of the affine
variety k[x_1,...,x_n]/f,
and the number of projective points
"""
function naivelypointcount(f, q)
    if !is_prime(q)
        error("naive point counting not implemented for non-prime fields yet")
    end

    rationalpoints = []

    n = length(gens(parent(f)))
    iter = Iterators.product(fill(0:(q-1), n)...)
    for xs in iter
        val = evaluate(f, collect(xs))
        if val == 0
            push!(rationalpoints, xs)
        end
    end

    n_aff_points = length(rationalpoints)
    (rationalpoints, divexact(n_aff_points - 1, q-1))
end


"""
    vp(a, p)
    Returns the p-adic valuation of a

    INPUTS:
    * "a" -- integer
    * "p" -- integer, a prime number
"""
function vp(a, p)
    @assert is_prime(p)
    if a == 0
        return Inf
    end

    e = 0
    while mod(a, p) == 0
        a = div(a, p)
        e = e + 1
    end

    return e
end

"""
    _power_sums_from_L_polynomial(coeffs, k)

Given the raw (`ZZRingElem`) descending coefficients `coeffs` of an
L-polynomial `P(T) = prod_j (1 - alpha_j T)`, returns `[p_1, ..., p_k]`
where `p_i = sum_j alpha_j^i` is the `i`-th power sum of the inverse roots
`alpha_j`, computed via Newton's identities in exact integer arithmetic
(no division).
"""
function _power_sums_from_L_polynomial(coeffs::AbstractVector{ZZRingElem}, k::Integer)
    d = length(coeffs) - 1
    e(i) = i == 0 ? ZZ(1) : (i <= d ? ((-1)^i) * coeffs[d+1-i] : ZZ(0))
    p = Vector{ZZRingElem}(undef, k)
    for kk = 1:k
        s = ZZ(0)
        for i = 1:(kk-1)
            s += (-1)^(i - 1) * e(i) * p[kk-i]
        end
        s += (-1)^(kk - 1) * kk * e(kk)
        p[kk] = s
    end
    return p
end

"""
    point_counts(coeffs, q, m, k)

Given the raw, [`tate_twist`](@ref) or [`normalized_tate_twist`](@ref)
descending coefficients `coeffs` of the middle (`H^m`) L-polynomial of a
smooth `m`-dimensional hypersurface over `F_q`, returns
`[#X(F_q), #X(F_{q^2}), ..., #X(F_{q^k})]` via
`#X(F_{q^i}) = sum_{t=0}^{m} q^{t i} + (-1)^m p_i`, where `p_i` is the
`i`-th power sum of the middle-cohomology Frobenius eigenvalues (Newton
identities on `coeffs`, in exact integer arithmetic).
"""
function point_counts(coeffs::AbstractVector, q, m::Integer, k::Integer)
    raw = _normalize_raw_coefficients(coeffs, q)
    p = _power_sums_from_L_polynomial(raw, k)
    qq = ZZ(q)
    return [sum(qq^(t * i) for t = 0:m) + (-1)^m * p[i] for i = 1:k]
end

"""
    point_counts(f, k; kwargs...)

[`point_counts`](@ref) for the smooth hypersurface defined by the
homogeneous polynomial `f` (with `q` the characteristic of `parent(f)` and
`m = nvars(parent(f)) - 2` its dimension), computed by calling
`zeta_coefficients(f; kwargs...)` and delegating to the coefficient method.
Returns `false` when `f` is not smooth, matching `zeta_coefficients`.
"""
function point_counts(f::MPolyRingElem, k::Integer; kwargs...)
    q = Int64(characteristic(parent(f)))
    m = nvars(parent(f)) - 2
    coeffs = zeta_coefficients(f; kwargs...)
    coeffs == false && return false
    return point_counts(coeffs, q, m, k)
end
