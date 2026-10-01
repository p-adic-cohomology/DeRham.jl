"""
    _twist_coefficient(c, q, e)

Returns `c * q^e` as a `QQFieldElem`, computed exactly (dividing rather than
raising `q` to a negative power when `e < 0`).
"""
function _twist_coefficient(c, q, e::Integer)
    if e >= 0
        return QQ(c) * QQ(q)^e
    else
        return QQ(c) // QQ(q)^(-e)
    end
end

"""
    _to_exact_ZZ(val, errmsg)

Converts a `QQFieldElem` that is exactly an integer to `ZZRingElem`, or
throws `ArgumentError(errmsg)` if it is not integral.
"""
function _to_exact_ZZ(val::QQFieldElem, errmsg)
    isone(denominator(val)) || throw(ArgumentError(errmsg))
    return numerator(val)
end

"""
    tate_twist(coeffs, q; j = 1)

Given the descending `zeta_coefficients`-style coefficients `coeffs` of an
L-polynomial `P(T)`, returns the coefficients (same orientation, as
`QQFieldElem`, ending in `1`) of the Tate twist `P(T / q^j)`.
"""
function tate_twist(coeffs::AbstractVector, q; j::Integer = 1)
    n = length(coeffs)
    deg = n - 1
    return [_twist_coefficient(coeffs[i], q, -j * (deg + 1 - i)) for i = 1:n]
end

"""
    tate_untwist(coeffs, q; j = 1, s = 0)

Inverts both a Tate twist and a normalization: divides `coeffs` by `q^s`,
then returns the coefficients (as `ZZRingElem`) of `P(q^j T)`, where `P` is
the L-polynomial with coefficients `coeffs / q^s`. Throws `ArgumentError` if
the result is not integral.
"""
function tate_untwist(coeffs::AbstractVector, q; j::Integer = 1, s::Integer = 0)
    n = length(coeffs)
    deg = n - 1
    return [
        _to_exact_ZZ(
            _twist_coefficient(coeffs[i], q, j * (deg + 1 - i) - s),
            "tate_untwist: result is not integral at coefficient index $i",
        ) for i = 1:n
    ]
end

"""
    normalized_tate_twist(coeffs, q, s; j = 1)

Given the descending `zeta_coefficients`-style coefficients `coeffs` of an
L-polynomial `P(T)`, returns the coefficients (as `ZZRingElem`) of
`q^s * P(T / q^j)`. Throws `ArgumentError` if that is not integral. With
`s = 1, j = 1` this is the old "boat shape" L-polynomial. Note: a
normalized twist with `s = 0` is indistinguishable from raw coefficients and
is not supported by [`_normalize_raw_coefficients`](@ref)'s detection.
"""
function normalized_tate_twist(coeffs::AbstractVector, q, s::Integer; j::Integer = 1)
    n = length(coeffs)
    deg = n - 1
    return [
        _to_exact_ZZ(
            _twist_coefficient(coeffs[i], q, s - j * (deg + 1 - i)),
            "normalized_tate_twist: q^$s * P(T/q^$j) is not integral at coefficient index $i",
        ) for i = 1:n
    ]
end

"""
    _normalize_raw_coefficients(coeffs, q)

Detects which of the three supported input kinds `coeffs` is — raw
`zeta_coefficients` output, [`tate_twist`](@ref) output (with `j = 1`), or
[`normalized_tate_twist`](@ref) output (with `j = 1`, some `s >= 1`) — and
returns the underlying raw `ZZRingElem` coefficients. Every data form of a
coefficient-taking invariant calls this first, so it accepts all three
kinds of input. Detection: `QQFieldElem` element type means `tate_twist`
output; integer coefficients ending in `1` are already raw; integer
coefficients ending in `q^s` for some `s >= 1` are `normalized_tate_twist`
output with that `s`. A normalized twist with `s = 0` is indistinguishable
from raw and is not supported.
"""
function _normalize_raw_coefficients(coeffs::AbstractVector, q)
    if coeffs isa AbstractVector{<:QQFieldElem}
        return tate_untwist(coeffs, q; j = 1)
    end

    last = ZZ(coeffs[end])
    isone(last) && return ZZRingElem[ZZ(c) for c in coeffs]

    qq = ZZ(q)
    val = last
    s = 0
    while mod(val, qq) == 0
        val = divexact(val, qq)
        s += 1
    end
    isone(val) || throw(
        ArgumentError(
            "coefficients are not raw, tate_twist, or normalized_tate_twist(..., s) output for q=$q (trailing coefficient $last)",
        ),
    )
    return tate_untwist(coeffs, q; j = 1, s = s)
end
