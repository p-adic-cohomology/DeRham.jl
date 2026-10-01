"""
    _require_plane_curve(f, fname)

Throws `ArgumentError` unless `f` defines a plane curve (three variables,
i.e. `nvars(parent(f)) == 3`). Every `f` form of a Jacobian invariant calls
this first, since `H^1(J) = H^1(C)` is only meaningful for curves `C`.
"""
function _require_plane_curve(f::MPolyRingElem, fname::AbstractString)
    n = nvars(parent(f)) - 1
    n == 2 || throw(
        ArgumentError(
            "$fname: f must define a plane curve (3 variables), got $(n + 1) variables",
        ),
    )
end

"""
    _elementary_symmetric_from_power_sums(p, d)

Given the power sums `p[1], ..., p[d]` of a multiset of `d` algebraic
integers `beta_1, ..., beta_d`, returns `[e_0, ..., e_d]`, their elementary
symmetric functions, via the Newton-Girard recurrence `e_t = (1/t) *
sum_{i=1}^{t} (-1)^(i-1) e_{t-i} p_i` (with `e_0 = 1`), in exact arithmetic
(`divexact` throughout — exact because the `beta_j` here are always
Frobenius eigenvalues or products thereof, whose elementary symmetric
functions are rational integers).
"""
function _elementary_symmetric_from_power_sums(p::AbstractVector{ZZRingElem}, d::Integer)
    e = Vector{ZZRingElem}(undef, d + 1)
    e[1] = ZZ(1)
    for t = 1:d
        s = ZZ(0)
        for i = 1:t
            s += (-1)^(i - 1) * e[t-i+1] * p[i]
        end
        e[t+1] = divexact(s, ZZ(t))
    end
    return e
end

"""
    _raw_coeffs_from_elementary(e)

Given `e = [e_0, ..., e_d]`, the elementary symmetric functions of a
multiset of `d` algebraic integers `beta_j`, returns the descending
`zeta_coefficients`-style coefficients of `prod_j (1 - beta_j T)` (raw,
ending in `e_0 = 1`): `coeffs[d+1-t] = (-1)^t * e_t`.
"""
function _raw_coeffs_from_elementary(e::AbstractVector{ZZRingElem})
    d = length(e) - 1
    coeffs = Vector{ZZRingElem}(undef, d + 1)
    for t = 0:d
        coeffs[d+1-t] = (-1)^t * e[t+1]
    end
    return coeffs
end

"""
    _jacobian_wedge_charpolys(coeffs, q; ring = nothing)

Given the raw, [`tate_twist`](@ref) or [`normalized_tate_twist`](@ref)
descending coefficients `coeffs` of the `H^1` L-polynomial of a smooth
curve over `F_q` (with Frobenius eigenvalues `alpha_1, ..., alpha_{2g}`),
returns `[P_0, ..., P_{2g}]`, where `P_i` is the characteristic polynomial
`prod_{S} (1 - (prod_{j in S} alpha_j) T)` of Frobenius on `wedge^i H^1`
(the product over `i`-subsets `S` of `{1, ..., 2g}`), as an element of
`ring` (a fresh `ZZ[T]` if `ring === nothing`). Each `P_i` is built from
`alpha_j^k`-power sums via two applications of
[`_elementary_symmetric_from_power_sums`](@ref): first, for each `k`, the
elementary symmetric functions `e_i(alpha_1^k, ..., alpha_{2g}^k)` (whose
`i`-th entry is `p_k(wedge^i)`, the `k`-th power sum of the `wedge^i`
eigenvalues, since `(prod_{j in S} alpha_j)^k = prod_{j in S} alpha_j^k`);
second, those `p_k(wedge^i)` (`k = 1, ..., binomial(2g, i)`) are themselves
turned into the elementary symmetric functions of the `wedge^i`
eigenvalues, i.e. the coefficients of `P_i`. No subset is ever enumerated.
"""
function _jacobian_wedge_charpolys(coeffs::AbstractVector, q; ring = nothing)
    raw = _normalize_raw_coefficients(coeffs, q)
    d = length(raw) - 1
    max_d_i = maximum(binomial(d, i) for i = 0:d)
    p = _power_sums_from_L_polynomial(raw, d * max_d_i)

    e_by_k =
        [_elementary_symmetric_from_power_sums([p[k*t] for t = 1:d], d) for k = 1:max_d_i]

    P = ring === nothing ? polynomial_ring(ZZ, "T")[1] : ring
    return [
        begin
            d_i = binomial(d, i)
            wedge_power_sums = [e_by_k[k][i+1] for k = 1:d_i]
            e_wedge = _elementary_symmetric_from_power_sums(wedge_power_sums, d_i)
            L_polynomial(_raw_coeffs_from_elementary(e_wedge); ring = P)
        end for i = 0:d
    ]
end

"""
    jacobian_zeta_coefficients(f; kwargs...)

The `H^1` `zeta_coefficients` of the Jacobian `J` of the smooth plane
curve `C` defined by the homogeneous polynomial `f`, i.e.
`zeta_coefficients(f; kwargs...)`, since `H^1(J) = H^1(C)`. There is no
data form, since it would be the identity. Throws `ArgumentError` when `f`
does not define a plane curve.
"""
function jacobian_zeta_coefficients(f::MPolyRingElem; kwargs...)
    _require_plane_curve(f, "jacobian_zeta_coefficients")
    return zeta_coefficients(f; kwargs...)
end

"""
    jacobian_point_counts(coeffs, q, k)

Given the raw, [`tate_twist`](@ref) or [`normalized_tate_twist`](@ref)
descending coefficients `coeffs` of the `H^1` L-polynomial of a smooth
curve over `F_q` (`H^1(J) = H^1(C)` for its Jacobian `J`), returns
`[#J(F_q), #J(F_{q^2}), ..., #J(F_{q^k})]` via `#J(F_{q^i}) = L_i(1) =
prod_j (1 - alpha_j^i)`, with `L_i` built exactly from the `alpha_j^i`
power sums via Newton's identities (see
[`_elementary_symmetric_from_power_sums`](@ref)), in exact integer
arithmetic.
"""
function jacobian_point_counts(coeffs::AbstractVector, q, k::Integer)
    raw = _normalize_raw_coefficients(coeffs, q)
    d = length(raw) - 1
    p = _power_sums_from_L_polynomial(raw, d * k)
    return [
        begin
            e = _elementary_symmetric_from_power_sums([p[i*t] for t = 1:d], d)
            sum((-1)^t * e[t+1] for t = 0:d)
        end for i = 1:k
    ]
end

"""
    jacobian_point_counts(f, k; kwargs...)

[`jacobian_point_counts`](@ref) for the Jacobian `J` of the smooth plane
curve defined by the homogeneous polynomial `f` (with `q` the
characteristic of `parent(f)`), computed by calling `zeta_coefficients(f;
kwargs...)` and delegating to the coefficient method. Returns `false` when
`f` is not smooth, matching `zeta_coefficients`. Throws `ArgumentError`
when `f` does not define a plane curve.
"""
function jacobian_point_counts(f::MPolyRingElem, k::Integer; kwargs...)
    _require_plane_curve(f, "jacobian_point_counts")
    q = Int64(characteristic(parent(f)))
    coeffs = zeta_coefficients(f; kwargs...)
    coeffs == false && return false
    return jacobian_point_counts(coeffs, q, k)
end

"""
    jacobian_zeta_function(coeffs, q)

Given the raw, [`tate_twist`](@ref) or [`normalized_tate_twist`](@ref)
descending coefficients `coeffs` of the `H^1` L-polynomial of a smooth
genus-`g` curve over `F_q` (`H^1(J) = H^1(C)` for its Jacobian `J`),
returns `Z(J, T) = prod_{i=0}^{2g} P_i(T)^((-1)^(i+1))`, where `P_i` is the
characteristic polynomial of Frobenius on `wedge^i H^1` (see
[`_jacobian_wedge_charpolys`](@ref)), as an element of
`fraction_field(ZZ[T])`. Warns when `g > 10`, since the cost grows like
`binomial(2g, g)`.
"""
function jacobian_zeta_function(coeffs::AbstractVector, q)
    raw = _normalize_raw_coefficients(coeffs, q)
    d = length(raw) - 1
    g = d ÷ 2
    g > 10 && @warn "jacobian_zeta_function: g = $g (cost grows like binomial(2g, g))"

    P, T = polynomial_ring(ZZ, "T")
    FF = fraction_field(P)
    P_list = _jacobian_wedge_charpolys(coeffs, q; ring = P)
    num = FF(1)
    den = FF(1)
    for (i, P_i) in enumerate(P_list)
        if isodd(i - 1)
            num *= FF(P_i)
        else
            den *= FF(P_i)
        end
    end
    return num // den
end

"""
    jacobian_zeta_function(f; kwargs...)

[`jacobian_zeta_function`](@ref) for the Jacobian `J` of the smooth plane
curve defined by the homogeneous polynomial `f` (with `q` the
characteristic of `parent(f)`), computed by calling `zeta_coefficients(f;
kwargs...)` and delegating to the coefficient method. Returns `false` when
`f` is not smooth, matching `zeta_coefficients`. Throws `ArgumentError`
when `f` does not define a plane curve.
"""
function jacobian_zeta_function(f::MPolyRingElem; kwargs...)
    _require_plane_curve(f, "jacobian_zeta_function")
    q = Int64(characteristic(parent(f)))
    coeffs = zeta_coefficients(f; kwargs...)
    coeffs == false && return false
    return jacobian_zeta_function(coeffs, q)
end
