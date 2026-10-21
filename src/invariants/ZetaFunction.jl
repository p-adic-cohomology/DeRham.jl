"""
    zeta_function(coeffs, q, m)

Given the raw, [`tate_twist`](@ref) or [`normalized_tate_twist`](@ref)
descending coefficients `coeffs` of the middle (``H^m``) L-polynomial ``P(T)``
of a smooth ``m``-dimensional hypersurface ``X`` over ``\\mathbb{F}_q``, returns the zeta
function ``Z(X, T) = P(T)^{(-1)^{m+1}} / \\prod_{i=0}^{m} (1 - q^i T)`` as an
element of `fraction_field(ZZ[T])`.

- `coeffs`: descending coefficients of the middle L-polynomial.
- `q`: size of the base finite field.
- `m`: dimension of the hypersurface.

### Examples

If you already have the coefficients, pass them with the field size and
dimension:

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> coeffs = DeRham.zeta_coefficients(f);

julia> DeRham.zeta_function(coeffs, 7, 1)
(7*T^2 + T + 1)//(7*T^2 - 8*T + 1)
```
"""
function zeta_function(coeffs::AbstractVector, q, m::Integer)
    raw = _normalize_raw_coefficients(coeffs, q)
    P, T = polynomial_ring(ZZ, "T")
    L = L_polynomial(raw; ring = P)
    FF = fraction_field(P)
    denom = prod(1 - ZZ(q)^i * T for i = 0:m)
    if iseven(m)
        return FF(1) // (FF(L) * FF(denom))
    else
        return FF(L) // FF(denom)
    end
end

"""
    zeta_function(f; kwargs...)

[`zeta_function`](@ref) for the smooth hypersurface defined by the
homogeneous polynomial ``f`` (with ``q`` the characteristic of `parent(f)` and
`m = nvars(parent(f)) - 2` its dimension), computed by calling
`zeta_coefficients(f; kwargs...)` and delegating to the coefficient method.
Returns `false` when ``f`` is not smooth, matching `zeta_coefficients`.
The result is an element of `fraction_field(ZZ[T])`.

- `f`: homogeneous polynomial defining the hypersurface.
- `kwargs...`: keywords forwarded to [`zeta_coefficients`](@ref).

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> DeRham.zeta_function(f)
(7*T^2 + T + 1)//(7*T^2 - 8*T + 1)
```
"""
function zeta_function(f::MPolyRingElem; kwargs...)
    q = Int64(characteristic(parent(f)))
    m = nvars(parent(f)) - 2
    coeffs = zeta_coefficients(f; kwargs...)
    coeffs == false && return false
    return zeta_function(coeffs, q, m)
end
