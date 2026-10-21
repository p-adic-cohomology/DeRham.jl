"""
    L_polynomial(coeffs; ring = nothing)

Builds the L-polynomial as an element of `ZZ[T]` from its descending
`zeta_coefficients`-style coefficients `coeffs` (`coeffs[i]` is the
coefficient of ``T^{\\deg + 1 - i}``, where ``\\deg`` is `length(coeffs) - 1`; the
constant term `coeffs[end]` is `1` for raw input). Uses `ring` (an
`MPolyRing`/`PolyRing` over `ZZ` or a compatible coefficient ring) if given,
else constructs a fresh `polynomial_ring(ZZ, "T")`.

- `coeffs`: descending L-polynomial coefficients, such as the output of
  [`zeta_coefficients`](@ref).
- `ring`: optional polynomial ring in which to construct the result.

Returns the polynomial in the selected ring.

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> coeffs = DeRham.zeta_coefficients(f);

julia> DeRham.L_polynomial(coeffs)
7*T^2 + T + 1
```
"""
function L_polynomial(coeffs::AbstractVector; ring = nothing)
    if ring === nothing
        P, T = polynomial_ring(ZZ, "T")
    else
        P = ring
        T = gens(P)[1]
    end
    deg = length(coeffs) - 1
    return sum(T^(deg + 1 - i) * ZZ(coeffs[i]) for i = 1:(deg+1))
end

"""
    L_polynomial(f; ring = nothing, kwargs...)

[`L_polynomial`](@ref) for the smooth hypersurface defined by the
homogeneous polynomial ``f``, computed by calling `zeta_coefficients(f;
kwargs...)` and delegating to the coefficient method. Returns `false` when
``f`` is not smooth, matching `zeta_coefficients`. Otherwise, returns the
middle L-polynomial in a polynomial ring over `ZZ` (or the supplied `ring`).

- `f`: homogeneous polynomial defining the hypersurface.
- `ring`: optional polynomial ring for the result.
- `kwargs...`: keywords forwarded to [`zeta_coefficients`](@ref).

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> DeRham.L_polynomial(f)
7*T^2 + T + 1
```
"""
function L_polynomial(f::MPolyRingElem; ring = nothing, kwargs...)
    coeffs = zeta_coefficients(f; kwargs...)
    coeffs == false && return false
    return L_polynomial(coeffs; ring = ring)
end
