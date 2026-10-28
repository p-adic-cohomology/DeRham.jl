"""
    newton_polygon(coeffs, p)

Given the raw, [`tate_twist`](@ref), or [`normalized_tate_twist`](@ref)
descending coefficients `coeffs` of an L-polynomial, returns its Newton
polygon. The coefficients are normalized to the raw convention before the
lower convex hull of the points ``(i, v_p(c_i))`` is formed; the slopes are the
``p``-adic valuations of the reciprocal roots.

- `coeffs`: descending L-polynomial coefficients.
- `p`: prime characteristic used for the ``p``-adic valuation.

Returns a [`SlopesPolygon`](@ref) whose slopes and horizontal lengths encode
the Newton polygon.

If you already have `zeta_coefficients` output for a curve, pass it with its
characteristic:

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> coeffs = DeRham.zeta_coefficients(f);

julia> DeRham.newton_polygon(coeffs, 7).slopes == [0//1, 1//1]
true
```
"""
function newton_polygon(coeffs::AbstractVector, p)
    raw = _normalize_raw_coefficients(coeffs, p)
    padic_val = x -> x == 0 ? typemax(Int64) : valuation(x, p)
    SlopesPolygon(raw, padic_val)
end

"""
    newton_polygon(f; kwargs...)

[`newton_polygon`](@ref) for the smooth hypersurface defined by the
homogeneous polynomial ``f``. For a quartic K3 surface in four variables over
a prime field, it uses [`artin_mazur_height`](@ref) to construct the slopes
directly. Other inputs obtain the characteristic ``p`` from `parent(f)`, call
`zeta_coefficients(f; kwargs...)`, and delegate to the coefficient method.

- `f`: homogeneous polynomial defining the hypersurface.
- `kwargs...`: keywords forwarded to [`zeta_coefficients`](@ref).

Returns a [`SlopesPolygon`](@ref), or `false` when ``f`` is not smooth, matching
`zeta_coefficients`.

For a quartic K3 surface over a prime field, the polynomial form uses its
Artin–Mazur height to construct the slopes directly. For finite height ``h``,
the slopes are ``1 - 1/h``, ``1``, and ``1 + 1/h`` with multiplicities
``h``, ``21 - 2h``, and ``h``. Infinite height gives 21 slopes equal to 1.

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> DeRham.newton_polygon(f).slopes == [0//1, 1//1]
true
```
"""
function newton_polygon(f::MPolyRingElem; kwargs...)
    p = Int64(characteristic(parent(f)))
    if _is_prime_field_k3_quartic(f)
        h = artin_mazur_height(f; kwargs...)
        h === false && return false
        if isinf(h)
            slopes = Rational{Int}[1//1]
            lengths = [21]
        else
            slopes = Rational{Int}[(h-1)//h, 1//1, (h+1)//h]
            lengths = [h, 21 - 2h, h]
        end
        return SlopesPolygon(slopes, lengths, values(slopes, lengths)...)
    end
    coeffs = zeta_coefficients(f; kwargs...)
    coeffs == false && return false
    return newton_polygon(coeffs, p)
end
