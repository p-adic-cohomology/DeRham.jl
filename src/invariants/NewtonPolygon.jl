"""
    newton_polygon(coeffs, p)

Given the raw, [`tate_twist`](@ref), or [`normalized_tate_twist`](@ref)
descending coefficients `coeffs` of an L-polynomial, returns its Newton
polygon. The coefficients are normalized to the raw convention before the
lower convex hull of the points `(i, v_p(c_i))` is formed; the slopes are the
`p`-adic valuations of the reciprocal roots.

- `coeffs`: descending L-polynomial coefficients.
- `p`: prime characteristic used for the `p`-adic valuation.

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
homogeneous polynomial `f`. It obtains the characteristic `p` from
`parent(f)`, calls `zeta_coefficients(f; kwargs...)`, and delegates to the
coefficient method.

- `f`: homogeneous polynomial defining the hypersurface.
- `kwargs...`: keywords forwarded to [`zeta_coefficients`](@ref).

Returns a [`SlopesPolygon`](@ref), or `false` when `f` is not smooth, matching
`zeta_coefficients`.

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> DeRham.newton_polygon(f).slopes == [0//1, 1//1]
true
```
"""
function newton_polygon(f::MPolyRingElem; kwargs...)
    p = Int64(characteristic(parent(f)))
    coeffs = zeta_coefficients(f; kwargs...)
    coeffs == false && return false
    return newton_polygon(coeffs, p)
end
