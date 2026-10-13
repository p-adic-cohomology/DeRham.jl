"""
    newton_polygon(coeffs, p)

Given the raw, [`tate_twist`](@ref) or [`normalized_tate_twist`](@ref)
descending coefficients `coeffs` of an L-polynomial, returns its Newton
polygon: the `SlopesPolygon` of `p`-adic valuations of the raw coefficients.
"""
function newton_polygon(coeffs::AbstractVector, p)
    raw = _normalize_raw_coefficients(coeffs, p)
    padic_val = x -> x == 0 ? typemax(Int64) : valuation(x, p)
    SlopesPolygon(raw, padic_val)
end

"""
    newton_polygon(f; kwargs...)

[`newton_polygon`](@ref) for the smooth hypersurface defined by the
homogeneous polynomial `f` (with `p` the characteristic of `parent(f)`),
computed by calling `zeta_coefficients(f; kwargs...)` and delegating to the
coefficient method. Returns `false` when `f` is not smooth, matching
`zeta_coefficients`.
"""
function newton_polygon(f::MPolyRingElem; kwargs...)
    p = Int64(characteristic(parent(f)))
    coeffs = zeta_coefficients(f; kwargs...)
    coeffs == false && return false
    return newton_polygon(coeffs, p)
end
