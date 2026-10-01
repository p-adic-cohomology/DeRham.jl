"""
    L_polynomial(coeffs; ring = nothing)

Builds the L-polynomial as an element of `ZZ[T]` from its descending
`zeta_coefficients`-style coefficients `coeffs` (`coeffs[i]` is the
coefficient of `T^(deg + 1 - i)`, where `deg = length(coeffs) - 1`; the
constant term `coeffs[end]` is `1` for raw input). Uses `ring` (an
`MPolyRing`/`PolyRing` over `ZZ` or a compatible coefficient ring) if given,
else constructs a fresh `polynomial_ring(ZZ, "T")`.
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
homogeneous polynomial `f`, computed by calling `zeta_coefficients(f;
kwargs...)` and delegating to the coefficient method. Returns `false` when
`f` is not smooth, matching `zeta_coefficients`.
"""
function L_polynomial(f::MPolyRingElem; ring = nothing, kwargs...)
    coeffs = zeta_coefficients(f; kwargs...)
    coeffs == false && return false
    return L_polynomial(coeffs; ring = ring)
end
