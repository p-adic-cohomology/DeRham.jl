"""
    artin_mazur_height(coeffs, q)

Computes the Artin-Mazur height of a quartic K3 surface given the
coefficients `coeffs` of its primitive-`H^2` L-polynomial over `F_q`.
`coeffs` may be raw, `tate_twist` or `normalized_tate_twist` zeta
coefficients ([`_normalize_raw_coefficients`](@ref) decides) and must have
length 22.
"""
function artin_mazur_height(coeffs, q)
    @assert length(coeffs) == 22
    raw = _normalize_raw_coefficients(coeffs, q)
    boatshape_coeffs = normalized_tate_twist(raw, q, 1)
    for i = 2:11
        if mod(boatshape_coeffs[i], q) != 0
            return i - 1
        end
    end
    return Inf
end

"""
    artin_mazur_height(f; kwargs...)

[`artin_mazur_height`](@ref) for the quartic K3 surface in `P^3` defined by
the homogeneous polynomial `f`, computed by calling `zeta_coefficients(f;
kwargs...)` (with `q` the characteristic of `parent(f)`) and delegating to
the coefficient method. Returns `false` when `f` is not smooth, matching
`zeta_coefficients`.
"""
function artin_mazur_height(f::MPolyRingElem; kwargs...)
    q = Int64(characteristic(parent(f)))
    coeffs = zeta_coefficients(f; kwargs...)
    coeffs == false && return false
    return artin_mazur_height(coeffs, q)
end
