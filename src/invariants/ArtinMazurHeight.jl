"""
    artin_mazur_height(coeffs::AbstractVector, q)
    artin_mazur_height(f::MPolyRingElem; kwargs...)

Artin–Mazur height of the quartic K3 surface in ``\\mathbb{P}^3`` defined by ``f``, or
given by the coefficients of its primitive-``H^2`` L-polynomial over ``\\mathbb{F}_q``.
The data form requires 22 coefficients, which may be raw,
`tate_twist` or `normalized_tate_twist` zeta coefficients. The polynomial
form computes them with `zeta_coefficients`.

# Arguments
- `coeffs::AbstractVector`: 22 primitive-``H^2`` L-polynomial coefficients.
- `f::MPolyRingElem`: homogeneous quartic in four variables defining a K3
  surface in ``\\mathbb{P}^3``.
- `q`: size of the finite base field for the data form.

# Keywords
- `kwargs...`: forwarded to `zeta_coefficients` by the polynomial
  form.

# Returns
An `Int` from 1 to 10, or `Inf::Float64` for infinite height. The polynomial
form also returns `false::Bool` when ``f`` is not smooth, matching
`zeta_coefficients`.

# Examples

The Fermat quartic over ``\\mathbb{F}_5`` is ordinary and has height 1:

```jldoctest
julia> f = DeRham.fermat_hypersurface(4, 4, 5);

julia> DeRham.artin_mazur_height(f)
1
```

If you already have `zeta_coefficients` output, reuse it with ``q = 5``:

```jldoctest
julia> f = DeRham.fermat_hypersurface(4, 4, 5);

julia> coeffs = DeRham.zeta_coefficients(f);

julia> DeRham.artin_mazur_height(coeffs, 5)
1
```

The Fermat quartic over ``\\mathbb{F}_7`` is supersingular and has infinite height:

```jldoctest
julia> f = DeRham.fermat_hypersurface(4, 4, 7);

julia> DeRham.artin_mazur_height(f)
Inf
```
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

[`artin_mazur_height`](@ref) for the quartic K3 surface in ``\\mathbb{P}^3`` defined by
the homogeneous polynomial ``f``, computed by calling `zeta_coefficients(f;
kwargs...)` (with ``q`` the characteristic of `parent(f)`) and delegating to
the coefficient method. Returns `false` when ``f`` is not smooth, matching
`zeta_coefficients`.
"""
function artin_mazur_height(f::MPolyRingElem; kwargs...)
    q = Int64(characteristic(parent(f)))
    coeffs = zeta_coefficients(f; kwargs...)
    coeffs == false && return false
    return artin_mazur_height(coeffs, q)
end
