"""
    _cyclotomic_orders(coeffs, q)

Factors the boat-shape L-polynomial `q*P(T/q)` over `QQ` and returns a
`Dict{Int,Int}` mapping each cyclotomic order `m` found among its roots to the
multiplicity of `Phi_m` in the factorization; non-cyclotomic irreducible
factors are ignored. `coeffs` may be raw, `tate_twist` or
`normalized_tate_twist` zeta coefficients ([`_normalize_raw_coefficients`](@ref)
decides). Built on Oscar/Hecke's `is_cyclotomic_polynomial_with_data` applied
to the monic form of each irreducible factor.
"""
function _cyclotomic_orders(coeffs, q)
    raw = _normalize_raw_coefficients(coeffs, q)
    bshape = normalized_tate_twist(raw, q, 1)
    Qx, Tq = polynomial_ring(QQ, "T")
    poly = L_polynomial(bshape; ring = Qx)
    orders = Dict{Int,Int}()
    for (fac, mult) in factor(poly)
        monic = fac * inv(leading_coefficient(fac))
        is_cyc, m = is_cyclotomic_polynomial_with_data(monic)
        if is_cyc
            orders[m] = get(orders, m, 0) + mult
        end
    end
    return orders
end

function _check_k3_length(coeffs)
    length(coeffs) == 22 || throw(
        ArgumentError(
            "k3_* Picard functions require 22 coefficients (the quartic K3 primitive H^2 L-polynomial), got $(length(coeffs))",
        ),
    )
end

"""
    picard_rank_bound(coeffs, q)

Geometric upper bound on the Picard rank of the smooth surface in `P^3` whose
primitive-`H^2` L-polynomial (from `zeta_coefficients`, with the `+1` below
accounting for the hyperplane class) has coefficients `coeffs` over `F_q`:
`1` plus the degree of the cyclotomic part of `L(T/q)`. This is an upper
bound for the Picard rank in general, with equality under the Tate
conjecture. `coeffs` may be raw or boat-shape zeta coefficients.
"""
function picard_rank_bound(coeffs::AbstractVector, q)
    orders = _cyclotomic_orders(coeffs, q)
    return 1 + sum(euler_phi(m) * mult for (m, mult) in orders; init = 0)
end

"""
    k3_geometric_picard_rank(coeffs, q)

Geometric Picard rank of a quartic K3 surface whose primitive-`H^2`
L-polynomial (from `zeta_coefficients`, with the `+1` below accounting for
the hyperplane class) has coefficients `coeffs` over `F_q`. Unlike
[`picard_rank_bound`](@ref), this is exact because the Tate conjecture is
known for K3 surfaces over finite fields. `coeffs` may be raw or boat-shape
zeta coefficients and must have length 22.
"""
function k3_geometric_picard_rank(coeffs::AbstractVector, q)
    _check_k3_length(coeffs)
    return picard_rank_bound(coeffs, q)
end

"""
    k3_picard_rank(coeffs, q, k)

Picard rank of a quartic K3 surface over `F_{q^k}`, given the primitive-`H^2`
L-polynomial (from `zeta_coefficients`, with the `+1` below accounting for
the hyperplane class) coefficients `coeffs` over `F_q`: `1` plus the degree
of the cyclotomic factors of `L(T/q)` whose order divides `k`. Exact because
the Tate conjecture is known for K3 surfaces over finite fields. `coeffs` may
be raw or boat-shape zeta coefficients and must have length 22.
"""
function k3_picard_rank(coeffs::AbstractVector, q, k)
    _check_k3_length(coeffs)
    orders = _cyclotomic_orders(coeffs, q)
    return 1 + sum(euler_phi(m) * mult for (m, mult) in orders if k % m == 0; init = 0)
end

"""
    k3_picard_realization_degree(coeffs, q)

Least `k` such that the Picard rank of the quartic K3 surface over
`F_{q^k}` (see [`k3_picard_rank`](@ref)) equals its geometric Picard rank
(see [`k3_geometric_picard_rank`](@ref)): the lcm of the cyclotomic orders
present in `L(T/q)`, or `1` if none are. `coeffs` may be raw or boat-shape
zeta coefficients and must have length 22.
"""
function k3_picard_realization_degree(coeffs::AbstractVector, q)
    _check_k3_length(coeffs)
    orders = _cyclotomic_orders(coeffs, q)
    return isempty(orders) ? 1 : reduce(lcm, keys(orders))
end

function _check_k3_quartic_fourvar(f)
    (total_degree(f) == 4 && nvars(parent(f)) == 4) || throw(
        ArgumentError(
            "k3_* Picard functions require a quartic surface in 4 variables (a hypersurface in P^3), got degree $(total_degree(f)) in $(nvars(parent(f))) variables",
        ),
    )
end

"""
    picard_rank_bound(f; kwargs...)

[`picard_rank_bound`](@ref) for the smooth surface in `P^3` defined by the
homogeneous polynomial `f`, computed by calling `zeta_coefficients(f;
kwargs...)` (with `q` the characteristic of `parent(f)`) and delegating to
the coefficient method. Returns `false` when `f` is not smooth, matching
`zeta_coefficients`.
"""
function picard_rank_bound(f::MPolyRingElem; kwargs...)
    q = Int64(characteristic(parent(f)))
    zf = zeta_coefficients(f; kwargs...)
    zf == false && return false
    return picard_rank_bound(zf, q)
end

"""
    k3_geometric_picard_rank(f; kwargs...)

[`k3_geometric_picard_rank`](@ref) for the quartic K3 surface in `P^3`
defined by the homogeneous polynomial `f`, computed by calling
`zeta_coefficients(f; kwargs...)` (with `q` the characteristic of
`parent(f)`) and delegating to the coefficient method. Returns `false` when
`f` is not smooth, matching `zeta_coefficients`.
"""
function k3_geometric_picard_rank(f::MPolyRingElem; kwargs...)
    _check_k3_quartic_fourvar(f)
    q = Int64(characteristic(parent(f)))
    zf = zeta_coefficients(f; kwargs...)
    zf == false && return false
    return k3_geometric_picard_rank(zf, q)
end

"""
    k3_picard_rank(f, k; kwargs...)

[`k3_picard_rank`](@ref) for the quartic K3 surface in `P^3` defined by the
homogeneous polynomial `f`, computed by calling `zeta_coefficients(f;
kwargs...)` (with `q` the characteristic of `parent(f)`) and delegating to
the coefficient method. Returns `false` when `f` is not smooth, matching
`zeta_coefficients`.
"""
function k3_picard_rank(f::MPolyRingElem, k; kwargs...)
    _check_k3_quartic_fourvar(f)
    q = Int64(characteristic(parent(f)))
    zf = zeta_coefficients(f; kwargs...)
    zf == false && return false
    return k3_picard_rank(zf, q, k)
end

"""
    k3_picard_realization_degree(f; kwargs...)

[`k3_picard_realization_degree`](@ref) for the quartic K3 surface in `P^3`
defined by the homogeneous polynomial `f`, computed by calling
`zeta_coefficients(f; kwargs...)` (with `q` the characteristic of
`parent(f)`) and delegating to the coefficient method. Returns `false` when
`f` is not smooth, matching `zeta_coefficients`.
"""
function k3_picard_realization_degree(f::MPolyRingElem; kwargs...)
    _check_k3_quartic_fourvar(f)
    q = Int64(characteristic(parent(f)))
    zf = zeta_coefficients(f; kwargs...)
    zf == false && return false
    return k3_picard_realization_degree(zf, q)
end
