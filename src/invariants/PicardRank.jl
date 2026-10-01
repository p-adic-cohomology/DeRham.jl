"""
    _is_boat_shape(coeffs, q)

Decides whether `coeffs` is already boat-shape L-polynomial data (`q*P(T/q)`,
trailing coefficient `q`) rather than raw `zeta_coefficients` output (trailing
coefficient `1`, since the raw L-polynomial is always monic-at-the-constant-term).
"""
function _is_boat_shape(coeffs, q)
    return coeffs[end] == q
end

"""
    _cyclotomic_orders(coeffs, q)

Factors the boat-shape L-polynomial `q*P(T/q)` over `QQ` and returns a
`Dict{Int,Int}` mapping each cyclotomic order `m` found among its roots to the
multiplicity of `Phi_m` in the factorization; non-cyclotomic irreducible
factors are ignored. `coeffs` may be raw or boat-shape zeta coefficients
(`_is_boat_shape` decides, normalizing raw input to boat shape first). Built
on Oscar/Hecke's `is_cyclotomic_polynomial_with_data` applied to the monic
form of each irreducible factor.
"""
function _cyclotomic_orders(coeffs, q)
    n = length(coeffs)
    bshape = _is_boat_shape(coeffs, q) ? coeffs : boat_shape_Lpoly(coeffs, n - 1, q)
    Qx, Tq = polynomial_ring(QQ, "T")
    poly = sum(Tq^(n - i) * Qx(bshape[i]) for i = 1:n)
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
