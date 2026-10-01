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
