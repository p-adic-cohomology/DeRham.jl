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
    picard_rank_bound(coeffs::AbstractVector, q)
    picard_rank_bound(f::MPolyRingElem; kwargs...)

Geometric upper bound on the Picard rank of a smooth surface in `P^3`.
The data form takes coefficients of its primitive-`H^2` L-polynomial over
`F_q` and returns `1` plus the degree of the cyclotomic part of `L(T/q)`;
the added `1` accounts for the hyperplane class. The polynomial form computes
those coefficients with `zeta_coefficients`. This is an upper bound
in general and need not be exact. The coefficients may be raw,
`tate_twist` or `normalized_tate_twist` zeta coefficients.

# Arguments
- `coeffs::AbstractVector`: primitive-`H^2` L-polynomial coefficients.
- `f::MPolyRingElem`: homogeneous polynomial defining a smooth surface in
  `P^3`.
- `q`: size of the finite base field for the data form.

# Keywords
- `kwargs...`: forwarded to `zeta_coefficients` by the polynomial
  form.

# Returns
An `Int` giving the geometric upper bound. The polynomial form returns
`false::Bool` when `f` is not smooth, matching `zeta_coefficients`.

# Examples

The polynomial form calculates the bound for the Fermat quartic over `F_7`:

```jldoctest
julia> f = DeRham.fermat_hypersurface(4, 4, 7);

julia> DeRham.picard_rank_bound(f)
22
```

If you already have `zeta_coefficients` output, pass it with the base-field
size:

```jldoctest
julia> f = DeRham.fermat_hypersurface(4, 4, 7);

julia> coeffs = DeRham.zeta_coefficients(f);

julia> DeRham.picard_rank_bound(coeffs, 7)
22
```
"""
function picard_rank_bound(coeffs::AbstractVector, q)
    orders = _cyclotomic_orders(coeffs, q)
    return 1 + sum(euler_phi(m) * mult for (m, mult) in orders; init = 0)
end

"""
    k3_geometric_picard_rank(coeffs::AbstractVector, q)
    k3_geometric_picard_rank(f::MPolyRingElem; kwargs...)

Geometric Picard rank of the quartic K3 surface in `P^3` defined by `f`, or
given by its primitive-`H^2` L-polynomial coefficients over `F_q`. Unlike
[`picard_rank_bound`](@ref), this is exact because the Tate conjecture is
known for K3 surfaces over finite fields. The result includes the
hyperplane class. The data form requires 22 coefficients, which may be raw,
`tate_twist` or `normalized_tate_twist` zeta coefficients. The polynomial
form computes the coefficients with `zeta_coefficients`.

# Arguments
- `coeffs::AbstractVector`: 22 primitive-`H^2` L-polynomial coefficients.
- `f::MPolyRingElem`: homogeneous quartic in four variables defining a K3
  surface in `P^3`.
- `q`: size of the finite base field for the data form.

# Keywords
- `kwargs...`: forwarded to `zeta_coefficients` by the polynomial
  form.

# Returns
An `Int` giving the exact geometric Picard rank. The polynomial form returns
`false::Bool` when `f` is not smooth; it throws `ArgumentError` unless `f`
is a quartic in four variables.

# Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(4, 4, 7);

julia> DeRham.k3_geometric_picard_rank(f)
22
```

If you already have `zeta_coefficients` output, reuse it with `q = 7`:

```jldoctest
julia> f = DeRham.fermat_hypersurface(4, 4, 7);

julia> coeffs = DeRham.zeta_coefficients(f);

julia> DeRham.k3_geometric_picard_rank(coeffs, 7)
22
```
"""
function k3_geometric_picard_rank(coeffs::AbstractVector, q)
    _check_k3_length(coeffs)
    return picard_rank_bound(coeffs, q)
end

"""
    k3_picard_rank(coeffs::AbstractVector, q, k)
    k3_picard_rank(f::MPolyRingElem, k; kwargs...)

Picard rank of the quartic K3 surface in `P^3` defined by `f`, or given by
its primitive-`H^2` L-polynomial coefficients over `F_q`, after extension to
`F_{q^k}`. The result is exact because the Tate conjecture is known for K3
surfaces over finite fields and includes the hyperplane class. The data form
requires 22 coefficients, which may be raw, `tate_twist` or
`normalized_tate_twist` zeta coefficients. The polynomial form computes the
coefficients with `zeta_coefficients`.

# Arguments
- `coeffs::AbstractVector`: 22 primitive-`H^2` L-polynomial coefficients.
- `f::MPolyRingElem`: homogeneous quartic in four variables defining a K3
  surface in `P^3`.
- `q`: size of the finite base field for the data form.
- `k`: positive extension degree.

# Keywords
- `kwargs...`: forwarded to `zeta_coefficients` by the polynomial
  form.

# Returns
An `Int` giving the Picard rank over `F_{q^k}`. The polynomial form returns
`false::Bool` when `f` is not smooth; it throws `ArgumentError` unless `f`
is a quartic in four variables.

# Examples

The Fermat quartic over `F_7` reaches its geometric rank over `F_{7^2}`:

```jldoctest
julia> f = DeRham.fermat_hypersurface(4, 4, 7);

julia> DeRham.k3_picard_rank(f, 2)
22
```

If you already have `zeta_coefficients` output, reuse it with the extension
degree:

```jldoctest
julia> f = DeRham.fermat_hypersurface(4, 4, 7);

julia> coeffs = DeRham.zeta_coefficients(f);

julia> DeRham.k3_picard_rank(coeffs, 7, 2)
22
```
"""
function k3_picard_rank(coeffs::AbstractVector, q, k)
    _check_k3_length(coeffs)
    orders = _cyclotomic_orders(coeffs, q)
    return 1 + sum(euler_phi(m) * mult for (m, mult) in orders if k % m == 0; init = 0)
end

"""
    k3_picard_realization_degree(coeffs::AbstractVector, q)
    k3_picard_realization_degree(f::MPolyRingElem; kwargs...)

Least extension degree `k` for which the Picard rank of the quartic K3
surface in `P^3` defined by `f` equals its geometric Picard rank, or the same
degree computed from its primitive-`H^2` L-polynomial coefficients over
`F_q`. The data form requires 22 coefficients, which may be raw,
`tate_twist` or `normalized_tate_twist` zeta coefficients. The polynomial
form computes the coefficients with `zeta_coefficients`.

# Arguments
- `coeffs::AbstractVector`: 22 primitive-`H^2` L-polynomial coefficients.
- `f::MPolyRingElem`: homogeneous quartic in four variables defining a K3
  surface in `P^3`.
- `q`: size of the finite base field for the data form.

# Keywords
- `kwargs...`: forwarded to `zeta_coefficients` by the polynomial
  form.

# Returns
An `Int` giving the least extension degree. The polynomial form returns
`false::Bool` when `f` is not smooth; it throws `ArgumentError` unless `f`
is a quartic in four variables.

# Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(4, 4, 7);

julia> DeRham.k3_picard_realization_degree(f)
2
```

If you already have `zeta_coefficients` output, reuse it with `q = 7`:

```jldoctest
julia> f = DeRham.fermat_hypersurface(4, 4, 7);

julia> coeffs = DeRham.zeta_coefficients(f);

julia> DeRham.k3_picard_realization_degree(coeffs, 7)
2
```
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
