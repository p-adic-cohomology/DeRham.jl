"""
    _hodge_polygon(basis::Array, n)

Calculates the hodge polygon of the cohomology module with
griffiths-dwork basis basis

basis -- an array of "polynomials with pole" as descirbed in PolynomialWithPole.jl
"""
function _hodge_polygon(basis::Vector, n)
    #WRONG: n = highestpoleorder(basis)
    hodgenumbers = zeros(Int, n)
    for i = 0:(n-1)
        h = length(termsoforder(basis, n-i))
        hodgenumbers[i+1] = h
    end

    SlopesPolygon(hodgenumbers)
end

"""
    hodge_polygon(n, d)

Hodge polygon of the primitive middle cohomology of a smooth hypersurface
of degree `d` in `P^n`: the `SlopesPolygon` of Hodge numbers
`h^{n-1-p, p}`, given in closed form by the coefficient of
`t^((p+1)d - n - 1)` in the Jacobian ring's Hilbert series
`(1 + t + ... + t^(d-2))^(n+1)`, for `p = 0, ..., n - 1`.

- `n`: dimension of the ambient projective space.
- `d`: degree of the hypersurface.

Returns a [`SlopesPolygon`](@ref). The slope `p` has horizontal length
`h^{n-1-p,p}`.

For example, a smooth plane cubic has Hodge numbers `(1, 1)`:

```jldoctest
julia> DeRham.hodge_polygon(2, 3).slopelengths == [1, 1]
true
```
"""
function hodge_polygon(n::Integer, d::Integer)
    P, t = polynomial_ring(ZZ, "t")
    hilbert_series = sum(t^i for i = 0:(d-2))^(n + 1)
    hilbert_coeff(e) = e < 0 ? 0 : Int(coeff(hilbert_series, e))
    hodgenumbers = [hilbert_coeff((p + 1) * d - n - 1) for p = 0:(n-1)]
    return SlopesPolygon(hodgenumbers)
end

"""
    hodge_polygon(f; basis=nothing, params=default_params())

Hodge polygon of the primitive middle cohomology of the smooth hypersurface
defined by the homogeneous polynomial `f`.

- `f`: homogeneous polynomial defining the hypersurface.
- `basis`: optional precomputed cohomology basis. If omitted, a basis is
  constructed from `f`.
- `params`: series-precision parameters used when constructing the basis;
  ignored when `basis` is supplied.

Returns a [`SlopesPolygon`](@ref), with one slope for each Hodge number.

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> DeRham.hodge_polygon(f).slopelengths == [1, 1]
true
```
"""
function hodge_polygon(f::MPolyRingElem; basis = nothing, params = default_params())
    n = nvars(parent(f)) - 1
    PR = parent(f)
    R = coefficient_ring(parent(f))

    if basis == nothing
        d = total_degree(f)
        S = collect(0:n) #collect(0:d-1)
        cache = controlled_reduction_cache(n, d, S, params)
        basis = compute_monomial_bases(f, params, cache) # basis of cohomology
    end

    Basis = []
    for i = 1:n
        for j in basis[i]
            push!(Basis, [j, i])
        end
    end

    _hodge_polygon(Basis, n)
end
