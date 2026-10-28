"""
    hasse_witt_matrix(FM, h)

Given a Frobenius matrix `FM` computed at relative precision 1 and the curve's
genus ``h``, returns the Hasse-Witt matrix: the last ``h \\times h`` block of `FM`,
reduced modulo ``p``. Here ``p`` is the unique prime factor of the modulus of
the first entry's parent ring.

- `FM`: Frobenius matrix computed at relative precision 1. Each column has a
  prime-power parent modulus, and all columns use the same prime.
- `h`: genus of the curve, which is the size of the Hasse-Witt matrix.

Returns an ``h \\times h`` matrix over ``\\mathbb{F}_p``.

If you already have a Frobenius matrix, pass it with the curve's genus:

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> FM = DeRham.frobenius_matrix_with_precision(f, [1, 0]);

julia> size(DeRham.hasse_witt_matrix(FM, 1))
(1, 1)
```
"""
function hasse_witt_matrix(FM::Union{MatElem,AbstractMatrix}, h::Integer)
    n = modulus(parent(FM[1, 1]))
    (p, _) = first(factor(ZZ(n)))

    highprec = FM[(end-h+1):end, (end-h+1):end]

    K = GF(Int(p))

    result = matrix_space(K, size(highprec, 1), size(highprec, 2))()
    for column in axes(highprec, 2), row in axes(highprec, 1)
        result[row, column] = K(lift(ZZ, highprec[row, column]))
    end
    result
end

"""
    hasse_witt_matrix(f; verbose = 0, kwargs...)

Hasse-Witt matrix of the smooth plane curve defined by the homogeneous
polynomial ``f``, computed via `frobenius_matrix_with_precision(f, r_m;
kwargs...)` with `r_m = [1, 0]`.

- `f`: homogeneous polynomial in three variables defining the curve.
- `verbose`: diagnostic verbosity forwarded to the Frobenius computation.
- `kwargs...`: additional keywords forwarded to
  [`frobenius_matrix_with_precision`](@ref) and its coefficient computation.

Returns a matrix over ``\\mathbb{F}_p``, where ``p`` is the characteristic of the field.
For a curve of genus ``g``, the matrix is ``g \\times g``.

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> rank(DeRham.hasse_witt_matrix(f))
1
```
"""
function hasse_witt_matrix(f::MPolyRingElem; verbose = 0, kwargs...)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    if d < n
        throw(
            ArgumentError(
                "Hasse-Witt matrices are no supported for varieties of negative Kodaira dimensin (i.e. degree < number of variables",
            ),
        )
    end

    r_m = fill(0, n)
    r_m[1] = 1

    basis = get_basis_of_cohomology(f)

    hodge_polygon = _hodge_polygon(basis, n)
    hodge_numbers = hodge_polygon.slopelengths
    fh = hodge_numbers[1]

    F = frobenius_matrix_with_precision(f, r_m, basis = basis, verbose = verbose, kwargs...)

    hasse_witt_matrix(F, fh)
end
