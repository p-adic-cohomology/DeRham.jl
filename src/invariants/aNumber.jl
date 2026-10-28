"""
    a_number(FM, h)

For a smooth curve, returns the dimension of the kernel of its Hasse-Witt
matrix ``H``: ``h - \\operatorname{rank}(H)``.

- `FM`: Frobenius matrix computed at relative precision 1.
- `h`: genus of the curve.

Returns the curve's a-number as an integer.

When `FM` comes from `frobenius_matrix_with_precision`, its relative precision
request should be `[1, 0]` or stronger. The matrix-level method assumes this
precision contract; see the [Frobenius matrix documentation](@ref).

If you already have a Frobenius matrix, pass it with the curve's genus:

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 5);

julia> FM = DeRham.frobenius_matrix_with_precision(f, [1, 0]);

julia> DeRham.a_number(FM, 1)
1
```
"""
function a_number(FM::MatElem, h::Integer)
    return h - rank(hasse_witt_matrix(FM, h))
end

"""
    a_number(f; verbose = 0, kwargs...)

Computes the a-number of the smooth plane curve defined by the homogeneous
polynomial ``f``, via its Hasse-Witt matrix. This method is currently supported
only for curves.

- `f`: homogeneous polynomial in three variables defining the curve.
- `verbose`: diagnostic verbosity forwarded to the Frobenius computation.
- `kwargs...`: additional keywords forwarded to
  [`frobenius_matrix_with_precision`](@ref) and its coefficient computation.

Returns the curve's a-number as an integer.

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 5);

julia> DeRham.a_number(f)
1
```
"""
function a_number(f::MPolyRingElem; verbose = 0, kwargs...)
    n = nvars(parent(f)) - 1

    if n != 2
        throw(ArgumentError("a-numbers are only supported for curves so far"))
    end

    if 0 < verbose
        println("Calculating Hasse-Witt matrix...")
    end

    r_m = fill(0, n)
    r_m[1] = 1

    basis = get_basis_of_cohomology(f)

    hodge_polygon = _hodge_polygon(basis, n)
    hodge_numbers = hodge_polygon.slopelengths
    g = hodge_numbers[1]

    FM =
        frobenius_matrix_with_precision(f, r_m; basis = basis, verbose = verbose, kwargs...)
    FM === false && return false

    a_number(FM, g)
end
