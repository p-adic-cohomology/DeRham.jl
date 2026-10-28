"""
    a_number(FM, h; column_moduli=nothing)

For a smooth curve, returns the dimension of the kernel of its Hasse-Witt
matrix ``H``: ``h - \\operatorname{rank}(H)``.

- `FM`: Frobenius matrix computed at relative precision 1.
- `h`: genus of the curve.

Returns the curve's a-number as an integer.

When `FM` comes from `frobenius_matrix_with_precision`, pass its per-column
known moduli if they differ; only the selected columns need precision divisible
by ``p``. Request them with `return_column_moduli=true` on
`frobenius_matrix_with_precision`.

If you already have a Frobenius matrix, pass it with the curve's genus:

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 5);

julia> FM, column_moduli = DeRham.frobenius_matrix_with_precision(
           f, [1, 0]; return_column_moduli = true);

julia> DeRham.a_number(FM, 1; column_moduli = column_moduli)
1
```
"""
function a_number(FM::Union{MatElem,AbstractMatrix}, h::Integer; column_moduli = nothing)
    return h - rank(hasse_witt_matrix(FM, h; column_moduli = column_moduli))
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

    result = frobenius_matrix_with_precision(
        f,
        r_m;
        basis = basis,
        verbose = verbose,
        return_column_moduli = true,
        kwargs...,
    )
    result === false && return false
    FM, column_moduli = result

    a_number(FM, g; column_moduli = column_moduli)
end
