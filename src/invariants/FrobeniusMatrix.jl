"""
    frobenius_matrix_with_precision(f, r::Integer; basis=nothing, verbose=0, kwargs...)

Compute the Frobenius matrix for the smooth hypersurface defined by ``f`` to
the requested relative ``p``-adic precision ``r``. The scalar request is applied
to every Hodge-filtered piece. The returned matrix meets at least this relative
precision for each piece; an existing automatic table floor can provide more.
The finite-tail certificate is documented in the
[Frobenius matrix precision proof](@ref sec-certified-series-truncation), based
on the full proof report, “Certified series precision: proof of the fixed algorithm.”
Unless `basis` is supplied,
the matrix is expressed in the internal cohomology basis computed by
`get_basis_of_cohomology(f)`;
that helper and its basis ordering are not public API. Additional keywords
are forwarded to [`zeta_coefficients`](@ref).

- `f`: homogeneous polynomial defining the hypersurface.
- `r`: requested relative ``p``-adic precision for each Hodge-filtered piece.
- `basis`: optional basis to use instead of computing the default basis.
- `verbose`: diagnostic verbosity, forwarded to the computation.
- `kwargs...`: additional keywords forwarded to [`zeta_coefficients`](@ref).

Returns the Frobenius matrix, or `false` if the input does not meet the
required smoothness condition.

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> size(DeRham.frobenius_matrix_with_precision(f, 2))
(2, 2)
```
"""
function frobenius_matrix_with_precision(
    f,
    r::Integer;
    basis = nothing,
    verbose = 0,
    kwargs...,
)
    r_m = fill(Int(r), nvars(parent(f)) - 1)
    return frobenius_matrix_with_precision(
        f,
        r_m;
        basis = basis,
        verbose = verbose,
        kwargs...,
    )
end

"""
    frobenius_matrix_with_precision(f, r_m::AbstractVector; basis=nothing, verbose=0, kwargs...)

Vector-precision overload of [`frobenius_matrix_with_precision`](@ref).
`r_m` gives one requested relative ``p``-adic precision for each
Hodge-filtered piece and must have the length required by the selected basis.
The returned matrix meets at least these relative precisions, using the
certificate in the
[Frobenius matrix precision proof](@ref sec-certified-series-truncation); an
existing automatic table floor can provide more. An all-zero request returns the Hodge-bound zero matrix; zero
requests require `p ≥ nvars(parent(f)) - 1`. The matrix is expressed in that
basis, which defaults to DeRham's internal
`get_basis_of_cohomology(f)` basis. Additional keywords are forwarded to
[`zeta_coefficients`](@ref).

Returns an ordinary Oscar matrix. See the Frobenius matrix documentation for
its per-column precision contract.

- `f`: homogeneous polynomial defining the hypersurface.
- `r_m`: vector of requested relative precisions, one per Hodge-filtered piece.
- `basis`: optional basis to use instead of computing the default basis.
- `verbose`: diagnostic verbosity, forwarded to the computation.
- `kwargs...`: additional keywords forwarded to [`zeta_coefficients`](@ref).

"""
function frobenius_matrix_with_precision(
    f,
    r::AbstractVector;
    basis = nothing,
    verbose = 0,
    kwargs...,
)
    p = Int(characteristic(parent(f)))
    n = nvars(parent(f)) - 1
    r_m = [Int(x) for x in r]
    _validate_series_precision_request(p, n, r_m)
    if all(iszero, r_m)
        check_smoothness(f) || return false
        cohomology_basis = basis === nothing ? get_basis_of_cohomology(f) : basis
        hodge_polygon = _hodge_polygon(cohomology_basis, n)
        hodge_numbers = hodge_polygon.slopelengths
        first_nonzero_piece = findfirst(x -> x > 0, hodge_numbers)
        absolute_precision = first_nonzero_piece === nothing ? 1 : first_nonzero_piece
        precision_ring, = residue_ring(ZZ, BigInt(p)^absolute_precision)
        dimension = sum(hodge_numbers)
        return zero_matrix(precision_ring, dimension, dimension)
    end
    precision_info =
        precision_information_user_r_m(f, r_m; basis = basis, verbose = verbose)
    result = zeta_coefficients(
        f;
        precision_info = precision_info,
        givefrobmat = true,
        verbose = verbose,
        kwargs...,
    )
    if result == false
        return false
    end
    return result[1]
end
