"""
    frobenius_matrix_with_precision(f, r::Integer; basis=nothing, verbose=0, kwargs...)

Compute the Frobenius matrix for the smooth hypersurface defined by `f` to
the requested relative `p`-adic precision `r`. The scalar request is applied
to every Hodge component. Unless `basis` is supplied, the matrix is expressed
in the internal cohomology basis computed by `get_basis_of_cohomology(f)`;
that helper and its basis ordering are not public API. Additional keywords
are forwarded to [`zeta_coefficients`](@ref).

- `f`: homogeneous polynomial defining the hypersurface.
- `r`: requested relative `p`-adic precision for each Hodge component.
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
`r_m` gives one requested relative `p`-adic precision for each Hodge
component and must have the length required by the selected basis. The matrix
is expressed in that basis, which defaults to DeRham's internal
`get_basis_of_cohomology(f)` basis. Additional keywords are forwarded to
[`zeta_coefficients`](@ref).

- `f`: homogeneous polynomial defining the hypersurface.
- `r_m`: vector of requested relative precisions, one per Hodge component.
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
    precision_info = precision_information_user_r_m(f, r; basis = basis, verbose = verbose)
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
