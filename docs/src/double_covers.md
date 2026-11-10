# Double covers (experimental)

This page describes experimental and provisional functionality. The frontend
`DeRham.DoubleCover` and its methods are not exported or declared public.

Given a smooth double cover of projective space branched along a homogeneous
polynomial $F$ of even degree over a finite field of odd characteristic,
`DeRham.zeta_coefficients` computes the L-polynomial of the primitive,
anti-invariant middle cohomology. Its weight is the dimension $n$ of the base
projective space. The supported input has even branch degree, odd
characteristic, and a smooth branch polynomial over $\mathbb{F}_p$.

The double-cover computation requires explicit `series_precision` and
`algorithm_precision` values. `relative_precision` is optional and defaults to
the existing Weil-bound formula. These precision choices are not certified:
there is no precision proof for this frontend yet. Only the depth-first CPU
algorithm is supported.

Each precision may be given as a scalar for all pieces or, where supported, as
a vector. `relative_precision` vectors have one entry per Hodge polygon piece
in the order used by `calculate_relative_precision`. `series_precision`
vectors have one entry per pole order in ascending order.
`algorithm_precision` is a single integer.

The only supported entry points for this frontend are
`DeRham.zeta_coefficients` and `DeRham.get_basis_of_cohomology`. Other
invariants are not implemented for double covers.

## A hyperelliptic example

The quartic branch polynomial below defines a smooth genus-one double cover over
$\mathbb{F}_3$. The coefficient vector is in descending order.

```jldoctest double_cover_curve
julia> using DeRham, Oscar

julia> R, (x, z) = Oscar.polynomial_ring(Oscar.GF(3), [:x, :z]);

julia> F = x^4 + z^4;

julia> DeRham.zeta_coefficients(
           DeRham.DoubleCover(F);
           algorithm = :depthfirst,
           series_precision = 2,
           algorithm_precision = 5,
       )
3-element Vector{ZZRingElem}:
 3
 0
 1
```

## Frontend and computation docstrings

```@docs
DeRham.DoubleCover
DeRham.zeta_coefficients(::DeRham.ZetaFrontend)
```
