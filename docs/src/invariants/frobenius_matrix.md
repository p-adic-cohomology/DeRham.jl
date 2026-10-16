# Frobenius matrix

For a smooth hypersurface over a finite field of characteristic `p`, the
Frobenius matrix represents the `p`-power Frobenius action on its middle
primitive de Rham cohomology. The matrix is computed to a requested
`p`-adic precision. By default, DeRham constructs the cohomology basis
internally with `get_basis_of_cohomology`; that helper and the basis ordering
are internal details. A supplied `basis` selects the basis used by the
computation.

```@docs
DeRham.frobenius_matrix_with_precision
```
