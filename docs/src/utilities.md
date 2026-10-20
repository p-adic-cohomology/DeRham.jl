# Utilities

These helpers are public but are not exported. Call them with the `DeRham.`
prefix.

## Hypersurfaces

```@docs
DeRham.fermat_hypersurface
DeRham.random_hypersurface
DeRham.check_smoothness
```

## Tate twists

For a coefficient vector describing
$P(T) = c_0 T^d + \cdots + c_d$, `tate_twist` returns the coefficients of
$P(T/q^j)$, `tate_untwist` reverses that substitution (and an optional
normalization), and `normalized_tate_twist` returns the coefficients of
$q^s P(T/q^j)$. These helpers keep the descending coefficient order.

```@docs
DeRham.tate_twist
DeRham.tate_untwist
DeRham.normalized_tate_twist
```
