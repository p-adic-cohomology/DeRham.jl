# Utilities

These helpers are public but are not exported. Call them with the `DeRham.`
prefix.

## Hypersurfaces

```@docs
DeRham.fermat_hypersurface
DeRham.random_hypersurface
DeRham.check_smoothness
```

## Slope polygons

`SlopesPolygon` stores the slopes and lengths of a polygon starting at
`(0, 0)`. Its integer-vector constructor uses slopes `0, 1, 2, ...` with the
corresponding entries giving their horizontal lengths. The two-vector
constructor takes points `(xs[i], yvals[i])` and forms their lower convex
hull.

```@docs
DeRham.SlopesPolygon
```

## Tate twists

For a coefficient vector describing
`P(T) = c₀T^d + ⋯ + c_d`, `tate_twist` returns the coefficients of
`P(T/q^j)`, `tate_untwist` reverses that substitution (and an optional
normalization), and `normalized_tate_twist` returns the coefficients of
`q^s P(T/q^j)`. These helpers keep the descending coefficient order.

```@docs
DeRham.tate_twist
DeRham.tate_untwist
DeRham.normalized_tate_twist
```
