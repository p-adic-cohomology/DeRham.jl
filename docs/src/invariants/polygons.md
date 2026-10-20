# Newton & Hodge polygons

For a smooth hypersurface over a finite field, the Newton polygon records the
$p$-adic valuations of the Frobenius eigenvalues on primitive middle
cohomology. If its L-polynomial is
$P(T) = c_d T^d + \cdots + c_0 = \prod_i (1 - \alpha_i T)$, the Newton
polygon is the lower convex hull of the coefficient points
$(i, v_p(c_i))$; its slopes are the valuations of the reciprocal roots
$\alpha_i$. Coefficient vectors passed to DeRham are in descending order.

The Hodge polygon records the Hodge numbers of that cohomology. For a smooth
degree-$d$ hypersurface in $\mathbb{P}^n$, the horizontal length at slope $r$ is
$h^{n-1-r,r}$, for $r = 0, \ldots, n-1$. Both functions return a
[`SlopesPolygon`](@ref DeRham.SlopesPolygon), which stores the slopes and
their horizontal lengths starting from $(0, 0)$.

```@docs
DeRham.newton_polygon
DeRham.hodge_polygon
```

## The SlopesPolygon type

`SlopesPolygon` stores the slopes and lengths of a polygon starting at
$(0, 0)$. Its integer-vector constructor uses slopes $0, 1, 2, \ldots$ with the
corresponding entries giving their horizontal lengths. The two-vector
constructor takes points `(xs[i], yvals[i])` and forms their lower convex
hull. `SlopesPolygon` is public but not exported, so call it as
`DeRham.SlopesPolygon`.

```@docs
DeRham.SlopesPolygon
```
