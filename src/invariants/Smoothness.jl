"""
    check_smoothness(f)

Check whether the homogeneous Oscar polynomial ``f`` defines a smooth
hypersurface over its coefficient field.

``f`` must already be homogeneous; Oscar may throw an error otherwise. The
result is a `Bool`.

This wrapper is named `check_smoothness` because Oscar already provides
`is_smooth`.

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> DeRham.check_smoothness(f)
true
```
"""
function check_smoothness(f)
    p = characteristic(parent(f))
    nVars = length(gens(parent(f)))

    graded, _ = grade(parent(f))

    R, _ = quo(graded, ideal(graded, [graded(f)]))

    V = proj(R)

    is_smooth(V)
end
