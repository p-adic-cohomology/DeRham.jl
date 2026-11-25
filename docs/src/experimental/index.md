# Experimental algorithms

Algorithms documented in this section have no proof of correctness, and their
precision choices are not certified. Their symbols are public but not exported.

## AKR pole reduction

Projective hypersurfaces accept three experimental AKR algorithm symbols:
`:akr`, `:akr_naive`, and `:akr_levelbylevel`. The first two select the
original method, which reduces each monomial's pole-order chain separately.
`:akr_levelbylevel` combines terms at each pole order before applying the
reduction map. It reuses the original polynomial-vector conversions for
derivatives.

These algorithms are selected with `zeta_coefficients(...; algorithm=...)`.
They are not part of the default reduction heuristic.

## Dense reduction handoff

Projective hypersurfaces also accept the experimental keywords
`dense_reduction` and `handoff_pole_order`. The controlled policy reduces to the
handoff order, then the selected AKR dense method reduces to the terminal pole
order. The existing precomputed final reduction maps that result to the
cohomology basis. The handoff order must be at least the terminal pole order
and no greater than the maximum pole order among the generated Frobenius source
terms.

`dense_reduction` accepts `:akr_naive` or `:akr_levelbylevel`; `:akr` is an
alias for `:akr_naive`. If the method is omitted, the existing precomputed
final reduction is used directly.

| `algorithm` | Controlled policy | Handoff | Dense method | Final step |
| --- | --- | --- | --- | --- |
| `:default` | Existing automatic policy | Terminal order | None | Precomputed |
| `:pchunk`, `:depthfirst`, `:varbyvar` | Selected policy | Terminal order | None | Precomputed |
| `:akr`, `:akr_naive` | None | Top source order | `:akr_naive` | Precomputed |
| `:akr_levelbylevel` | None | Top source order | `:akr_levelbylevel` | Precomputed |

For example, `zeta_coefficients(f; algorithm=:depthfirst,
dense_reduction=:akr_levelbylevel, handoff_pole_order=4)` runs controlled
reduction to pole order 4, followed by level-by-level AKR reduction. These
methods have no proof of correctness and their precision choices are not
certified.
