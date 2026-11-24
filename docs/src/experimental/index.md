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
