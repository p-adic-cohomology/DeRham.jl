# Experimental algorithms

Algorithms documented in this section have no proof of correctness, and their
precision choices are not certified. Their symbols are public but not exported.

## Support-cone reduction

Projective hypersurfaces can select the experimental CPU policy with
`algorithm=:supportcone` in `DeRham.zeta_coefficients`. The policy builds exact
plans from the exponent vectors in the defining polynomial's support, sorting
those vectors lexicographically. To bound construction cost, it always includes
the lexicographically first full-rank subcone when one exists, searches at most
the first 256 support-index combinations for more candidates, and retains at
most 64. Candidate subcones may overlap and do not form a triangulation. For a
chain key it takes the first retained subcone containing the key, uses floored
Cramer coordinates for chunk lengths, and lowers the largest length first
(lowest index breaks ties) until the plan reaches the terminal pole order.
Planned chunks are checked against the frontend's `V_S(u,k)` predicate. Keys
without a valid retained support plan fall back to the depth-first direction
rule. An impossible fallback direction raises an error before calling the
shared direction chooser.

The support-cone policy preserves the caller's `S`; the default `S=[-1]`
expands to all ambient indices through the ordinary hypersurface setup. The
policy does not support GPU execution, double covers, or custom reduction
contexts. It has no proof of correctness, so compare its results with a proven
policy such as `:depthfirst` for the varieties and primes of interest.

## AKR pole reduction

Projective hypersurfaces accept three experimental AKR algorithm symbols:
`:akr`, `:akr_naive`, and `:akr_levelbylevel`. The first two select the
original method, which reduces each monomial's pole-order chain separately.
`:akr_levelbylevel` combines terms at each pole order before applying the
reduction map. It reuses the original polynomial-vector conversions for
derivatives.

These algorithms are selected with `zeta_coefficients(...; algorithm=...)`.
They are not part of the default reduction heuristic.
