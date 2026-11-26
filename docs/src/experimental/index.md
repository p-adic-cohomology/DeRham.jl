# Experimental algorithms

Algorithms documented in this section have no proof of correctness, and their
precision choices are not certified. Their symbols are public but not exported.

## AKR pole reduction

Projective hypersurfaces accept three experimental AKR algorithm symbols:
`:akr`, `:akr_naive`, and `:akr_levelbylevel`. The first two select the
original method, which reduces each monomial's pole-order chain separately.
`:akr_levelbylevel` combines terms at each pole order before applying the
reduction map. Both schedules use coefficient maps for derivatives when they
use the cutoff table and UInt64 arithmetic. The BigInt and pseudo-inverse
paths differentiate polynomials.

These algorithms are selected with `zeta_coefficients(...; algorithm=...)`.
They are not part of the default reduction heuristic.

When `S` contains all variables and the characteristic does not divide the
hypersurface degree, AKR builds one Jacobian relation table at degree
`(n + 1)(d - 2) + 1` and obtains higher-degree relations by multiplying those
cutoff relations by monomials. For a proper variable subset or a characteristic
dividing the degree, AKR keeps the per-pole-order pseudo-inverse path.
