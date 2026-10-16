# Picard rank

For a smooth projective surface, the Picard rank is the rank of its group of
divisor classes modulo algebraic equivalence. The geometric Picard rank is
computed over an algebraic closure; the Picard rank over `F_{q^k}` counts
classes defined over that finite field.

For a surface in `P^3`, the hyperplane class is always algebraic. The
primitive `H^2` L-polynomial omits this class, so the ranks returned by the
functions below include it.

`picard_rank_bound` gives an upper bound for the geometric Picard rank; for a
general surface this bound need not be exact. The three K3-specific functions
give exact geometric or finite-field ranks because the Tate conjecture is
known for K3 surfaces over finite fields. They apply to quartic K3 surfaces
in `P^3`: the polynomial form requires a quartic in four variables, and the
data form requires 22 primitive-`H^2` L-polynomial coefficients over `F_q`.
`k3_picard_realization_degree` is the least extension degree over which the
finite-field Picard rank reaches the geometric Picard rank.

```@docs
DeRham.picard_rank_bound
DeRham.k3_geometric_picard_rank
DeRham.k3_picard_rank
DeRham.k3_picard_realization_degree
```
