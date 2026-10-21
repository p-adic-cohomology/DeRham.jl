# Jacobians

Let $C$ be a smooth plane curve over $\mathbb{F}_q$, with genus $g$, and let
$J = \operatorname{Jac}(C)$. The Jacobian's first cohomology agrees with the curve's, so
both have the same $H^1$ L-polynomial
$P_C(T) = \prod_{j=1}^{2g}(1 - \alpha_j T)$. This is the coefficient
vector returned by `jacobian_zeta_coefficients`.

The Jacobian's higher cohomology is formed from exterior powers of
$H^1(C)$. For each $i$, the Frobenius eigenvalues on
$H^i(J) = \wedge^i H^1(C)$ are products of $i$ distinct
$\alpha_j$, so the zeta function is

```math
Z(J,T) = \prod_{i=0}^{2g} P_i(T)^{(-1)^{i+1}},
\qquad
P_i(T) = \prod_{1 \leq j_1 < \cdots < j_i \leq 2g}
(1 - \alpha_{j_1} \cdots \alpha_{j_i} T).
```

The point counts are obtained from the curve's Frobenius eigenvalues:

```math
\#J(\mathbb{F}_{q^k}) = \prod_{j=1}^{2g}(1 - \alpha_j^k).
```

For example, a Fermat quartic in three variables over $\mathbb{F}_5$ defines a
smooth genus-three curve:

```jldoctest
julia> using DeRham, Oscar

julia> f = DeRham.fermat_hypersurface(3, 4, 5);

julia> length(DeRham.jacobian_zeta_coefficients(f))
7

julia> length(DeRham.jacobian_point_counts(f, 2))
2
```

## API

```@docs
DeRham.jacobian_zeta_coefficients
DeRham.jacobian_point_counts
DeRham.jacobian_zeta_function
```
