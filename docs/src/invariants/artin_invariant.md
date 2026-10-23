# Artin invariant

Let $X$ be a supersingular K3 surface over $\mathbb{F}_p$, i.e. one of infinite
[Artin–Mazur height](artin_mazur_height.md). Its geometric Picard rank is
$22$, and the discriminant of the Néron–Severi lattice is $-p^{2\sigma}$ for an
integer $1 \le \sigma \le 10$, the Artin invariant of $X$. By Ogus's theory of
supersingular K3 crystals, $\sigma$ is determined by the Frobenius matrix on
$H^2_{\mathrm{crys}}(X/\mathbb{Z}_p)$.

`artin_invariant` computes $\sigma$ for a smooth quartic K3 surface in
$\mathbb{P}^3$ from its $21 \times 21$ Frobenius matrix $F$ on primitive $H^2$,
using only $F \bmod p^2$ and linear algebra over $\mathbb{F}_p$. The polynomial form computes
$F$ with [`frobenius_matrix_with_precision`](@ref) at relative precision
`[0, 1, 2]`, which is absolute precision $p^2$ in every column; see the
[Frobenius matrix](frobenius_matrix.md) page for the column ordering and the
divisibility it relies on. The matrix form accepts any such $F$ over
$\mathbb{Z}/p^M$ with $M \ge 2$.

Supersingularity is a precondition and is not checked; check it first with
[`artin_mazur_height`](@ref), which returns `Inf` exactly in the supersingular
case. The prime $p$ must be odd.

## Example

The quartic below is supersingular over $\mathbb{F}_3$ with Artin invariant
$2$ (one of the examples in Table 2 of Takamatsu–Yoshikawa, arXiv:2605.06252):

```jldoctest
julia> R, (x, y, z, w) = polynomial_ring(GF(3), ["x", "y", "z", "w"]);

julia> f = x^4 + y^3*z + y*z^3 + x*y^2*w + w^4;

julia> DeRham.artin_invariant(f)
2
```

## Validation

The values computed by `artin_invariant` have been cross-checked against the
independent Macaulay2 computation of Takamatsu–Yoshikawa (arXiv:2605.06252)
on supersingular quartics at $p = 3, 5, 7, 11$.

## Reference

- A. Ogus, *Supersingular K3 crystals*, Astérisque 64 (1979), 3–86.

```@docs
DeRham.artin_invariant
```
