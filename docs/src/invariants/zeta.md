# Zeta function

Let $X$ be a smooth projective hypersurface of dimension $m$ over $\mathbb{F}_q$. Its
zeta function records the point counts over finite extensions of $\mathbb{F}_q$:

```math
Z(X,T) = \exp\left(\sum_{k\geq 1} \#X(\mathbb{F}_{q^k})\frac{T^k}{k}\right).
```

## Projective hypersurfaces

For a projective hypersurface,
the middle-cohomology L-polynomial $P_m(T)$ is the only interesting factor.
This factor is computed by
`zeta_coefficients`. With the convention used here,

```math
Z(X,T) = \frac{P_m(T)^{(-1)^{m+1}}}{\prod_{i=0}^{m}(1-q^iT)}.
```

Coefficient vectors are written in descending powers of $T$; the last entry
of the raw vector is the constant coefficient `1`. For example, the output
`[7, 1, 1]` corresponds to the polynomial
```math
7T^2 + T + 1.
```
For the methods below which take as input a homogeneous polynomial $f$, a
nonsmooth hypersurface currently produces `false`.

## API

The relative precision passed to `DeRham.zeta_coefficients_with_precision` is
the requested relative precision of the Frobenius matrix. The intermediate
Frobenius matrix, or the returned matrix when `givefrobmat=true`, is computed to
at least that precision on each Hodge-filtered piece and may use greater
precision when needed for the zeta coefficients to be correct. The certified
series-precision selector is proved in dr-eevm.1, “Certified series precision:
proof of the fixed algorithm.”

```@docs
DeRham.zeta_coefficients
DeRham.zeta_coefficients_with_precision
DeRham.zeta_function
DeRham.L_polynomial
DeRham.point_counts
```
