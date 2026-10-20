# Zeta function

Let `X` be a smooth projective hypersurface of dimension `m` over `F_q`. Its
zeta function records the point counts over finite extensions of `F_q`:

```math
Z(X,T) = \exp\left(\sum_{k\geq 1} \#X(\mathbb{F}_{q^k})\frac{T^k}{k}\right).
```

## Projective hypersurfaces


For a projective hypersurface, 
the middle-cohomology L-polynomial `P_m(T)` is the only interesting factor. This factor is computed by
`zeta_coefficients`. With the convention used here,

```math
Z(X,T) = \frac{P_m(T)^{(-1)^{m+1}}}{\prod_{i=0}^{m}(1-q^iT)}.
```

Coefficient vectors are written in descending powers of `T`; the last entry
of the raw vector is the constant coefficient `1`. For example, the output 
`[7, 1, 1]` corresponds to the polynomial 
```math
7t^2 + t + 1.
```
For the methods below which take as input a homogenous polynomial f, a nonsmooth hypersurface currently produces `false`.

## API

```@docs
DeRham.zeta_coefficients
DeRham.zeta_coefficients_with_precision
DeRham.zeta_function
DeRham.L_polynomial
DeRham.point_counts
```
