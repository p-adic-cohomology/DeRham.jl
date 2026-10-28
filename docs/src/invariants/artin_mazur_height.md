# Artin–Mazur height

The Artin–Mazur height of a K3 surface is the height of its one-dimensional
formal Brauer group. It is an integer from 1 to 10, or infinity. Height 1 is
the ordinary case; infinite height is the supersingular case.

The polynomial form is for a smooth quartic K3 surface in $\mathbb{P}^3$. The data form
takes the 22 coefficients of its primitive $H^2$ L-polynomial over $\mathbb{F}_q$;
it accepts raw `zeta_coefficients` output, `tate_twist` output, or
`normalized_tate_twist` output.

For a quartic in four variables over F_p, the polynomial form calls
`frobenius_matrix_with_precision(f, [1, 1, 1])` and uses the characteristic
polynomial of that matrix. It does not compute the full L-polynomial. Inputs
outside this case retain the `zeta_coefficients` path.

The vector returned by `coefficients(charpoly(...))` has the same order as
DeRham's descending L-polynomial coefficient vector. Write its entries as
c_1, ..., c_22. The normalized Tate twist tests c_i / p^(21-i) modulo p for
i = 2, ..., 11. Therefore c_i must be known modulo p^(22-i):

| Index i | Divisor before reduction | Required modulus for c_i |
| ---: | ---: | ---: |
| 2 | p^19 | p^20 |
| 3 | p^18 | p^19 |
| 4 | p^17 | p^18 |
| 5 | p^16 | p^17 |
| 6 | p^15 | p^16 |
| 7 | p^14 | p^15 |
| 8 | p^13 | p^14 |
| 9 | p^12 | p^13 |
| 10 | p^11 | p^12 |
| 11 | p^10 | p^11 |

The primitive H^2 Hodge polygon is [1, 19, 1]. Its coefficient precision
vector at relative request [1, 1, 1] gives absolute precisions
[20, 19, 18, 17, 16, 15, 14, 13, 12, 11] for these indices, exactly the
moduli in the table. Each Frobenius column keeps its own residue ring. For
pole order m, let N_m be the certified series cutoff, s_m = m + N_m - 1, and
q_m = 2 - v_p((p s_m - 1)!). The final column scaling is p^q_m, so that
column is represented modulo p^(M + q_m), with common algorithm precision M.
For the quartic request [1, 1, 1], the certified values for the benchmark
characteristics are:

| p | N_1, N_2, N_3 | M | q_1, q_2, q_3 | Column modulus |
| ---: | ---: | ---: | ---: | ---: |
| 3 | [5, 4, 3] | 6 | [-3, -3, -3] | p^3 |
| 5 | [3, 2, 1] | 3 | [0, 0, 0] | p^3 |
| 7 | [3, 2, 1] | 3 | [0, 0, 0] | p^3 |
| 11 | [3, 2, 1] | 3 | [0, 0, 0] | p^3 |

This scale is included in the Hodge-filtered precision guarantee; the
calculation does not assume every entry has modulus p^M. The all-tail selector
certifies N_m as described in the [Frobenius matrix precision proof](@ref sec-certified-series-truncation).
The F5 quadric regression in
`test/Precision.jl` checks that the final per-column scale preserves the
expected residue after division by p.

```@docs
DeRham.artin_mazur_height
```
