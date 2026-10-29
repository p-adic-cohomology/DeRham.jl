# Artin–Mazur height

The Artin–Mazur height of a K3 surface is the height of its one-dimensional
formal Brauer group. It is an integer from 1 to 10, or infinity. Height 1 is
the ordinary case; infinite height is the supersingular case.

The polynomial form is for a smooth quartic K3 surface in $\mathbb{P}^3$. The data form
takes the 22 coefficients of its primitive $H^2$ L-polynomial over $\mathbb{F}_q$;
it accepts raw `zeta_coefficients` output, `tate_twist` output, or
`normalized_tate_twist` output.

For a quartic in four variables over $\mathbb{F}_p$, the polynomial form calls
`frobenius_matrix_with_precision(f, [1, 1, 1])` and uses the characteristic
polynomial of the returned plain matrix over $\mathbb{Z}/p^M$. It lifts the
matrix entries to integers before computing that polynomial; it does not
compute the full L-polynomial. Inputs outside this case retain the
`zeta_coefficients` path.

Write the entries returned by `coefficients(charpoly(...))` as
$c_1,\ldots,c_{22}$. The height test divides $c_i$ by $p^{21-i}$ and checks
the result modulo $p$ for $i=2,\ldots,11$. Therefore $c_i$ must be known
modulo $p^{22-i}$:

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

The primitive $H^2$ Hodge numbers are $[1,19,1]$: one column has offset $0$,
nineteen have offset $1$, and one has offset $2$. Thus at relative request
$[1,1,1]$, column $j$ is known to at least modulus $p^{1+o_j}$, where $o_j$ is
its Hodge offset. The returned matrix has one common parent ring
$\mathbb{Z}/p^\ell$; its modulus is at least the largest absolute precision,
while some columns contain fewer known digits than that common modulus. These
column guarantees include the final scaling in `FinalReduction.jl`.

To see that the tested coefficients meet their required precision, set
$k=22-i$. The coefficient $c_i$ is the coefficient of $T^{i-1}$ in
$\det(TI-F)$, hence a signed sum of principal minors of size $k$. Mazur
divisibility gives every entry in a column of offset $o_j$ valuation at least
$o_j$. In the difference between a minor of the exact matrix and the minor of
the lifted approximation, each term has at least one error factor from a
column $j$; that factor is divisible by $p^{1+o_j}$, and the other factors
contribute their column offsets. The error is therefore divisible by
$p^{1+\sum_{j\in I}o_j}$ for the selected column set $I$. For
$11\leq k\leq20$, the least possible offset sum among $k$ columns is $k-1$
(the offset-$0$ column and $k-1$ of the nineteen offset-$1$ columns). Hence
$c_i$ is determined modulo $p^k=p^{22-i}$, exactly the modulus required by
the normalized Tate-twist test above. This argument uses only the per-column
precision contract of the plain matrix, not a per-column ring or metadata.

For pole order $m$, let $N_m$ be the certified retained-term count,
$s_m=m+N_m-1$, and $q_m=2-v_p((p s_m-1)!)$. The computation applies the final
column scale $p^{q_m}$ before embedding all columns into one matrix over a
common residue ring. The selector and its all-tail proof are documented in
the [Frobenius matrix precision proof](@ref sec-certified-series-truncation)
and proved in `dr-ine1` / `dr-eevm.1`. For the quartic request $[1,1,1]$,
the certified values for the benchmark characteristics are:

| p | N_1, N_2, N_3 | Algorithm precision M | q_1, q_2, q_3 | Common matrix modulus |
| ---: | ---: | ---: | ---: | ---: |
| 3 | [5, 4, 3] | 6 | [-3, -3, -3] | p^3 |
| 5 | [3, 2, 1] | 3 | [0, 0, 0] | p^3 |
| 7 | [3, 2, 1] | 3 | [0, 0, 0] | p^3 |
| 11 | [3, 2, 1] | 3 | [0, 0, 0] | p^3 |

The F5 quadric regression in `test/Precision.jl` checks that the final
per-column scale preserves the expected residue after division by $p$.

```@docs
DeRham.artin_mazur_height
```
