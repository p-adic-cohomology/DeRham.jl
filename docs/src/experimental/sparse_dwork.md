# Sparse projective hypersurfaces

`DeRham.sparse_dwork_zeta(f)` computes the projective zeta function of a
smooth cubic or quartic in three or four variables over a prime field of
characteristic at least five. The equation must have all pure powers and at
most one mixed monomial. This CPU implementation uses Julia integer
arithmetic and Oscar. It is opt-in.

```julia
R, (x, y, z, w) = polynomial_ring(GF(7), [:x, :y, :z, :w])
f = x^4 + y^4 + z^4 + w^4 + 2*x*y*z*w
cs = DeRham.sparse_dwork_coefficients(f)    # all 22 coefficients, constant 1 last
Z = DeRham.sparse_dwork_zeta(f)
```

## Scope and correctness

We use the projective complex of [Sperber and Voight, §§1, 2, 5](https://jvoight.github.io/articles/sparse-dwork-073025.pdf).
The ordinary partial derivatives certify geometric smoothness: for a smooth
degree-$d$ form in $n$ variables with $p\nmid d$, their quotient has socle
degree $n(d-2)$. Surjectivity in degree $n(d-2)+1$ excludes every common
projective zero over the algebraic closure. We test this by unit pivots
modulo $p$. Under the stated support hypothesis every coordinate section is
also smooth. If setting some coordinates to zero removes the mixed term,
the section is diagonal. Otherwise each omitted variable occurs only as a
pure power and its derivative vanishes on the section. A singular point of
that section would therefore be singular on the full hypersurface.

For the splitting function write $\ell_e=[t^e]\exp(t+t^p/p)$ and
$F_e=v_p(e!)$. The Cartier kernel computes $p^{F_e}\ell_e$ by the integral
recurrence obtained from $e\ell_e=\ell_{e-1}+\ell_{e-p}$, with negative
indices zero. The weight bound (1.4) implies that a term with
$e(p^2-3p+1)\ge Np^2(p-1)$ may be omitted modulo $p^N$. The exact residue
conditions select the monomials surviving Cartier. The coefficient lifts
are Teichmüller representatives. After writing $W=\pi w$, reduction uses
Jacobian relations and differentiation with integral coefficients, so the
raw cutoff remains valid after reduction.

The resulting matrix $A$ satisfies $\det(1-TA)=P(pT)$ for the primitive
middle factor $P$. To reconstruct $P=1+\sum c_k T^k$, sort the minimum
entry valuations of the computed columns, capped at $N$, as
$v_1\le\cdots\le v_r$. In each principal minor of order $k$, changing a
column by $p^N$ changes the determinant by a multiple of
$p^{N+v_1+\cdots+v_{k-1}}$. This also covers terms with several changed
columns because $N\ge v_j$. After dividing by $p^k$ we have a certified
residue for $c_k$. No divisibility by pole degree is assumed.

Given the preceding coefficients, Newton's identity determines the $k$th
power sum as a known integer minus $k c_k$. Its absolute value is at most
$r p^{(n-2)k/2}$. The implementation uses integer square roots and accepts
a coefficient only if exactly one lift satisfies this bound. An ambiguous
lift triggers a new calculation with two more digits. For fixed rank and
prime the modulus grows without bound, so this criterion eventually decides
every coefficient. DeRham's existing coefficient-based zeta constructor then
adds the projective ambient factors. The algorithm does not infer the answer
from the equation's membership in a special pencil.

```@docs
DeRham.sparse_dwork_coefficients
DeRham.sparse_dwork_zeta
```
