# Dense Dwork quartic K3 derivation

This page records the original K3 range and derivation from PR #99. See
[Native dense hypersurfaces](native_hypersurfaces.md) for the current supported
ranges, threefold/fourfold extension, and the enlarged arithmetic bounds.
The argument below applies to the original 11–101 K3 subrange.

`DeRham.dense_dwork_zeta(f)` is an opt-in experimental calculator for smooth
quartic K3 surfaces over prime fields with **11 ≤ p ≤ 101**, requiring **all
35 degree-four monomials**. It returns the exact full projective zeta function
in `fraction_field(ZZ[T])`. `DeRham.dense_dwork_coefficients(f)` returns the
22 **descending** integer coefficients of the primitive degree-21 polynomial.
Both are public and unexported. Neither changes the default solver.

```julia
using Pkg
Pkg.build("DeRham")  # once: C++17 compiler (GCC or Clang) required
using DeRham, Oscar
R, xs = polynomial_ring(GF(37), ["x", "y", "z", "w"])
exponents = [[a,b,c,4-a-b-c] for a in 0:4 for b in 0:(4-a) for c in 0:(4-a-b)]
coefficients = [20, 17, 4, 15, 23, 2, 30, 36, 7, 21, 17, 9, 8, 20, 25, 14, 3, 32, 13, 36, 18, 33, 25, 1, 25, 30, 6, 29, 5, 31, 29, 29, 27, 36, 36]
f = sum(R(c)*prod(xs .^ e) for (e,c) in zip(exponents, coefficients))
Z = DeRham.dense_dwork_zeta(f)
```

This is the archived p=37 equation. See `benchmarks/dense_dwork/fixtures.jl`
for all four checked equations. Unsupported inputs and
singular surfaces raise `ArgumentError`. If the bounded precision attempts
cannot establish a unique lift, the call raises an error rather than returning
a guess. A supported prime/shape is not a promise that every surface will
resolve within the native precision limits.

The method is **dense Dwork–Cartier with Harvey window coefficient extraction**.
It uses batched Harvey windows for p ≥ 17 and the workshop's faster half-power
coefficient calculation at p=11,13, with Julia/Oscar validation and arbitrary-precision
integer reconstruction. There is no Python runtime or workshop dependency.
The optional binary is built once, not during solves. One process per matrix
preserves batching and isolates native global state. See
`deps/native/PROVENANCE.md` in the repository for the
process protocol, licensing, source attribution and build limitations.
A JLL artifact is a possible future packaging improvement. The tested backend
is serial; the implementation makes no square-root-in-p complexity claim and
is not controlled reduction.

## Mathematical argument tied to the implementation

Write `f` for the coefficientwise Teichmüller lift. Let `W = πw`,
`π^(p-1) = -p`, and use positive monomials `W^h x^u`, `|u| = 4h`.
The native rows represent images under the projective Cartier operator
`α = ψ_p ∏_a θ(f_a w x^a)`, with `θ(t) = exp(π(t-t^p))`.
They are not asserted to equal an existing DeRham matrix entry by entry.
The relevant normalization is **det(1 − Tα) = P(pT)**, where `P` is the
primitive Frobenius polynomial of the K3 surface.

### Integral quotient and geometric smoothness

`build_level` uses `D_i = x_i∂_i + Wx_i∂_i f`. On a monomial allowed to
omit `x_i`, the relation is

```math
 W^k x^m(x_i\partial_i f) = -m_i W^{k-1}x^m.
```

After removing `x₁x₂x₃x₄`, its leading term is the ordinary Jacobian relation.
Smoothness, Euler's identity and `p ∤ 4` imply that the four partials form a
regular sequence. Its Hilbert series is `(1+t+t²)^4`. Thus the complements
in W-degrees 1,2,3 have dimensions 1,19,1, and degree 4 has no complement.
The native code verifies these ranks over F_p, including surjectivity in
ordinary polynomial degree 12. Surjectivity excludes a common projective
zero over the algebraic closure: at such a point some `x_i^12` would be
nonzero yet lie in the Jacobian ideal. Julia also calls `check_smoothness`.

The degreewise Koszul complex is split exact away from degree zero over F_p.
Lifting its unit pivots gives free Jacobian quotients over Z_p; the filtered
Dwork complex has the same basis. `invert`, `Bas` and `Low` implement these
relations modulo p^N, dividing only by units. Reduction therefore loses no
absolute precision. In the Griffiths–Dwork dictionary the positive monomial
corresponds, up to an overall normalization, to
`(-1)^(h-1)(h-1)! x^(u-1) Ω/f^h`: the factorial has already been absorbed
in the lattice. There is no missing division by the pole order.

### Projective trace, including coordinate boundaries

This paragraph supplies the extension beyond the general-position hypothesis
of the projective specialization in Sperber–Voight §5. Full smoothness does
**not** imply smoothness of every coordinate section for a dense equation.
No such implication is used here.

For `J ⊆ {1,2,3,4}`, let `C_J` allow zero exponents in J and require positive
ones elsewhere. In homological degree j take the direct sum over `|J|=j`,
with differential D_i and chain Frobenius p^j α. Filter by h+j. Dividing by
the product of the required variables identifies the associated graded with
the ordinary Jacobian Koszul complex, which is exact in positive degrees as
above. Its degree-zero homology is the quotient computed by `build_level`.
The splitting estimate below allows a finite truncation at every precision;
unit-pair cancellation equates its alternating trace with the quotient trace.

For q=p^r, let S_I be the exponential sum of wf on the stratum with precisely
the coordinates I nonzero and w nonzero. Teichmüller averaging gives raw trace
`S_I/(q-1)^(|I|+1)` on that closed face. Inclusion–exclusion for required positive
coordinates and the alternating factors `(-q)^j` give

```math
 \operatorname{Tr}(\alpha^r\mid H_0)
 =\frac1{q-1}\sum_I S_I
 =q\bigl(\#X(\mathbf F_q)-(1+q+q^2)\bigr).
```

Indeed the strata partition affine four-space; summing first over nonzero w
counts `q #{f=0 in A⁴} − q⁴`, and the affine cone has `1+(q-1)#X` points.
The empty face is included. The logarithmic determinant identity proves the
normalization for every r. This uses smoothness for the Jacobian complex,
not for the boundary strata. This comparison is a central review obligation;
finite polynomial agreement is separate evidence, not its proof.

### Cutoff and every precision saving used here

[Sperber–Voight (1.4)](https://jvoight.github.io/articles/sparse-dwork-073025.pdf)
gives `v_p(λ_e) ≥ e(p-1)/p²`. After the W rescaling, a coefficient from input
W-degree h to output W-degree K has valuation at least

```math
 L_p(K,h)=\frac{(pK-h)(p-1)}{p^2}+\frac{h-K}{p-1}.
```

This is a bound for the **complete coefficient** of the splitting product,
not for individual terms of an exponential expansion. Its coefficients of
K and h are positive for p ≥ 11. The accepted pairs are N=3,4 for primes
11–101, and additionally N=5 for primes 11–47. For all these pairs,
`L_p(N,1)>N-1` (the smallest margin is 38/605, at p=11,N=5).
The W-scaled coefficients are integral rational p-adic numbers, so all
complete output degrees K ≥ N vanish modulo p^N. This proves `top=N-1`
without assuming extra divisibility after reduction.

For retained K, `compute` sums **all** k+b=K in

```math
 \alpha(W^h x^\mu)=
 \sum_{k,b}\frac{(-1)^b}{b!}(Wf)^b
 \frac{(-p)^k}{(pk-h)!}W^k
 \sum_\nu [x^{p\nu-\mu}]f^{pk-h}x^\nu.
```

Here h=1,2,3 and k≤4<p, so m=pk−h is positive, m<p², and
`v_p(m!)=k-1`. The multiplier supplies **exactly one p**. Consequently
coefficients of f^m modulo p^(N−1) suffice. The remaining b! has b≤3<p.
The lift `a^(p^(N-1)) mod p^N` is the Teichmüller lift to the needed precision.
Using that integer representative in the higher-precision window is safe:
its exact polynomial powers agree with the true lift modulo p^(N−1).
There is no assumption that extra Teichmüller digits have been computed.

### Window recurrence and guarded division

In `harvey_window.hpp`, take `G=Σx_i^4` and `H_j=f^jG^(m-j)`.
This is the deformation recurrence of
[Harvey, *Computing zeta functions of arithmetic schemes*, Theorem 4.1](https://doi.org/10.1112/plms/pdv056).
The product rule gives

```math
 \theta_i(fH_{j-1})-j(\theta_i f)H_{j-1}
 =4(m-j+1)x_i^4H_j.
```

A slot t reads the coefficient at `center−t`. Choose t_i≥4 and put
`base=t−4e_i`; the resulting row weights are exactly
`center_i−base_i−j a_i`, and successor slots are `base+a`.
All slots have degree 13, so such an axis always exists and closure is finite
(at most binomial(16,3)=560 slots). Requests for a common input pole order
share this recurrence. The initial values are the multinomial coefficients
of G^m, computed with factorials split into valuations and units.

Induction scales the output by `4^m m!`. Set e=v_p(m!); the window computes
modulo `q=p^(N-1+e)`, divides the canonical residue exactly by p^e (checking
its divisibility), then multiplies by the inverse of `4^m m!/p^e` modulo
p^(N−1). This division gives a well-defined residue at the required precision.
It does not invert p. Negative/infeasible multinomial indices are zero.

### Small-prime coefficient calculation

At p=11,13, `Pull` builds f^h modulo Q=p^(N−1) by ordinary homogeneous
polynomial multiplication, retaining the half powers needed by the requests.
`Box` extracts `[x^e] f^m` as the convolution of f^floor(m/2) and
f^ceil(m/2). Its recursive bounds enumerate exactly the exponent splits
u+v=e of the required degrees. In a square, the first-coordinate slices u₁
and e₁−u₁ agree; doubling their paired sum and adding the central slice once
preserves the coefficient. All requests still share the power tables within
one native matrix call. No p-adic division is used in this extraction; the
same one-digit multiplier proved above makes precision Q sufficient.
This retained workshop path fixes the measured small-input regression caused
by forcing a window at every prime.

### Native word audit

These bounds cover **every accepted pair**, including the N=5 retry, not only
the four benchmark runs. The standalone executable repeats the prime, degree,
full-support, coefficient and precision guards; environment variables cannot
select an unaudited alternate path.

| Quantity / code | Bound and arithmetic |
| --- | --- |
| Matrix modulus M=p^N | at most 47^5=229345007 < 2^28; residues fit uint32 |
| Canonical modular multiply/add | `(M-1)^2+(M-1)<2^56`; uint64 suffices |
| Extended Euclid `invmod` | remainders and convergent coefficients bounded by M; even the coarse product bound M²+M < 2^57 fits int64 |
| Delayed `invert` sums | at most `floor((UINT64_MAX-(M-1))/(M-1)^2)` updates between reductions; pivot resets only decrease entries |
| Jacobian rank elimination | at most `floor((UINT32_MAX-p)/p²)` unnormalized updates; each adds <p²; initial entries <p; reduction occurs before exceeding that budget |
| `Bas` accumulation | at most `np*35` products; it reduces each addition if that exceeds the same uint64 product budget |
| All other matrix/slab sums | reduce after every addition; operands canonical modulo M |
| Dimensions / indices | ≤455 columns, ≤880 relation rows, rank21; slot count≤560, successor count≤19600; signed int indices/products fit |
| Exponents | m≤302; centers≤918, slot coordinates≤13; the recurrence factor magnitude is <2^11 |
| Guard modulus q | at most 47^7=506623120463 <2^39 (code additionally enforces ≤2^40) |
| Small-prime `Pull` sums | stored residues <13^4<2^15 fit uint16; ≤35 products total <2^36 in uint64 |
| Small-prime `Box` sums | half-power degree≤104, ≤binomial(107,3)=198485 terms; the complete dot product (including square pairing) is <2^48 in uint64 |
| Small-prime table indices | binomial table through 114, entries≤binomial(114,4)=6672876; half-power vectors have <2^18 entries; all integer indices fit |
| Window multiplication / modular power | uint128 product <2^80 under the code's guard |
| Window weights | magnitude <2^11·2^28=2^39, so signed int64 fits |
| Window accumulated dot product | ≤35(q−1)·max_weight <2^85, so signed int128 fits; the code checks the actual absolute-sum bound, using int64 only when it fits |

`choose` uses the small level table for windows, and the bounded half-power
table just described at p=11,13. Factorials use modular multiplication, never an integer
factorial. Julia characteristic polynomials, Newton sums, interval endpoints,
functional equations and zeta assembly all use arbitrary-precision integers.
The test harness also checks maximum-range runs under address/undefined-behavior
sanitizers; finite tests supplement, rather than establish, these bounds.

### Unique lift and full zeta

Let A be the true integral matrix and Ã its residues modulo p^N. Sort measured
column valuation lower bounds `s₁≤…≤s₂₁≤N`. A principal minor of size k with
j erroneous columns has valuation at least
`jN+Σ_{i≤k-j}s_i ≥ N+Σ_{i<k}s_i`. Thus `_dense_dwork_residues` obtains the
coefficient of P modulo `p^max(0,N+Σ_{i<k}s_i-k)` after removing the scale p^k.
It checks divisibility at the available precision. This argument uses observed
column bounds, not an assumed Hodge-adapted basis.

`_dense_dwork_lift` enumerates every prefix allowed by these congruences and
Newton's identity `t_k=−k c_k−Σ_{i<k}c_i t_{k-i}`, with `|t_k|≤21 p^k`.
The self-dual weight-two functional equation supplies the remaining half as
`c_(21-i)=s p^(21−2i)c_i`, trying both determinant signs. All coefficient
congruences and further trace bounds through degree 84 are then checked.
The true polynomial necessarily survives every test. **Exactly one survivor
therefore is the true polynomial**; zero is an inconsistency error; more than
one or exceeding the 4096-prefix work cap requests more precision. Exhausting
the supported precisions raises an explicit error. The finite trace checks
are not being used to prove purity of an arbitrary candidate polynomial.

Finally `zeta_function(coefficients,p,2)` constructs

```math
 Z(X,T)=\frac1{(1-T)(1-pT)(1-p^2T)P(T)}.
```

## Reproducibility and measurements

See `benchmarks/dense_dwork/README.md` in the repository for commands and results.
The four full-polynomial fixtures come from the pinned workshop; those are
regression evidence. Comparisons with upstream controlled reduction and direct
point enumeration are independent evidence, listed separately there.

```@docs
DeRham.dense_dwork_coefficients
DeRham.dense_dwork_zeta
```
