# Native dense hypersurfaces

The public, unexported functions `DeRham.dense_dwork_coefficients(f)` and
`DeRham.dense_dwork_zeta(f)` now accept these fully dense smooth hypersurfaces:

| Shape | Variables / degree | Primes | Primitive rank | Matrix digits |
| --- | --- | --- | --- | --- |
| Quartic K3 | 4 / 4 | 11–251 | 21 | 3, then 4, then 5 only for p≤47 |
| Cubic threefold | 5 / 3 | 11–23 | 10 | 5 |
| Cubic fourfold | 6 / 3 | 11–167 | 22 | 5 |

All degree-d monomials must have nonzero coefficients. The coefficient field
must be prime. Julia checks geometric smoothness with Oscar. The coefficient
method returns descending integers; the zeta method returns the exact full
projective rational function using DeRham's existing result type. Singular,
unsupported and ambiguous inputs fail explicitly. Precision is computed by the
wrapper; there is no user-supplied uncertified precision or automatic dispatch.
`ProjectiveHypersurface` is also accepted.

Run `Pkg.build("DeRham")` once with GCC/Clang C++17 and 128-bit integer support.
Installation creates one serial executable per variable count. Each call sends
one complete matrix request, preserving all coefficient batching. Solves require
neither Python nor a compiler nor the workshop. There is no shared mutable state
between processes. Prebuilt JLL packaging remains future work.

This extends the K3 integration with selected direct-kernel improvements from
workshop revision `a5f802454f425ed6eae2aea197a6b3fee65ddba8`: shared fourfold
windows, grouped small-prime residue sums, an integral divided matrix, and a
socle smoothness certificate. It uses the portable scalar recurrence. AVX-512 panels, Picard–Fuchs transport, toric
controlled reduction and double-plane kernels are separate routes and are not
part of this implementation. In particular a degree-2 K3 is not a quartic K3.
The method is **dense Dwork–Cartier with Harvey window coefficient extraction**;
no square-root-in-p runtime is claimed.

## Runnable example

This archived smooth cubic threefold over F₁₃ has full monomial support. The
coefficient order below is lexicographic in the exponent tuples.

```julia
using DeRham, Oscar
R, x = polynomial_ring(GF(13), ["x$i" for i in 1:5])
compositions(t, n) = n == 1 ? [[t]] :
    [vcat(a, e) for a in 0:t for e in compositions(t-a, n-1)]
c = [5,5,11,11,3,11,4,11,3,4,11,12,3,3,2,9,4,12,5,1,7,3,11,10,1,5,3,2,5,8,12,7,3,5,6]
f = sum(R(a)*prod(x .^ e) for (a, e) in zip(c, compositions(3, 5)))
Z = DeRham.dense_dwork_zeta(f)  # exact full projective rational function
```

The primitive polynomial has degree 10. It occurs in the numerator of Z;
the denominator is `(1-T)(1-13T)(1-13^2T)(1-13^3T)`.

## Matrix and trace normalization

The [quartic derivation](dense_dwork.md) describes the same positive-monomial
Dwork complex. In n variables and degree d its associated graded Jacobian
quotient has Hilbert series `(1+t+…+t^(d−2))^n`. At W-degree h its ordinary
polynomial degree is dh−n. The three shapes have respective graded dimensions
`(1,19,1)` at h=1,2,3; `(5,5)` at h=2,3; and `(1,20,1)` at h=2,3,4.
The same unit-pivot reduction is integral, with no pole-order division omitted.

The complex allowing zero exponents in a subset J has chain Frobenius
`p^|J| α`. For q=p^r, alternating its coordinate-face traces gives

```math
 \operatorname{Tr}(\alpha^r)
 =(-1)^n q\left(\#X(\mathbf F_q)-\sum_{i=0}^{n-2}q^i\right).
```

This is q times the primitive middle trace by the projective trace formula,
including its odd-dimensional sign. Hence `det(1−Tα)=P(pT)` for all three
shapes. The argument uses smoothness of X in the ordinary Jacobian complex;
it does not infer general position of its coordinate sections. That comparison
remains the central mathematical review point of this experimental integration.

## Conservative cutoff, including the wider range

For input degree h and output degree K the complete splitting coefficient has
valuation at least

```math
 L_p(K,h)=\frac{(pK-h)(p-1)}{p^2}+\frac{h-K}{p-1}.
```

It is increasing in K and h. For every accepted pair (p,N),
`L_p(N,h_min)>N−1`, with h_min=1 for the K3 and 2 for the cubics. Its rational
p-adic coefficient therefore has integral valuation at least N. Thus **all
complete degrees K≥N vanish modulo p^N**, before integral reduction. The short
unit tests check this exact rational inequality for the entire accepted range.
No sharper empirical reduced-slab bound is used, and primes below 11 are refused.

Retained k≤4 and h≤4 satisfy h<p and m=pk−h<p². Therefore
`v_p(m!)=k−1`: the multiplier `(-p)^k/m!` supplies exactly one p. Coefficients
of f^m modulo Q=p^(N−1) suffice. The remaining exponential factorials are units.
For the window recurrence the largest guard is `p^(N−1+v_p(m!))`.
The scaled recurrence computes `d^m m!` times the coefficient at this modulus;
checked division by its p-part and inversion of its unit part give precisely
Q digits, without inverting p. The proof and scalar recurrence are as in the
quartic derivation and Harvey's deformation recurrence.

## Fourfold changes tied to the code

**Divided matrix.** Every retained raw coefficient contains p, and every
reduction is integral. The fourfold computes α/p modulo p^4 throughout, using
a Teichmüller lift to that precision, and restores p on output. Changing the
lift or the unit-pivot reduction by p^4 changes α/p by p^4 and hence α by p^5.
The complete cutoff still uses absolute precision N=5, not N=4. This lowers
the working modulus while preserving all five output digits. Output entries
use 64-bit words because p^5 need not fit 32 bits.

**Shared readouts.** For input poles h=2,3,4, compute f^(pk−4) once and convolve
its readouts with f^(4−h). A readout with input exponent μ and exponent e from
the latter power asks for `[x^(pν−μ−e)]f^(pk−4)`. The slot μ+e lies in the full
positive level-4 monomial space. Summing these exact products modulo Q gives
the same requested coefficient as separate windows. There are at most 462
products per readout; accumulation uses unsigned 128-bit integers. The library
still advances the recurrence linearly in m; sharing is a constant-factor change.

**Small-prime residue grouping.** At p=11,13 write a monomial in the first
half power as e=pU+r with 0≤r_i<p. For the target pν−μ, its complementary
monomial is pV+(p−1−r)−(μ−1), with U+V=ν−1. The packed arrays enumerate
exactly these complementary residue classes; reversing their lexicographic
order pairs r with p−1−r. Negative exponents give zero. Each exponent split
occurs once, so the complete sum is bounded by the smaller half-power table
size times (Q−1)². `residue_coefficients` uses this packing only for the highest
fourfold slab, retaining ordinary convolution for the smaller slabs.

**Socle certificate.** The leading Jacobian quotient in degree 6 is a line
spanned by b. Write `[u x_i]=L(u x_i)b` for each degree-5 monomial u. Commutativity
gives `L(u x_i)x_j b=L(u x_j)x_i b`. Two nonproportional vectors
`(L(u x_i))_i` force every x_i b to vanish. Thus degree 7 and every higher
piece vanish, excluding a common projective zero of the partials. The native
`zero_above_line` tests a nonzero 2×2 minor modulo p; an inconclusive certificate
falls back to the full top-degree rank test. Julia independently checks smoothness.

## Native arithmetic bounds

These supersede the smaller numerical bounds in the original K3 audit, while
retaining the same modular reduction discipline. Every executable repeats its
shape, prime, support, coefficient and precision checks. Environment variables
cannot select extra algorithms or enlarge the range.

| Operation | Bound / storage |
| --- | --- |
| Working modulus M | at most 251^4=3969126001 <2^32; fourfold uses p^4 |
| One modular product plus a canonical residue | `(M−1)^2+(M−1)<2^64`; unsigned 64 bits |
| Delayed Gauss–Jordan / basis sums | explicit budget `floor((UINT64_MAX−(M−1))/(M−1)^2)`; reduce before exhausting it |
| Extended Euclid | signed 128-bit coefficient intermediates; even the coarse M²+M bound fits; remainder products ≤M |
| Rank elimination modulo p | budget `floor((UINT32_MAX−p)/p²)` with canonical pivots |
| All other slab/reduction sums | reduce each addition; canonical operands |
| Fourfold shared readouts | ≤462 products <2^64 each; sum <2^73 in unsigned 128 bits |
| Guard modulus q | at most 167^7 <2^52; enforced before multiplication |
| Window multiplication | unsigned 128-bit product <2^104 |
| Recurrence factors and weights | factor magnitude <2^12; coefficient <2^32; signed weight <2^44 |
| Window dot product | ≤56 terms, each <2^96; total <2^102 in signed 128 bits; actual sum bound checked at runtime |
| Window indices | ≤binomial(18,5)=8568 slots and 56 successors each; signed indices fit |
| Jacobian dimensions | at most 2002 columns and 4752 relation rows, including the fourfold fallback |
| Restored output | <167^5<2^37 in unsigned 64 bits |

Half-power convolution is used for K3s and fourfolds at p=11,13, and for all
accepted cubic threefolds. Every raw power coefficient is a canonical residue
modulo Q=p^(N−1). The complete multiplication sums contain at most the support
size products; complete convolution sums contain at most the stated table size.

| Half-power path | Largest Q | Largest half table | Multiplication sum | Full convolution sum |
| --- | ---: | ---: | --- | --- |
| Quartic K3 | 13^4 | binomial(107,3)=198485 | <2^36 | <2^48 |
| Cubic threefold | 23^4 | binomial(139,4)=14891626 | <2^42 | <2^61 |
| Small cubic fourfold | 13^4 | binomial(80,5)=24040016 | <2^36 | <2^55 |

Threefold power tables use 16-bit words at p=11,13 and 32-bit words at
p=17,19,23; the other small-prime tables use 16-bit words. All complete sums fit unsigned 64 bits; the residue-grouped
fourfold path additionally checks the product budget. The largest binomial
table index is 148, and its entries (through choose(148,5)) fit 32 bits;
all half-table indices fit signed 32 bits. Memory use of the small-prime
fourfold tables is measured separately from the Julia host.
Characteristic polynomials, reconstruction, trace bounds and zeta assembly use
arbitrary precision in Julia. Sanitizers and matrix comparisons supplement
these bounds; they do not establish them.

## Integer reconstruction and ambient factors

The observed column valuation lower bounds give the same principal-minor
precision `N+sum(s_i,i<k)` as in the quartic proof. Removing p^k produces a
certified congruence for the k-th coefficient of P. For dimension m=n−2 and
primitive rank r, the exact trace bound is `t_k²≤r² p^(mk)`; integer interval
endpoints use `isqrt(r² p^(mk))`, avoiding floating-point square roots for
threefolds. Every possible coefficient prefix survives enumeration.

The functional equation is `c_(r−i)=s p^(m(r−2i)/2)c_i`. In odd dimension the
alternating pairing forces s=+1; in even dimension both signs are tried.
For even rank the middle coefficient is checked as well, so s=−1 forces it
to vanish. All coefficient congruences and further trace bounds are checked.
Exactly one candidate certifies the answer. More than one, or exhaustion of
the prefix cap, requests another allowed precision; at the limit it is an
explicit error. Finite trace inequalities alone are not asserted to prove purity.

Finally the existing DeRham assembly returns

```math
 Z(X,T)=\frac{P(T)^{(-1)^{m+1}}}{\prod_{i=0}^m(1-p^iT)}.
```

Thus the primitive factor is in the numerator for the cubic threefold and the
denominator for the K3 and fourfold. Tests compare complete rational functions,
not just their first trace. Reproduction commands, recorded polynomials and
benchmarks are in `benchmarks/native_hypersurfaces/`.
