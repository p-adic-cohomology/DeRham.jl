# PR #97: sparse Dwork reduction, precision, and dependence on p

Review target: commit 21f0ba6309e0649b4b710cb22943815b9d4461d6.
Assignment: dr-flnzq. Date: 2026-11-23.

## Invariant result metadata

~~~json
{
  "assignment": "dr-flnzq: answer the five mathematical questions about PR #97 at head 21f0ba6309e0649b4b710cb22943815b9d4461d6, with derivations and paper/code references; no package-code changes or package loading",
  "result_status": "mathematical report complete and self-audit passed; requested canonical report copy remains a publication-location obligation, dr-q37iq",
  "epistemic_status": "informal proofs and source-level semantic analysis, conditional on the cited Dwork cohomology theorem; exact rational sanity checks are separate computational evidence; no formal verification or Julia execution",
  "assumptions": [
    "n is the number of homogeneous variables, so X lies in P^(n-1), unlike the thesis convention",
    "p is prime, p >= 5, p does not divide d, and N >= 1 is absolute p-adic precision",
    "the public API restricts n to 3 or 4, d to 3 or 4, prime-field coefficients, every nonzero pure-power term, and at most one additional monomial",
    "geometric smoothness is certified by the Jacobian test analyzed in this report",
    "the projective/general-position Dwork comparison and characteristic-polynomial identity from Sperber-Voight section 5 are imported mathematical results"
  ],
  "dependencies": [
    "PR head 21f0ba6309e0649b4b710cb22943815b9d4461d6 and the exact-head code links below",
    "Sperber-Voight, sparse-dwork-073025.pdf, equations (1.4), (1.9), (2.13)-(2.16), Remark 2.17, section 4, and projective section 5",
    "authorized thesis draft contents/cr/main.tex, lines 913-1104",
    "arXiv:2602.24155v1, section 3.4"
  ],
  "checks_performed": [
    "read the source bead, its empty comments, formula contract, exact-head implementation and documentation",
    "derived the reduction identity, factorial dictionary, pole bound, Jacobian cutoff, integrality, Teichmueller lifts, minor precision bound, and Newton/Weil reconstruction",
    "Python 3.9.6 exact rational checks: 735 splitting coefficients, including 525 retained units, for p=5,7,11,13 and N=1,...,5",
    "six exact retained-tuple checks for cubic and quartic pure-power-plus-product supports, p=5,7,11 and N=3; an exact finite-field table transcription also exhibits signed displacement (4,-1,0)",
    "no Julia package loading, package tests, benchmarks, or E2E changes"
  ],
  "findings": [
    "Dwork reduction is integral because the pole factorial is incorporated into the basis; the finite Jacobian table is sufficient",
    "retained output pole k is strictly less than beta(p) N - (beta(p)-1) h/p, with beta(p)<=20/11 for p>=5",
    "the implemented total cost is soft-O(p) for fixed n,d,N in its supported scope; reduction itself has no positive power of p in its operation count",
    "the monomial reducer shares the controlled-reduction mechanism but is not literally the thesis's fixed-block nonnegative-direction game",
    "the six specifically requested correctness concerns have consistent mathematical justifications in the supported scope; documentation and definition issues are identified, rather than a demonstrated program counterexample"
  ],
  "open_obligations": [
    "dr-543md: add factorial normalization, post-Cartier pole bounds, graded freeness, and guard-precision derivations to PR exposition",
    "dr-h5a8g: reconcile the thesis/arXiv V_S support predicate with V_full and the implementation",
    "dr-q37iq: have an authorized owner place the report at the requested canonical-checkout destination, or explicitly resolve the output-file exception"
  ]
}
~~~

The missing per-direction matrices are not a missing mathematical ingredient.
The PR uses finite Jacobian decompositions directly. Its Cartier expansion
already has bounded output degree at fixed precision. The present implementation
is nevertheless linear in \(p\), up to logarithms, because of its tables and,
with one mixed term, its residue enumeration. Accelerating a long pole-reduction
chain cannot remove those costs.

## Sources, conventions, and exact scope

All code line numbers below refer to the review head, not this worktree's older
checkout:

| Abbreviation | Exact source |
| --- | --- |
| R | [SparseDworkReduction.jl][R] |
| C | [SparseCartier.jl][C] |
| I | [invariants/SparseDwork.jl][I] |
| D | [experimental/sparse_dwork.md][D] |
| CR | [reduction/ComputeRuv.jl][CR] |

The external theorem inputs are the coefficient estimate (1.4), the Dwork
operators (1.9), the cohomological reduction (2.15)-(2.16), and the projective
comparison \(\det(1-T\alpha)=P(pT)\) in
[Sperber–Voight, §§1,2,5][SV]. Its complexity claims are Theorems A and B and
§4, especially (4.7). Below are derivations for the PR's particular loops
and representation.

Write
\[
 f=\sum_{j=1}^s \widehat a_j x^{\nu_j},\qquad
 S=(\nu_1\ \cdots\ \nu_s),\qquad |\nu_j|=d,\qquad Q=x_1\cdots x_n.
\]
The coefficients are Teichmüller lifts. The public domain is the prime-field
cubic/quartic domain at I:82-95, with the support restriction at I:53-56.
Statements about larger \(n,d,s\) are local analyses of the internal mathematical
constructions, not extensions of the public API.

The positive exponent \(b\) represents \(W^h x^b\), with \(h=|b|/d\).
Its ordinary numerator is \(g=x^{b-\mathbf1}\), of degree \(hd-n\).
Thus R:1 is shorthand after factoring out \(Q\): the actual projective-complex
monomial is \(QW^h g\). The thesis's ambient dimension is our \(n-1\).
Its target norm \(N_{\rm target}\) is unrelated to precision \(N\).

## 1. Pole reduction and the Griffiths–Dwork dictionary

Take \(\pi^{p-1}=-p\), \(W=\pi w\), and keep \(W\) independent when
differentiating in \(x_i\). The logarithmic Dwork operator is
\[
 D_i=x_i\partial_i+W x_i f_i,\qquad f_i=\partial_i f.                 \tag{A}
\]
The projective complex allows a preimage divisible by every variable except
\(x_i\). For homogeneous \(H_i\) of degree \((k-1)d-n+1\), applying (A) to
\(W^{k-1}(Q/x_i)H_i\) gives
\[
 QW^{k-1}\partial_i H_i+QW^kH_i f_i.
\]
This is zero in the top quotient. Suppress \(Q\) to obtain
\[
 \boxed{W^kH_i f_i\equiv-W^{k-1}\partial_iH_i.}                     \tag{B}
\]
This is the ordinary-partial version of (2.16) in the projective complex of
[SV §5][SV-projective]. In particular differentiating \(Q/x_i\) contributes
zero, which explains the absence of an extra \(+1\).

If the stored Jacobian decomposition is
\[
 x^\eta=r_\eta+\sum_{i,g}c_{\eta,i,g}x^g f_i,                       \tag{C}
\]
then, writing \(u=\xi+\eta\), (B) gives
\[
 W^k x^u\equiv W^k x^\xi r_\eta-
 W^{k-1}\sum_{i,g}c_{\eta,i,g}(\xi_i+g_i)x^{\xi+g-e_i}.             \tag{D}
\]
R:183-203 implements this identity, including the minus sign and
\(\xi_i+g_i\). A zero derivative is discarded before constructing an
exponent that could be negative (R:200-202). Each child has numerator
degree \(kd-n-d\), hence pole \(k-1\).

### Where the denominator went

Define
\[
 \Omega=\sum_{i=1}^n(-1)^{i-1}x_i\,
       dx_1\wedge\cdots\widehat{dx_i}\cdots\wedge dx_n.
\]
Normalize the overall comparison sign by sending a pole-one class to
\(g\Omega/f\). The dictionary on finite-polynomial representatives is
\[
 \boxed{[W^k x^b]\longmapsto(-1)^{k-1}(k-1)!
       \left[\frac{x^{b-\mathbf1}\Omega}{f^k}\right].}              \tag{E}
\]
Indeed the Griffiths–Dwork relation, for \(k>1\), is
\[
 \left[\frac{H_i f_i\Omega}{f^k}\right]
   =\frac1{k-1}\left[\frac{\partial_iH_i\Omega}{f^{k-1}}\right].
                                                                    \tag{F}
\]
For \(c_k=(-1)^{k-1}(k-1)!\), one has
\(c_k/(k-1)=-c_{k-1}\). Multiplying (F) by \(c_k\) gives precisely (B).
The relation \(W^{k+1}fg\equiv-kW^kg\), from \(D_0\), likewise agrees with
cancellation of \(f\) under (E). These computations fix the factorial and sign.
This comparison is for finite representatives and their cohomology classes;
it does not map arbitrary infinite series termwise to rational differentials.

The inverse dictionary divides by \((k-1)!\), up to sign, whose valuation is
\[
 v_p((k-1)!)=\sum_{\ell\ge1}\left\lfloor(k-1)/p^\ell\right\rfloor.  \tag{G}
\]
Denominators have moved into the lattice choice and the splitting coefficients.
Dwork reduction preserves precision in its own integral lattice, while
conversion back to an arbitrary high-pole unscaled Griffiths–Dwork
representative can lose the digits in (G).

For the cohomology basis,
\[
 \lceil n/d\rceil\le h\le H:=\lfloor n(d-1)/d\rfloor.               \tag{H}
\]
In the public domain \(H\le3<p\), so factorials of the basis elements are
units. High intermediate poles need not have this property. The existing
DeRham code already handles factorial rescaling:
[ZetaFrontend.jl:18-24][prescale] and
[CohomologyQuotient.jl:47-51][scalar] give
\(\texttt{reduction_prescale}(m)=(m-1)!\), and
[FinalReduction.jl:149-175][final-scale] separates its valuation and unit.

Equation (E) is not an assertion that the PR matrix equals \(p\) times an
existing DeRham matrix in the same ordered basis. Comparing matrices also
requires fixing Cartier versus pullback, duality, and the residue twist.
The needed normalization is the characteristic-polynomial identity
\(\det(1-TA)=P(pT)\) from [SV §5][SV-projective], used at I:1-3 and D:39.

## 2. Maximum pole order and degree after Cartier

Put
\[
 \ell_e=[t^e]\exp(t+t^p/p),\quad a_p=\frac{2p-1}{p^2(p-1)},\quad
 c_p=\frac1p-a_p=\frac{p^2-3p+1}{p^2(p-1)},\quad
 \beta(p)=\frac1{pc_p}=\frac{p(p-1)}{p^2-3p+1}.                    \tag{I}
\]
As \(\lambda_e=\pi^e\ell_e\), (1.4) implies
\(v_p(\ell_e)\ge-a_pe\), as in (2.13). Thus the integer weight satisfies
\[
 w_e=e+p\,v_p(\ell_e)\ge pc_pe\ge0.                               \tag{J}
\]
For input \(W^h x^b\), select indices \(e=(e_1,\ldots,e_s)\), \(E=\sum e_j\).
Cartier retains exactly the terms for which \(b+Se\equiv0\pmod p\).
Homogeneity and \(p\nmid d\) imply \(h+E\equiv0\pmod p\). The output is
\[
 k=\frac{h+E}{p},\quad b'=\frac{b+Se}{p},\quad
 \text{coefficient}=(-p)^k\prod_j\ell_{e_j}\widehat a_j^{e_j}.       \tag{K}
\]
The remaining factor is
\(\pi^{h+E-k}=\pi^{(p-1)k}=(-p)^k\) when Cartier divides the
\(w\)-exponent. This derives C:142-145 and 223-238 from (2.14).

The exact coefficient valuation and its bound are
\[
 k+\sum_jv_p(\ell_{e_j})=\frac{h+\sum_jw_{e_j}}p
             \ge h/p+c_pE.                                      \tag{L}
\]
The budget requires this valuation to be strictly less than \(N\)
(C:148-159,234). Therefore
\[
 E<\beta(p)(pN-h),\qquad
 \boxed{k<\beta(p)N-(\beta(p)-1)h/p.}                             \tag{M}
\]
Cancellation of tuples can remove output keys but cannot increase this bound.
Since \(h>0\), \(p\ge5\), and \(\beta(p)\le20/11\),
\[
 k<20N/11,\qquad
 k\le K_h:=\left\lceil\beta(p)N-(\beta(p)-1)h/p\right\rceil-1
       \le\lceil20N/11\rceil-1.                                 \tag{N}
\]
For an empty image the bound is vacuous. Output positive exponents have
degree \(dk\); ordinary numerators have degree \(dk-n\). These bounds are
uniform in \(p\) at fixed \(n,d,N\); (H) gives the input-pole range.

C:12-18 uses the weaker input-independent index bound
\[
 L=\lceil N/c_p\rceil-1=\lceil p\beta(p)N\rceil-1.                  \tag{O}
\]
Equation (L) also bounds the sum \(E\le L\) at every accepted prefix,
justifying the accumulated-index and coordinate bound at C:193-198.
The splitting indices are \(O(pN)\), whereas the poles after Cartier are
\(O(N)\). Confusing these two sizes incorrectly attributes a factor \(p\)
to reduction.

For \(K=\lceil20N/11\rceil-1\), one image has at most
\[
 B_K=\sum_{\substack{1\le k\le K\\dk\ge n}}\binom{dk-1}{n-1}        \tag{P}
\]
distinct keys: count positive compositions of each \(dk\).
Thus dictionaries passed to reduction have size bounded independently of \(p\),
even if many tuples contribute to one key.

There is a strengthening relevant to the lift: (L) is strictly positive,
and its left side is an integer. Every raw coefficient is therefore
divisible by \(p\). Integral reduction preserves this, so \(A/p\) is integral
in the Dwork basis.

## 3. Why one finite Jacobian table suffices

Let \(R_N=\mathbb Z/p^N\mathbb Z\), \(J=(f_1,\ldots,f_n)\).
Euler's relation \(df=\sum x_i f_i\), with \(p\nmid d\), identifies common
projective zeros of the partials over \(\mathbb F_p\) with singular points.
When there are none, the partials are a homogeneous regular sequence of
degree \(d-1\). The quotient Hilbert series is
\[
 (1+t+\cdots+t^{d-2})^n,                                        \tag{Q}
\]
with socle degree \(\delta=n(d-2)\). The lifted graded quotients are free
over \(\mathbb Z_p\) and \(R_N\). One justification is degreewise lifting
of the exact mod-\(p\) Koszul complex: in each degree its modules are finite
free, its mod-\(p\) acyclicity lifts, and its cokernel is free.
Unit pivots therefore lift a complement without \(p\)-division. This is
the ordinary-partial specialization of the freeness input in
[SV Lemma 2.1][SV-basis].

At \(C_0=\delta+1\), every monomial is in \(J\).
For any higher-degree monomial, select a divisor \(x^\eta\) of degree \(C_0\)
and multiply its stored identity (C) by the remaining factor \(x^\xi\).
This generates decompositions in all higher degrees from one table.
R:186-193 selects that divisor greedily.

Below \(C_0\), only numerator degrees \(t=kd-n\ge0\), \(t\equiv-n\pmod d\),
can occur. R:156-157 tabulates exactly those degrees. Their remainders are
basis monomials. At high degree the remainder is zero; at low degree
\(\eta=u,\xi=0\). Therefore R:194-196 can output a remainder directly;
there is no missing multiplication by \(x^\xi\).
Each other child drops the pole by one, so descending layers R:178-205
terminate. At a would-be pole below \(\lceil n/d\rceil\), derivatives
vanish or cofactor degree precludes a term.

The dimension, summing coefficients of (Q) in this congruence class, is
\[
 r=\frac{(d-1)^n+(-1)^n(d-1)}d.                                 \tag{R}
\]
A roots-of-unity filter gives this directly: at a nontrivial \(d\)-th root,
\(1+t+\cdots+t^{d-2}=-t^{-1}\). This derives R:159-160 and (H).

Conversely, surjectivity at \(C_0\) proves smoothness over the algebraic
closure. At a putative common projective zero \(z\) of the partials, take
\(z_i\ne0\). A representation of \(x_i^{C_0}\) in \(J\) would evaluate
to a nonzero value equal to zero. Thus R:153-155 is not merely checking
for rational singularities.

### The abstract reduction problem

The thesis, main.tex:913-1104, and [arXiv §3.4][abstract-game] use a fixed
coefficient block and an outer exponent: a nonnegative direction
\(v\), \(|v|=d\), moves \(u\) to \(u-kv\); \(R_{u,v}\) transforms block
coefficients. Setup costs \(T\), each multiplication costs \(M\), and
merging duplicates costs \(D\). The arXiv text calls the setup constant
\(S\), rather than \(T\).

For a fixed table index \(\eta\), (D) has coefficients
\[
 -c_{\eta,i,g}\xi_i-c_{\eta,i,g}g_i,                              \tag{S}
\]
affine-linear in the outer exponent \(\xi\).
Stored cofactors supply both constant and linear pieces on demand.
They need not be separately assembled into square matrices for each direction.
CR:59-100 performs that assembly for a fixed block and direction;
R:193-203 instead evaluates a local transition. Dictionary additions at
R:175,196,203 immediately combine like terms. Jacobian-table construction
is a shared setup cost, not a new \(T\) per chunk.

There is, however, no literal identification of the raw monomial reducer
with that game. The selected \(\eta\) has norm \(C_0\), not \(d\), and one
monomial can branch. Each child's displacement is
\[
 \Delta=\eta-g+e_i,\qquad \sum_i\Delta_i=d,                        \tag{T}
\]
but coordinates of \(\Delta\) can be negative. It need not belong to
\(V_{\rm all}(u,1)\), much less a chosen \(V_S(u,1)\).
For a concrete smooth supported example, take
\(f=x^3+y^3+z^3+x^2y\) over \(\mathbb F_5\). Its partials have no
common projective zero: \(z=0\), and \(3x^2+2xy=x^2+3y^2=0\)
forces \(x=y=0\). Here \(C_0=4\), and the identity
\[
 x^4=\frac{9x^2-6xy}{31}f_x+\frac{4x^2}{31}f_y
\]
is valid over \(\mathbb Z_5\). Multiplying by \(x^2\) and applying (B)
produces
\[
 W^3x^6\equiv-\frac{W^2}{31}(36x^3-18x^2y).
\]
The child \(x^2y\) has displacement \((4,-1,0)\), whose coordinate sum
is 3 but whose \(\ell_1\)-norm is 5. A finite-field transcription of
the PR's table order also gives children \(4x^3+3x^2y\) modulo 5,
so this issue occurs in its actual table policy, not only in an arbitrary
decomposition. The original positive input exponent is \((7,1,1)\).
For diagonal \(f\), a nonzero step has direction \(d e_i\), which fits
directly. In general, literal game equivalence requires a fixed-block
regrouping and suitable logarithmic cofactors/boundary conditions.
The common reduction mechanism alone does not identify state spaces or
justify assigning the score \(Mk+T\) to every raw dictionary operation.

There is also a definition inconsistency in the cited thesis.
Main.tex:969-972 defines \(V_{\rm full}\) by
\(u_i\ne0\Rightarrow v_i\ne0\), whereas 985-991 defines \(V_S\) by
\(u_i=0\Rightarrow v_i=0\), yet asserts their equality for all indices.
For \(u=(d,1)\), \(v=(d,0)\), \(k=1\), \(v\) is in \(V_{\rm all}\),
preserves zero coordinates, and fails full support. These are different
predicates. In fact the stated zero-preservation condition is already
implied by \(u-kv\ge0\), so that version of \(V_S\) adds no restriction.
The actual predicate
[ZetaFrontend.jl:46-56][directions] uses the nonzero-support rule on \(S\).
The thesis also uses an index range involving the target norm in place of
the variable index range. These definitions need correction before a precise
cross-document \(V_S\) dictionary; the PR reducer does not depend on them.

## 4. Cost in p and the square-root question

Fix \(n,d,N\), and hence rank \(r\). Let
\(\rho=\operatorname{rank}(S\bmod p)\), \(t=s-\rho\).
The operation count differs from the bit cost of arithmetic on
\(O(N\log p)\)-bit integers. A tilde hides logarithms in \(p\).

| Stage | Code | Dependence on \(p\) at fixed \(n,d,N,s\) |
| --- | --- | --- |
| Lifts and Jacobian/basis setup | R:137-164; C:34 | No positive power; fixed-size linear algebra and modular powering |
| Splitting tables | C:21-87 | \(\widetilde O(p)\); \(L+1=\Theta(pN)\) |
| Residue-system elimination | C:92-125 | No positive power; fixed dimensions |
| Residue and tuple enumeration | C:209-239 | \(\widetilde O(p^t)\) for the fixed basis |
| Dictionary accumulation | C:142-145 | Same enumeration cost, with at most \(B_K\) distinct keys |
| Cohomological reduction | R:167-208 | No positive power; bounded degree and number of keys |
| Integer charpoly and lift | I:3-37 | Polynomial in \(\log p\) at fixed \(r,N\) |

The compatible residue system \(Sr=-b\) has \(p^t\) solutions; C:207
skips incompatible right-hand sides. Each allowed integer index in a fixed
residue bucket lies between 0 and \(L\), so each bucket has \(O(N)\) terms.
An explicit enumeration bound is
\[
 O\bigl(r\,s\,p^t(\lceil\beta(p)N\rceil+1)^s\bigr),                \tag{U}
\]
before modular-arithmetic bit costs. Pruning can reduce the tuple work, but
the outer loop visits every free-residue parameter choice.

The pure-power columns give the submatrix \(dI_n\), so \(\rho=n\).
Pure powers alone give \(s=n\), \(t=0\); one mixed term gives \(s=n+1\),
\(t=1\). The table loop is linear in either case. In fact
\(\ell_e=1/e!\), \(w_e=e<pN\) for every \(0\le e<p\), so the tables
explicitly store at least \(p\) entries per support term. They also allocate
\(sp\) buckets and a history buffer of length \(p\).
The current implementation therefore cannot be sublinear in \(p\),
even on a diagonal equation. With unrestricted support the corresponding
cost would be \(\widetilde O(p^{\max(1,s-\rho)})\), not necessarily linear.

The exponent \(p^{N-1}\) in the lift at R:138/C:34 does not mean
\(p^{N-1}\) multiplications: binary modular powering uses an exponent of
\((N-1)\log_2p+O(1)\) bits.

[SV Theorem A and (4.7)][SV] state the bit bound
\[
 \widetilde O\!\left(s^{\lceil m/2\rceil}+pN^3\log q+
      p^{s-\rho}(6N+m)^s\,v^4N\log^2q\right),                     \tag{V}
\]
with \(m\) the paper's torus dimension and \(v\) its normalized volume.
At fixed dimension and support size, Theorem B states
\[
 \widetilde O\!\left(p^{\max(1,s-\rho)}v^{s+5}\log^{s+3}q\right).
                                                                    \tag{W}
\]
For homogeneous support, the augmented columns \((1,\nu_j)\) have the same
rank as \(S\), since their first row is \(d^{-1}\) times the sum of the
remaining rows modulo \(p\). These are consistent with the PR analysis,
and are not square-root claims.

The public routine increases \(N\) by two when lifting is ambiguous (I:58-66).
This introduces no hidden positive power of \(p\) at fixed \(n,d\).
A conservative sufficient precision is
\[
 N>r(1+(n-2)/2)+\log_5(2r).                                      \tag{X}
\]
For \(N\ge r\), the coefficient modulus is at least \(p^{N-k}\), while
the trace interval has width \(2r p^{(n-2)k/2}\). Bound (X) makes the
interval shorter than the spacing \(k p^{N-k}\) for every \(k\le r\).
There are consequently boundedly many precision attempts, uniformly for
\(p\ge5\). Column valuations usually improve this conservative bound.

**Square-root assessment: a research possibility, not an established bound.**
A Bostan–Gaudry–Schost-style product method can plausibly accelerate an
affine-linear recurrence of length \(\Theta(p)\). Here cohomological
reduction has only \(O(N)\) steps, so accelerating it does not improve
the exponent of \(p\). A sublinear algorithm would have to avoid table
materialization and, for \(t=1\), replace the explicit \(p\)-term sum by
a structured recurrence or block summation. Faster isolated coefficients
alone would leave that scan intact.

Moreover the current coefficient recurrence has lag \(p\), and its unit
normalization changes at multiples of \(p\). A fixed-dimensional recurrence
suitable for acceleration must be derived; using an order-\(p\) companion
matrix would not establish a square-root algorithm. Such research is
consistent with [SV's discussion of acceleration, introduction p.4][SV-acceleration].
No square-root bound for this implementation or a replacement is proved here.

## 5. Audit of the requested correctness concerns

### Integrality with unit pivots

The freeness argument after (Q) justifies R:71: the mod-\(p\) rank gives
unit minors, and there is no hidden \(p\)-torsion in a lifted graded quotient.
R:91-93 rejects residual rows with nonzero polynomial columns.
Negative-index columns track cofactors, so (C) is an identity modulo \(p^N\),
not just a dimension count. Reduction then only adds, multiplies, and
differentiates in \(R_N\), with no division by pole order.

### The (1.4) cutoff and coefficient working precision

Equations (I)-(O) prove that an omitted splitting index cannot contribute
modulo \(p^N\), including after reduction. Nonnegative weights justify both
prefix pruning and dropping factors before multiplication.

For the coefficient computation let \(F_e=v_p(e!)\),
\(C_e=p^{F_e}\ell_e\), and \(e=p^a u\) with \(u\) a unit. The exact
scaled recurrence is
\[
 C_e=u^{-1}\bigl(C_{e-1}+p^{F_{e-1}-F_{e-p}}C_{e-p}\bigr),         \tag{Y}
\]
omitting the second summand for \(e<p\). Its exponent is nonnegative.
This is C:43-58, including the buffer indexing.
If the coefficient survives, then \(w_e<pN\) implies
\[
 v_p(C_e)<N+F_e-e/p
          \le N+v_p(\lfloor L/p\rfloor!).                         \tag{Z}
\]
The work modulus is \(p^{2N+g}\),
\(g=v_p(\lfloor L/p\rfloor!)\) (C:27-33).
After extracting the \(p\)-power, at least \(N\) correct digits of every
retained unit remain. A true \(C_e\) equal to zero modulo the work modulus
cannot survive (Z). Thus zero detection and unit extraction are both safe.
Computing an unscaled rational recurrence at only \(N\) digits would not
justify the same claim.

### Teichmüller lifts

For any integer unit \(a=\omega(a)\langle a\rangle\), where
\(\omega(a)^{p-1}=1\), \(\langle a\rangle\in1+p\mathbb Z_p\),
\[
 a^{p^{N-1}}\equiv\omega(a)\pmod {p^N}.
\]
Indeed \(p^{N-1}\equiv1\pmod {p-1}\), while that power of
\(\langle a\rangle\) is 1 modulo \(p^N\).
R:138 and C:34 therefore use the same valid lifts. The API is prime-field
only; no extension-field coefficient Frobenius or norm is missing.
Ordinary integer lifts in that kernel would need a different comparison.

### Characteristic polynomial and Weil-bound reconstruction

Let \(\widetilde A\) be the integer matrix formed from residues and
\(A\equiv\widetilde A\pmod {p^N}\) the true operator in a lifted basis.
Sort the column lower bounds
\[
 v_j=\min\bigl(N,\min_i v_p(\widetilde A_{ij})\bigr).
\]
An all-zero residue column has bound \(N\) (I:5-7). Expanding principal
minors columnwise proves
\[
 [T^k]\det(1-T\widetilde A)\equiv[T^k]\det(1-TA)
       \pmod {p^{N+v_1+\cdots+v_{k-1}}}.                          \tag{AA}
\]
With \(j\ge1\) changed columns, the valuation is at least
\(jN+\sum_{\ell=1}^{k-j}v_\ell\), which is no smaller than the exponent
in (AA), since every \(v_\ell\le N\). Thus multiple changed columns do
not invalidate the bound.

Write \(P(T)=1+\sum c_kT^k\). By the projective comparison,
\([T^k]\det(1-TA)=p^k c_k\).
I:14,33 maintains \(\texttt{digits}=N+\sum_{\ell<k}v_\ell\);
I:18 checks divisibility, and I:20-21 yields the valid residue modulus
\(m_k=p^{\max(0,\texttt{digits}-k)}\).
If digits is less than \(k\), the modulus is one; integer division of
the representative then asserts no nonexistent precision.
No filtration-adapted ordering is required.

For the power sums \(s_k\), Newton's identity is
\[
 s_k=-\sum_{i=1}^{k-1}c_i s_{k-i}-k c_k.                           \tag{AB}
\]
Given exact earlier coefficients, the first sum is known.
Purity of weight \(n-2\) gives \(|s_k|\le r p^{(n-2)k/2}\).
For integer \(s_k\) the exact bound is
\(\lfloor r p^{(n-2)k/2}\rfloor
=\operatorname{isqrt}(r^2p^{(n-2)k})\), as in I:23.
Substitute \(c_k=\texttt{residue}+m_kj\) into (AB); the allowed
integers \(j\) are exactly the interval at I:25-26.
The true coefficient lies there. A singleton proves the reconstructed
coefficient by induction; ambiguity calls for a new computation.
I:35's determinant-square test is consistent with primitive self-duality,
but is a check, not a replacement for this proof.

Thus the comment “\(A=pF\)” should be understood through the
characteristic-polynomial normalization and a specified comparison, not as
entrywise equality with another routine's matrix.

### Coordinate-section smoothness

Write \(f=\sum a_i x_i^d+\lambda x^m\), with the mixed summand optional.
Set some variables to zero. If the mixed term vanishes, the section is
diagonal: the partials \(d a_i x_i^{d-1}\) cannot all vanish at a
projective point. If it survives, every omitted variable has exponent zero
in it. Its normal derivatives, and those of omitted pure powers, vanish
on the section. A singular point of the section would then be singular
on the full hypersurface. This proves D:23-27 over the algebraic closure.
Empty coordinate sections cause no exception.

The nonzero pure powers make the Newton polytope the full degree-\(d\)
simplex. Its faces are coordinate faces. These section arguments supply
the general-position/nondegeneracy hypothesis for the projective comparison.

### At most one mixed monomial

This restriction has two uses: it validates that section argument and forces
\(s-\rho\le1\). It is not required for identity (B), the ordinary-Jacobian
socle bound, or the integral reducer considered separately.
With several mixed terms, a section can lose some and retain others; a
normal derivative of a lost term can rescue a section-singular point on the
full hypersurface. Full smoothness alone then need not imply the required
general position. Generalizing the API would require checking that condition
and accepting a possibly higher power of \(p\) in enumeration.

**Verdict within scope.** The requested reduction, truncation, lift, and
support claims have the preceding mathematical justifications. No
counterexample to the guarded algorithm was found. The short PR documentation
omits the factorial dictionary, the uniform output-pole bound, the graded
flatness argument, and the guard-precision proof (Z). These are exposition
obligations, not evidence of a failed calculation. Runtime behavior,
Oscar/Julia integration, empirical scaling, and correctness outside the
public input domain are not claimed to have been tested.

## Checks, limitations, and publication record

This is a self-audited informal report, not independent verification or a
machine-checked proof. It imports the projective Dwork theorem and the Weil
bound; it does not reprove p-adic cohomology.

Independent exact arithmetic used Python 3.9.6 standard-library rational
numbers and the closed formula
\[
 \ell_e=\sum_{j=0}^{\lfloor e/p\rfloor}
            \frac1{(e-pj)!\,j!\,p^j}.                            \tag{AC}
\]
For \(p=5,7,11,13\), \(N=1,\ldots,5\), all 735 coefficients for
\(0\le e\le L\) agreed with the scaled recurrence modulo the work modulus.
All 525 retained units passed the precision and weight checks.
The Teichmüller formula was checked for four integer-unit representatives
at each parameter pair. These finite computations are indexing/normalization
sanity evidence, not the proofs of the general bounds.

For pure powers plus the product monomial, input \(b=\mathbf1\), \(N=3\):

| \(p\) | \((n,d)=(3,3)\): retained tuples | \((n,d)=(4,4)\): retained tuples | Output poles in both cases |
| --- | --- | --- | --- |
| 5 | 8 | 7 | 1,2 |
| 7 | 14 | 13 | 1,2 |
| 11 | 20 | 17 | 1,2 |

These are tuple counts before cancellation, not benchmarks or zeta tests.
All accepted tuples satisfied the pole bound, the total-index cutoff, and
positive integral coefficient valuation.

Reproducibility baseline: owned branch gc-mathcats.mathcat-7a655e6bdaef,
initial checkout 4ff65f3a4c80dd3076e9e619ba828ea96d1221aa.
Target files were read using git show at the review commit.
No checkout, package loading, package tests, or E2E changes were needed.
Pre-existing .gitignore changes and .task-dr-eevm-* directories were preserved.

The tracked artifact is this owned-worktree file. The requested additional
destination is
/Users/Shared/Developer/ai/city/math-code-projects/DeRham.jl/.gc/maintenance/pr97-sparse-dwork-20261123/math-justification.md.
Base worker instructions prohibit writes to that canonical checkout.
The conflict was escalated by durable mail mc-wisp-trqzo56; its resolution
and the final commit/branch are recorded on the source bead.
No canonical source files were changed.
Follow-ups dr-543md, dr-h5a8g, and dr-q37iq were filed with discovered-from
links to dr-flnzq. Self-audit passed: all five assignment questions and all
eight metadata fields are present, references and math delimiters resolve,
and proof, imported theorems, finite computation, and unperformed runtime
checks are distinguished. The canonical-copy obligation remains explicit.

[R]: https://github.com/UCSD-computational-number-theory/DeRham.jl/blob/21f0ba6309e0649b4b710cb22943815b9d4461d6/src/finalreduction/SparseDworkReduction.jl
[C]: https://github.com/UCSD-computational-number-theory/DeRham.jl/blob/21f0ba6309e0649b4b710cb22943815b9d4461d6/src/frobenius/SparseCartier.jl
[I]: https://github.com/UCSD-computational-number-theory/DeRham.jl/blob/21f0ba6309e0649b4b710cb22943815b9d4461d6/src/invariants/SparseDwork.jl
[D]: https://github.com/UCSD-computational-number-theory/DeRham.jl/blob/21f0ba6309e0649b4b710cb22943815b9d4461d6/docs/src/experimental/sparse_dwork.md
[CR]: https://github.com/UCSD-computational-number-theory/DeRham.jl/blob/21f0ba6309e0649b4b710cb22943815b9d4461d6/src/reduction/ComputeRuv.jl
[SV]: https://jvoight.github.io/articles/sparse-dwork-073025.pdf
[SV-projective]: https://jvoight.github.io/articles/sparse-dwork-073025.pdf#page=24
[SV-basis]: https://jvoight.github.io/articles/sparse-dwork-073025.pdf#page=10
[SV-acceleration]: https://jvoight.github.io/articles/sparse-dwork-073025.pdf#page=4
[abstract-game]: https://arxiv.org/html/2602.24155v1#S3.SS4
[prescale]: https://github.com/UCSD-computational-number-theory/DeRham.jl/blob/21f0ba6309e0649b4b710cb22943815b9d4461d6/src/frontends/ZetaFrontend.jl#L18-L24
[scalar]: https://github.com/UCSD-computational-number-theory/DeRham.jl/blob/21f0ba6309e0649b4b710cb22943815b9d4461d6/src/quotients/CohomologyQuotient.jl#L47-L51
[final-scale]: https://github.com/UCSD-computational-number-theory/DeRham.jl/blob/21f0ba6309e0649b4b710cb22943815b9d4461d6/src/finalreduction/FinalReduction.jl#L149-L175
[directions]: https://github.com/UCSD-computational-number-theory/DeRham.jl/blob/21f0ba6309e0649b4b710cb22943815b9d4461d6/src/frontends/ZetaFrontend.jl#L46-L56
