# Frobenius matrix

## Projective Hypersurfaces

Let $X$ be a projective hypersurface over $\mathbb{F}_p$ of dimension $n$. We compute the Frobenius matrix on crystalline cohomology:

```math
F \colon H^n_{\mathrm{crys}}(X/\mathbb{Z}_p) \to H^n_{\mathrm{crys}}(X/\mathbb{Z}_p)
```

The absolute Frobenius of $X$ induces $F$ by functoriality. DeRham computes the Frobenius matrix on the primitive part of middle cohomology; for a K3 surface, this gives a $21\times 21$ matrix. The module $H^n_{\mathrm{crys}}$ is a finitely generated $\mathbb{Z}_p$-module, free for a smooth hypersurface with rank equal to the middle Betti number; the matrix is the matrix of $F$ in a $\mathbb{Z}_p$-basis.

Because this matrix lives over $\mathbb{Z}_p$, we can only compute an approximation modulo $p^N$ for some $N$.

Our algorithm actually gives something much nicer; we give a basis for $H^n_{\mathrm{crys}}$ which is ordered by the Hodge filtration. Furthermore,
for a basis element $e_{i,m} \in \mathrm{Fil}^m H^n_{\mathrm{crys}}$, we have the strong Mazur divisibility property:

```math
F(e_{i,m}) \in p^{m} H^n_{\mathrm{crys}}.
```

This is very strong and uses many special properties, including the fact that hypersurfaces are liftable and Mazur–Ogus, the computation of de Rham cohomology of hypersurfaces due to Griffiths–Dwork, various comparison theorems in characteristic zero, and a delicate argument in Abbott–Kedlaya–Roe which compares the integral lattice in the Griffiths–Dwork basis to the lattice of crystalline cohomology in rigid cohomology.

However, what we get is that our Frobenius matrix has the strong Mazur divisibility. For example, for a genus 3 curve, the Frobenius matrix on $H^1_{\mathrm{crys}}$ has the form

```math
\begin{pmatrix}
p\ast & p\ast & p\ast & \ast & \ast & \ast \\
p\ast & p\ast & p\ast & \ast & \ast & \ast \\
p\ast & p\ast & p\ast & \ast & \ast & \ast \\
p\ast & p\ast & p\ast & \ast & \ast & \ast \\
p\ast & p\ast & p\ast & \ast & \ast & \ast \\
p\ast & p\ast & p\ast & \ast & \ast & \ast
\end{pmatrix}
```
where each $\ast$ denotes an element of $\mathbb{Z}_p$.

The result is an ordinary Oscar matrix over $\mathbb{Z}/p^m$, where $m$ is at
least the largest absolute precision needed by its columns. A requested
relative precision applies to each Hodge-filtered piece, and every column in
that piece is guaranteed at least that relative precision. The common matrix
ring can represent more digits than are known for some columns; those higher
digits are not meaningful. This contract applies to both all-zero and nonzero
precision requests and is determined by the request and chosen basis.
By default, the matrix carries no per-column precision metadata. Set
`return_column_moduli=true` to also return `(FM, column_moduli)`, where each
column's value is the modulus to which it is known. Pass that vector to
matrix-level invariant methods when columns have different known precision;
without it, they assume every column is known to the full base-ring modulus.

We also provide the option to choose a relative precision for each basis vector.
`DeRham.frobenius_matrix_with_precision` computes the returned matrix with
at least the requested relative precision for each Hodge-filtered piece. An
existing automatic table floor can provide more. The certified truncation
argument is summarized below. The full proof, “Certified series precision:
proof of the fixed algorithm,” is not included in this repository; see
`dr-eevm.1` for the proof report.

### [Certified series truncation](@id sec-certified-series-truncation)

Here $n$ is DeRham's ambient projective dimension, $p$ is prime, $m$ is a pole
order from $1$ to $n$, and $r_m$ is its requested relative precision. Assume a
smooth hypersurface and smooth lift with the AKR monomial basis and Frobenius
normalization. The repaired AKR bound $A(t)$ is an upper bound for the H-lattice
defect and is pointwise at most
$f_0(t) = \sum_{j=1}^{n} \lfloor \log_p(\max(1,t-j)) \rfloor$. For the
coordinate lattice, set $c = 0$ when $p \ge n$ and
$c = v_p((n-1)!)$ when $p < n$; the lattice comparison gives the upper bound
$B_V(t) = A(t) + c$.

For a proposed tail cutoff $T$, the implementation checks, using integer
arithmetic, that $m+T \ge 2n+1$, that
$E = m-1+T-r_m-c \ge 0$, and that $(p(m+T)-1)^n \le p^E$. The AKR bound and
the floor-log estimate imply $B_V(p(m+i)) \le c +
n\log_p(p(m+i)-1)$. Define
$h(u) = m-1+u-r_m-c-n\log_p(p(m+u)-1)$. The integer checks give
$h(T) \ge 0$; for every $u \ge T$ with $m+u \ge 2n+1$,
$h'(u) = 1 - n/((m+u-1/p)\ln p) > 0$. Therefore
$m-1+i-B_V(p(m+i)) \ge r_m$ for every integer $i \ge T$.

The selector checks the same inequality exactly for every $i < T$ using the
repaired bound. This gives a candidate retained-term count: one more than the
largest failing index, or zero if no index fails. The public count $N_m$ is
raised to the Hodge-bound minimum $n-m+1$, so it can be positive even when
there is no failing index. Every omitted term $i \ge N_m$ then satisfies the
inequality, which is Costa's truncation condition. Such a $T$ always exists;
the proof gives
$T \le \max\{16n^2, 2(r_m+c+1), m\}-m$. Raising $N_m$ to an existing table
floor preserves the certificate. The AKR lattice correction requires the
specified integral basis hypotheses; an arbitrary rational basis is not
covered. For zero requests with $p < n$, the API raises `ArgumentError`; for
$p \ge n$, an all-zero request returns the Hodge-bound zero matrix.

```@docs
DeRham.frobenius_matrix_with_precision
```
