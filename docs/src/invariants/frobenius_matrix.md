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

When we compute the Frobenius matrix, we always compute each vector with a precision relative to the p-adic valuations guaranteed to us by strong Mazur divisibility. Thus, if we choose relative precision $r = 2$ for our curve of genus 3, the first three vectors will be computed modulo $p^3$ and the second three will be computed modulo $p^2$. If what you cared about was the *absolute* precision (i.e. the whole matrix modulo $p^2$), you can reduce the answer modulo $p^2$.

We also provide the option to choose a relative precision for each basis vector.
`DeRham.frobenius_matrix_with_precision` computes the returned matrix with
exactly the requested relative precision for each Hodge-filtered piece, using
the certificate proved in dr-eevm.1, “Certified series precision: proof of the
fixed algorithm.” It does not automatically increase the requested precision.

```@docs
DeRham.frobenius_matrix_with_precision
```
