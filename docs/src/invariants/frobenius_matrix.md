# Frobenius matrix

## Projective Hypersurfaces

Let $X$ be a projective hypersurface over $\mathbb{F}_p$ of dimension $n$. We compute the Frobenius matrix on crystalline cohomology

```math
F --insert definition of crystalline cohomology here--
```

Because this matrix lives over $\mathbb{Z}_p$, we can only compute an approximation modulo $p^N$ for some $N$.

Our algorithm actually gives something much nicer; we give a basis for $H^2_crys$ which is ordered by the hodge filtration. Furthermore,
for a basis element $e_{i,m} \in Fil^m(H^2_crys)$, our we have the strong Mazur divisibility property:

```math
F(e_{i,m}) \in p^{m} H^2_crys.
```

This is very strong and uses many special properties, including the fact that hypersurfaces are liftable and Mazur-Ogus, the computation of de Rham cohomology of hypersurfaces due to Griffiths-Dwork, various comparison theorems in characteristic zero, and a delicate argument in Abbott-Kedlaya-Roe which compares the integral lattice in the Griffiths-Dwork basis to the lattice of crystalline cohomology in rigid cohomology.

However, what we get is that our Frobenius matrix has the strong Mazur divisibility. For example, for a genus 3 curve

```math
Insert matrix of a genus three curve here. Should be three columns of * ... * and then three colums of p* ... p*, the *'s denoting any number
```

When we compute the Frobenius matrix, we always compute each vector with a precision relative to the p-adic valuations guaranteed to us by strong Mazur divisibility. Thus, if we choose relative precision $r = 2$ for our curve of genus 3, the first three vectors will be computed modulo $p^2$ and the second three will be computed modulo $p^3$. If what you cared about was the *absolute* precision (i.e. the whole matrix modulo $p^2$), you can reduce the answer modulo $p^2$.

We also provide the option to choose a relative precision for each basis vector.

```@docs
DeRham.frobenius_matrix_with_precision
```
