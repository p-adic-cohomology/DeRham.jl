# DeRham.jl

[![CI](https://github.com/UCSD-computational-number-theory/DeRham.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/UCSD-computational-number-theory/DeRham.jl/actions/workflows/CI.yml)

This package aims to implement cutting-edge algorithms that calculate

- point counts of equations mod $p$ (algebraic varieties over $\mathbb{F}_p$)
- the cohomology of algebraic varieties
- the extra structure on this cohomology

Long term goals include having fast versions of all the various forms of Kedlaya's algorithm, Harvey's algorithm, Harvey's Trace Formula, Moonen's algorithm, etc.

**This package is work-in-progress research software.** It may have rough edges. PRs and contributions are welcome.

## Current functionality

We currently implement a variant of Kedlaya's algorithm that computes the zeta functions (and consequently the newton polygons) of smooth projective hypersurfaces of degree $d$ over $\mathbb{F}_p$ (for $p\nmid d$) using a method known as controlled reduction, see Proposition 1.15 [here](https://edgarcosta.org/assets/articles/EdgarCosta-PhDthesis.pdf).

## Getting Started

See [Getting started](https://p-adic-cohomology.github.io/DeRham.jl/dev/getting_started/) in the documentation for installation instructions.

## Sample code

<!-- examples -->
#### Zeta Functions and Newton polygons

Execute the following lines in Julia REPL:

```julia
p = 7
R, (x,y,z) = GF(p)[:x,:y,:z]
f = y^2*z - x^3 - x*z^2 - z^3
zeta_coeffs = DeRham.zeta_coefficients(f)
DeRham.L_polynomial(zeta_coeffs)
DeRham.newton_polygon(f)

g = DeRham.random_hypersurface(3,4,p)
DeRham.zeta_coefficients(g,givefrobmat=true)

h = DeRham.fermat_hypersurface(4,4,p)
DeRham.newton_polygon(h)
```

#### Nondegeneracy Conditions (advanced)

We support using the nondegeneracy condition in [Costa's Thesis](https://edgarcosta.org/assets/articles/EdgarCosta-PhDthesis.pdf) called "S-smoothness", where "S" is a subset of `{0, ..., n-1}` (where `n` is the number of variables in the polynomial ring R).

By default, S is the full set `{0, ..., n-1}`, which is equivalent to the smoothness of the hypersufaces. Most users will not have
to worry about this.

There are two situations (currently) which require the use of S. First, if the degree `d` of the hypersurface is less than the number of varibles `n`. In this case, S can have size at most 3.

```julia
p = 7
R, (x,y,z,w) = GF(p)[:x,:y,:z,:w]
f = x^3 + y^3 + z^3 + w^3

DeRham.zeta_coefficients(f,S=[0,1,2])
```

Secondly, if one uses that `:varbyvar` reduction policy, then S must be `{n-1}`

```julia
p = 11
R, (x,y,z,w) = GF(p)[:x,:y,:z,:w]
f = x^4 + y^4 + z^4 + w^4 + x*y*z*w

DeRham.zeta_coefficients(f,S=[3], algorithm=:varbyvar)
```

#### Parallel Programming

For cubic threefolds and cubic fourfolds, if you have the CUDA driver installed, you can run the following code to use your Nvidia GPU:

```julia
p = 7
R, (x,y,z,w,v) = GF(p)[:x,:y,:z,:w,:v]
f = x^4 + y^4 + z^4 + w^4 + v^4

DeRham.zeta_coefficients(f, S=[4], algorithm=:varbyvar, use_gpu=true)
```

If you start Julia with `julia --threads n` or `julia --threads auto`, you can also use Julia's mutlithreading. This is most effective for quartic and quintic surfaces:

```julia
p = 7
R, (x,y,z,w) = GF(p)[:x,:y,:z,:w]
f = x^4 + y^4 + z^4 + w^4

DeRham.zeta_coefficients(f, use_threads=true)
```

Simultaneous use of `use_gpu=true` and `use_threads=true` is not currently supported.
<!-- /examples -->

## Citing Our Work

You are welcome to use the code in this repository for your own research, but we ask that you please cite our paper [Newton strata realization for hypersurfaces via explicit p-adic cohomology](https://arxiv.org/abs/2602.24155) (by Ryan Batubara, Jack J Garzella, Yongyuan Huang, and Maximus Mellberg, arxiv:2602.24155 (2026)) if and when you publish your results.
