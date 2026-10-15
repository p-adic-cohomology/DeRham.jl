# DeRham.jl

[![CI](https://github.com/p-adic-cohomology/DeRham.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/p-adic-cohomology/DeRham.jl/actions/workflows/CI.yml)
[![Documentation (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://p-adic-cohomology.github.io/DeRham.jl/dev/)

DeRham.jl is research software in development for computing point counts and
cohomological invariants of algebraic varieties over finite fields. See the
[development documentation](https://p-adic-cohomology.github.io/DeRham.jl/dev/)
or [stable documentation](https://p-adic-cohomology.github.io/DeRham.jl/stable/)
for an overview, and [Getting started](https://p-adic-cohomology.github.io/DeRham.jl/dev/getting_started/)
for installation instructions. Contributions and feedback are welcome.

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

## Citing

If you use this package in published research, please cite [*Newton strata realization for hypersurfaces via explicit p-adic cohomology*](https://arxiv.org/abs/2602.24155) by Ryan Batubara, Jack J. Garzella, Yongyuan Huang, and Maximus Mellberg (arXiv:2602.24155, 2026).
