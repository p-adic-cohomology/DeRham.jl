# DeRham.jl

[![CI](https://github.com/UCSD-computational-number-theory/DeRham.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/UCSD-computational-number-theory/DeRham.jl/actions/workflows/CI.yml)

DeRham.jl is research software for computing point counts, cohomology, and
additional structure on the cohomology of algebraic varieties over finite
fields.

**This package is work in progress.** Contributions and feedback are welcome.

## Current functionality

The package implements a variant of Kedlaya's algorithm using controlled
reduction to compute zeta functions, and consequently Newton polygons, of
smooth projective hypersurfaces of degree `d` over `𝔽ₚ` when `p` does not divide
`d`. See Proposition 1.15 in
[Edgar Costa's thesis](https://edgarcosta.org/assets/articles/EdgarCosta-PhDthesis.pdf).

## Documentation

See the [DeRham.jl documentation](https://ucsd-computational-number-theory.github.io/DeRham.jl/dev/)
for installation instructions and examples.

## Citing

If you use this package in published research, please cite [*Newton strata
realization for hypersurfaces via explicit p-adic cohomology*](https://arxiv.org/abs/2602.24155)
by Ryan Batubara, Jack J. Garzella, Yongyuan Huang, and Maximus Mellberg
(arXiv:2602.24155, 2026).
