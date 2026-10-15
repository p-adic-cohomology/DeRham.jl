# Getting started

DeRham.jl requires Julia 1.11 or newer. The documentation workflow tests Julia 1.12.

## Install

In a Julia session, install DeRham.jl with:

```julia
using Pkg
Pkg.add(url="https://github.com/p-adic-cohomology/DeRham.jl")
```

Installing DeRham.jl downloads the CUDA stack, but a CUDA GPU is not needed for
installation or ordinary computations. Only computations with `use_gpu=true`
require a CUDA GPU.

Load DeRham.jl together with Oscar:

```julia
using DeRham, Oscar
```

## First computation

TODO: T3a will add a doctested Fermat example after the exports land in T2.
