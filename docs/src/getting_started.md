# Getting started

DeRham.jl is research software. Use Julia 1.10 or newer.

## Install

Clone the repository and instantiate its project dependencies:

```sh
git clone https://github.com/UCSD-computational-number-theory/DeRham.jl.git
cd DeRham.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The project includes Oscar and the other required Julia dependencies. The
initial instantiation downloads and compiles them and may take a while.

## First calculation

Start Julia from the repository directory and load DeRham.jl together with
Oscar:

```sh
julia --project=.
```

At the Julia prompt, compute the zeta function coefficients and Newton polygon
of a smooth cubic curve over `𝔽₇`:

```julia
using DeRham, Oscar

p = 7
R, (x, y, z) = GF(p)[:x, :y, :z]
f = y^2 * z - x^3 - x * z^2 - z^3

zeta_coeffs = DeRham.zeta_coefficients(f)
DeRham.L_polynomial(zeta_coeffs)
DeRham.newton_polygon(f)
```

For computations that use a GPU, install a compatible CUDA driver first; see
the [CUDA.jl installation guide](https://cuda.juliagpu.org/stable/installation/overview/).
