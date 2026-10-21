# Getting started

## Installing for the first time

First, install Julia 1.12 and Oscar. If you're new to Julia and Oscar, you can find a tutorial for non-experts [here](https://jjgarzella.github.io/blog/how-to-install-julia/). In short:

- On Windows, install Ubuntu through WSL and use the Ubuntu terminal.
- On Mac, open Terminal and run `xcode-select --install`.
- On Linux, open your usual terminal.

Then install Julia with `curl -fsSL https://install.julialang.org | sh`, restart your terminal if necessary, open a Julia REPL with `julia`, press `]` to enter package mode, and run `add Oscar`. Press backspace to return to the Julia prompt and check that `using Oscar` works.

Then, clone this repo with `git clone https://github.com/UCSD-computational-number-theory/DeRham.jl.git`.
If you'd like to use CUDA, you'll need to install the CUDA driver, see [the instructions in the CUDA.jl docs](https://cuda.juliagpu.org/stable/installation/overview/).

`Revise.jl` is recommended for hot reloading. From a new Julia REPL, run
`using Pkg; Pkg.add("Revise")` (or, using Julia's package prompt, run `] add Revise`) once, and then run `using Revise` upon opening every new REPL.

The first time you install DeRham.jl, it will take a while for Julia to download and compile all of the dependincies. To do this, make sure your Julia REPL is in the DeRham.jl project folder. Then, in the package prompt (i.e. after pressing `]`), run

```
activate .
add https://github.com/UCSD-computational-number-theory/GPUFiniteFieldMatrices.jl.git
instantiate
```

The `instantiate` step will take a while as Julia's package manager downloads and compiles all dependencies.

Then, press backspace to go back to a julia prompt and run

```
using DeRham, Oscar
```

## Starting a new session after installation is complete

After starting a new Julia REPL, run

```
using Revise
] activate .
```

Then, press backspace to go back to a julia prompt and run

```
using DeRham, Oscar
```

## Troubleshooting

Please ensure that you're using the latest version of both Julia and Oscar. If the package manager gives you errors, try deleting Manifest.toml, upgrading Julia, updating the package registry, and updating Oscar.

## First computation

Build the Fermat cubic $x^3 + y^3 + z^3 = 0$ over $\mathbb{F}_7$, then compute its
middle-cohomology L-polynomial coefficients and zeta function:

```jldoctest
julia> using DeRham, Oscar

julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> coeffs = zeta_coefficients(f)
3-element Vector{ZZRingElem}:
 7
 1
 1

julia> L_polynomial(coeffs)
7*T^2 + T + 1

julia> zeta_function(coeffs, 7, 1)
(7*T^2 + T + 1)//(7*T^2 - 8*T + 1)
```
