# Square-root linear recurrence evaluation

The CPU DepthFirst and PChunk reduction policies can evaluate the linear
recurrence `v[k] = F(k) * v[k+1]` with affine matrices `F(k) = A*k + B`.
The square-root method uses a baby-step/giant-step product tree and fast
multipoint evaluation. The implementation uses exact polynomial and modular
arithmetic and remains experimental pending a general proof of correctness.

In [`zeta_coefficients`](@ref), set `algorithm=:depthfirst` or `algorithm=:pchunk`
alongside `linear_recurrence_algorithm=:sqrt`, `fastevaluation=true`, and a
positive `sqrt_linear_recurrence_min_length`. The cutoff is the smallest eligible
run length for square-root evaluation. Fast evaluation is enabled by default
for finite-field inputs. Pass
`fastevaluation=true` explicitly to `zeta_coefficients(f, p)`.

For DepthFirst reductions of projective hypersurfaces, omitting
`linear_recurrence_algorithm` or passing `nothing` inherits the supplied
pre-generated context's algorithm. Pass `:finitediff` explicitly to select
finite differences. Omitting the cutoff inherits the context's cutoff, and
passing `nothing` explicitly resets it to the default and selects finite
differences. PChunk uses the built-in context and accepts the recurrence
settings directly.

Unsupported matrix and vector types also use finite differences. The crossover
depends on the workload, so use representative timings to choose a cutoff when
enabling the experimental method.

```@docs
DeRham.sqrt_linear_recurrence!
```
