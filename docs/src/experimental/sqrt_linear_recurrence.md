# Square-root linear recurrence evaluation

The CPU DepthFirst and PChunk reduction policies can evaluate the linear
recurrence `v[k] = F(k) * v[k+1]` with affine matrices `F(k) = A*k + B`.
The square-root method uses a baby-step/giant-step product tree and fast
multipoint evaluation. The implementation uses exact polynomial and modular
arithmetic and remains experimental pending a general proof of correctness.

In [`zeta_coefficients`](@ref), set `algorithm=:depthfirst` or `algorithm=:pchunk`
alongside `linear_recurrence_algorithm=:sqrt` and `fastevaluation=true`.
Fast evaluation is enabled by default for finite-field inputs. Pass
`fastevaluation=true` explicitly to `zeta_coefficients(f, p)`.

With a pre-generated reduction context, omitted recurrence settings inherit
the context's algorithm and cutoff. Otherwise the defaults are `:finitediff`
and an automatic cutoff. Pass `sqrt_linear_recurrence_min_length=nothing`
explicitly to use the automatic cutoff, the greater of 64 and the square of
the matrix width. A positive integer sets a fixed minimum run length for
`:sqrt`. Unsupported matrix and vector types use finite differences.

The crossover still depends on the workload. Measure representative cases and
raise `sqrt_linear_recurrence_min_length` when the product tree remains slower.

```@docs
DeRham.sqrt_linear_recurrence!
```
