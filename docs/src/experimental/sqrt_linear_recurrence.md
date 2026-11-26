# Square-root linear recurrence evaluation

The CPU DepthFirst and PChunk reduction policies can evaluate the linear
recurrence `v[k] = F(k) * v[k+1]` with affine matrices `F(k) = A*k + B`.
The square-root method uses a baby-step/giant-step product tree and fast
multipoint evaluation. The implementation uses exact polynomial and modular
arithmetic and remains experimental pending a general proof of correctness.

Select it with `linear_recurrence_algorithm=:sqrt` in
[`zeta_coefficients`](@ref). The default is `:finitediff`. Use
`sqrt_linear_recurrence_min_length` to set the minimum run length for `:sqrt`. An explicit
positive integer is used as-is. When omitted, the automatic cutoff is the
maximum of 64 and the square of the matrix width, which delays the
cubic-cost polynomial matrix products for wider matrices. Unsupported matrix
and vector types continue through finite-difference evaluation.

The linear recurrence algorithm selection is used only when `fastevaluation=true` (the
default). With `fastevaluation=false`, the DepthFirst and PChunk reduction
loops use their stepwise path and ignore the linear recurrence settings.

The crossover still depends on the workload. Measure representative cases and
raise `sqrt_linear_recurrence_min_length` when the product tree remains slower.
