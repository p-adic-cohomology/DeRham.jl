# Square-root-time interval products

The CPU DepthFirst and PChunk reduction policies can evaluate a run of affine
matrix factors `F(k) = A*k + B` with a baby-step/giant-step product tree and
fast multipoint evaluation. The implementation uses exact polynomial and
modular arithmetic, but the method has not been proved correct for all inputs
and remains experimental.

Select it with `interval_product_algorithm=:sqrt` in
[`zeta_coefficients`](@ref). The default is `:finitediff`; set
`interval_product_min_length` to set the minimum run length. An explicit
positive integer is used as-is. When omitted, the automatic cutoff is the
maximum of 64 and the square of the matrix width, which delays the
cubic-cost polynomial matrix products for wider matrices. Unsupported matrix
and vector types continue through finite-difference evaluation.

The interval-product selection is used only when `fastevaluation=true` (the
default). With `fastevaluation=false`, the DepthFirst and PChunk reduction
loops use their stepwise path and ignore the interval-product settings.

The crossover still depends on the workload. Measure representative cases and
raise `interval_product_min_length` when the product tree remains slower.
