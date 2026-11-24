# Square-root-time interval products

The CPU DepthFirst and PChunk reduction policies can evaluate a run of affine
matrix factors `F(k) = A*k + B` with a baby-step/giant-step product tree and
fast multipoint evaluation. The implementation uses exact polynomial and
modular arithmetic, but the method has not been proved correct for all inputs
and remains experimental.

Select it with `interval_product_algorithm=:sqrt` in
[`zeta_coefficients`](@ref). The default is `:finitediff`; set
`interval_product_min_length` to set a lower bound on the run length. The
dispatcher also requires the interval length to be at least the square of the
matrix width, which delays the cubic-cost polynomial matrix products for wider
matrices. Unsupported matrix and vector types continue through
finite-difference evaluation.

The crossover still depends on the workload. Measure representative cases and
raise `interval_product_min_length` when the product tree remains slower.
