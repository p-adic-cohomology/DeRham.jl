# Square-root-time interval products

The CPU DepthFirst and PChunk reduction policies can evaluate a run of affine
matrix factors `F(k) = A*k + B` with a baby-step/giant-step product tree and
fast multipoint evaluation. The implementation uses exact polynomial and
modular arithmetic, but the method has not been proved correct for all inputs
and remains experimental.

Select it with `interval_product_algorithm=:sqrt` in
[`zeta_coefficients`](@ref). The default is `:finitediff`; set
`interval_product_min_length` to control the shortest run that uses the new
method. Unsupported matrix and vector types continue through finite-difference
evaluation.

The crossover depends on matrix width. The default minimum is fixed, so measure
representative workloads and raise `interval_product_min_length` for wider
matrices when the product tree is slower.
