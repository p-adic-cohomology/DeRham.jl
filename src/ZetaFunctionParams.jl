struct ZetaFunctionParams
    verbose::Int
    givefrobmat::Bool
    algorithm::Symbol
    termorder::Symbol
    vars_reversed::Bool
    fastevaluation::Bool
    always_use_bigints::Bool
    use_gpu::Bool
    use_threads::Bool
    interval_product_algorithm::Symbol
    interval_product_min_length::Int
end

ZetaFunctionParams(
    verbose,
    givefrobmat,
    algorithm,
    termorder,
    vars_reversed,
    fastevaluation,
    always_use_bigints,
    use_gpu,
    use_threads,
) = ZetaFunctionParams(
    verbose,
    givefrobmat,
    algorithm,
    termorder,
    vars_reversed,
    fastevaluation,
    always_use_bigints,
    use_gpu,
    use_threads,
    :finitediff,
    64,
)

default_params() =
    ZetaFunctionParams(0, false, :default, :invlex, false, true, false, false, false)
