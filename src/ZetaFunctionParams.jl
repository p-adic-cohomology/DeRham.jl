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

function _validate_interval_product_options(algorithm, min_length)
    algorithm in (:finitediff, :sqrt) ||
        throw(ArgumentError("interval_product_algorithm must be :finitediff or :sqrt"))
    min_length isa Integer && !(min_length isa Bool) && min_length > 0 ||
        throw(ArgumentError("interval_product_min_length must be a positive integer"))
    min_length = try
        Int(min_length)
    catch error
        error isa InexactError || rethrow()
        throw(ArgumentError("interval_product_min_length must fit in an Int"))
    end
    return algorithm, min_length
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
