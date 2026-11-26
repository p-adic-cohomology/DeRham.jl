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
    linear_recurrence_algorithm::Symbol
    sqrt_linear_recurrence_min_length::Union{Nothing,Int}
end

const _UNSPECIFIED_SQRT_LINEAR_RECURRENCE_MIN_LENGTH = Val{:unspecified}()

function _validate_linear_recurrence_options(algorithm, min_length)
    (algorithm === :finitediff || algorithm === :sqrt) ||
        throw(ArgumentError("linear_recurrence_algorithm must be :finitediff or :sqrt"))
    if min_length !== nothing
        min_length isa Integer && !(min_length isa Bool) && min_length > 0 || throw(
            ArgumentError(
                "sqrt_linear_recurrence_min_length must be a positive integer or nothing",
            ),
        )
        min_length = try
            Int(min_length)
        catch error
            error isa InexactError || rethrow()
            throw(ArgumentError("sqrt_linear_recurrence_min_length must fit in an Int"))
        end
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
    nothing,
)

default_params() =
    ZetaFunctionParams(0, false, :default, :invlex, false, true, false, false, false)
