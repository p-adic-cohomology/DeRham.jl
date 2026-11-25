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
    dense_reduction::Union{Nothing,Symbol}
    handoff_pole_order::Union{Nothing,Int}
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
    nothing,
    nothing,
)

default_params() =
    ZetaFunctionParams(0, false, :default, :invlex, false, true, false, false, false)
