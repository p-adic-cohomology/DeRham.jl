
"""
    ControlledReductionContext{MatrixType,VectorType}

A controlled reduction context is used for in-place evaluation
chunks of controlled reduction.
It contains all of the major allocations necessary, so you can
create one context at the beginning (or one per thread)
and then not have to worry about allocating anything else.

MatrixType - the type of A and B, for example zzModMatrix
VectorType - the type of g, for exampel Vector{UInt64}

fields
------
Ruvs - an AbstractPEP which containes the R_uv necessary for the reduction.
  This is intended to be shared across all contexts (i.e. all threads)
  for a single example.
A - the matrix which will store A in Ay + B
B - the matrix which will store B in Ay + B
temp - temporary matrix used for intermediate calculations
g - the vector g which contains the vector which is being reduced
g_temp - temporary vector used for intermediate calculations
linear_recurrence_algorithm - default linear recurrence algorithm for reductions
sqrt_linear_recurrence_min_length - minimum run length for the square-root method,
  or `nothing` for the automatic matrix-width-based cutoff

Current examples that are in use:

ControlledReductionContext{zzModMatrix,Vector{UInt64}}
ControlledReductionContext{ZZModMatrix,Vector{ZZRingElem}}
ControlledReductionContext{CuModMatrix,CuModVector}

"""
struct ControlledReductionContext{MatrixType,VectorType}
    Ruvs::AbstractPEP{MatrixType}
    A::MatrixType
    B::MatrixType
    temp::MatrixType
    g::VectorType
    g_temp::VectorType
    linear_recurrence_algorithm::Symbol
    sqrt_linear_recurrence_min_length::Union{Nothing,Int}
end

ControlledReductionContext(Ruv, A, B, temp, g, g_temp) =
    ControlledReductionContext(Ruv, A, B, temp, g, g_temp, :finitediff, nothing)

function ControlledReductionContext{MatrixType,VectorType}(
    Ruv,
    A,
    B,
    temp,
    g,
    g_temp,
) where {MatrixType,VectorType}
    return ControlledReductionContext{MatrixType,VectorType}(
        Ruv,
        A,
        B,
        temp,
        g,
        g_temp,
        :finitediff,
        nothing,
    )
end

function _resolve_linear_recurrence_options(context, algorithm, min_length)
    if min_length === _UNSPECIFIED_SQRT_LINEAR_RECURRENCE_MIN_LENGTH
        min_length =
            context isa ControlledReductionContext ?
            context.sqrt_linear_recurrence_min_length : nothing
    end
    if context isa ControlledReductionContext
        algorithm = something(algorithm, context.linear_recurrence_algorithm)
    else
        algorithm = something(algorithm, :finitediff)
    end
    return _validate_linear_recurrence_options(algorithm, min_length)
end

"""
    PreconfiguredReductionEngine(context)

Wrap reusable controlled-reduction buffers and R_uv storage for constructing a
matrix engine for one frontend polynomial.
"""
struct PreconfiguredReductionEngine{C<:ControlledReductionContext}
    context::C
end

function default_context(matspace, Ruv, params)
    g_length = number_of_rows(matspace)
    B = base_ring(matspace)
    m = modulus(B)
    if params.use_gpu == true && (ZZ(2)^25 < m < ZZ(2)^106)
        println("using karatsuba")
        temp = collect(factor(m))[1]
        d = temp[2]
        p = temp[1]

        N1 = Int(p)^Int(round(d/2))
        N2 = Int(p)^Int((d - round(d/2)))

        A = GPUFiniteFieldMatrices.KaratsubaZeros(
            Float64,
            g_length,
            g_length,
            N1,
            N2,
            N1*N2,
            true,
        )
        B = GPUFiniteFieldMatrices.KaratsubaZeros(
            Float64,
            g_length,
            g_length,
            N1,
            N2,
            N1*N2,
            true,
        )
        temp = GPUFiniteFieldMatrices.KaratsubaZeros(
            Float64,
            g_length,
            g_length,
            N1,
            N2,
            N1*N2,
            true,
        )
        GPUFiniteFieldMatrices.initialize_plan!(A)
        GPUFiniteFieldMatrices.initialize_plan!(B)
        GPUFiniteFieldMatrices.initialize_plan!(temp)

        g = GPUFiniteFieldMatrices.KaratsubaZeros(Float64, g_length, N1, N2, N1*N2, true)
        g_temp =
            GPUFiniteFieldMatrices.KaratsubaZeros(Float64, g_length, N1, N2, N1*N2, true)
        GPUFiniteFieldMatrices.initialize_plan!(g)
        GPUFiniteFieldMatrices.initialize_plan!(g_temp)

        return ControlledReductionContext{KaratsubaMatrix{Float64},KaratsubaVector{Float64}}(
            Ruv,
            A,
            B,
            temp,
            g,
            g_temp,
            params.linear_recurrence_algorithm,
            params.sqrt_linear_recurrence_min_length,
        )
    elseif params.use_gpu == true

        M = Int(m)
        A = GPUFiniteFieldMatrices.zeros(Float64, g_length, g_length, M)
        B = GPUFiniteFieldMatrices.zeros(Float64, g_length, g_length, M)
        temp = GPUFiniteFieldMatrices.zeros(Float64, g_length, g_length, M)

        g = GPUFiniteFieldMatrices.zeros(Float64, g_length, M)
        g_temp = GPUFiniteFieldMatrices.zeros(Float64, g_length, M)

        return ControlledReductionContext{CuModMatrix{Float64},CuModVector{Float64}}(
            Ruv,
            A,
            B,
            temp,
            g,
            g_temp,
            params.linear_recurrence_algorithm,
            params.sqrt_linear_recurrence_min_length,
        )
    else # use the CPU
        A = matspace()
        B = matspace()
        temp = matspace()

        if params.always_use_bigints || ZZ(2)^64 <= ZZ(m)
            # Big modulus
            # NOTE: can't use `similar` for a pointer type
            g = [ZZ(0) for i = 1:g_length]
            g_temp = [ZZ(0) for i = 1:g_length]
        else
            # Small modulus
            g = zeros(UInt, g_length)
            g_temp = similar(g)
        end

        return ControlledReductionContext(
            Ruv,
            A,
            B,
            temp,
            g,
            g_temp,
            params.linear_recurrence_algorithm,
            params.sqrt_linear_recurrence_min_length,
        )
    end
end

"""
    pregen_default_context(n, d, p, S; kwargs...)

Build a precomputed reduction context for the parameters `n`, `d`, `p`, and
`S`. Here, `n` is the number of variables minus one, `d` is the degree, `p`
is the prime characteristic, and `S` is a vector of integers.

Keyword arguments and defaults:
`verbose = 0`, `givefrobmat = false`, `algorithm = :naive`,
`termorder = :invlex`, `vars_reversed = false`, `fastevaluation = false`,
`always_use_bigints = false`, `use_gpu = false`, `use_threads = false`,
`linear_recurrence_algorithm = :finitediff`,
`sqrt_linear_recurrence_min_length = nothing`, and `lazy = false`.

When `sqrt_linear_recurrence_min_length` is `nothing`, the automatic cutoff is the
maximum of 64 and the square of the current matrix width. An explicit positive
integer is used as-is.

The linear recurrence settings are retained in the context and supply the
defaults for `zeta_coefficients` when a context is passed and those keywords are
omitted.
"""
function pregen_default_context(
    n,
    d,
    p,
    S;
    verbose = 0,
    givefrobmat = false,
    algorithm = :naive,
    termorder = :invlex,
    vars_reversed = false,
    fastevaluation = false,
    always_use_bigints = false,
    use_gpu = false,
    use_threads = false,
    linear_recurrence_algorithm = :finitediff,
    sqrt_linear_recurrence_min_length = nothing,
    lazy = false,
)
    linear_recurrence_algorithm, sqrt_linear_recurrence_min_length =
        _validate_linear_recurrence_options(
            linear_recurrence_algorithm,
            sqrt_linear_recurrence_min_length,
        )
    params = ZetaFunctionParams(
        verbose,
        givefrobmat,
        algorithm,
        termorder,
        vars_reversed,
        fastevaluation,
        always_use_bigints,
        use_gpu,
        use_threads,
        linear_recurrence_algorithm,
        sqrt_linear_recurrence_min_length,
    )

    m = pregen_precision_info(n, d, p)
    M = Integer(p^m)
    g_length = binomial(d*n, d*n-n)
    matspace = matrix_space(residue_ring(ZZ, M)[1], g_length, g_length)

    Ruv = pregen_select_Ruv_PEP(n, d, S, params, matspace)

    if params.use_gpu == true && (ZZ(2)^25 < M < ZZ(2)^106)
        println("using karatsuba")

        N1 = Int(p)^Int(round(m/2))
        N2 = Int(p)^Int((m - round(m/2)))

        A = GPUFiniteFieldMatrices.KaratsubaZeros(
            Float64,
            g_length,
            g_length,
            N1,
            N2,
            N1*N2,
            true,
        )
        B = GPUFiniteFieldMatrices.KaratsubaZeros(
            Float64,
            g_length,
            g_length,
            N1,
            N2,
            N1*N2,
            true,
        )
        temp = GPUFiniteFieldMatrices.KaratsubaZeros(
            Float64,
            g_length,
            g_length,
            N1,
            N2,
            N1*N2,
            true,
        )
        GPUFiniteFieldMatrices.initialize_plan!(A)
        GPUFiniteFieldMatrices.initialize_plan!(B)
        GPUFiniteFieldMatrices.initialize_plan!(temp)

        g = GPUFiniteFieldMatrices.KaratsubaZeros(Float64, g_length, N1, N2, N1*N2, true)
        g_temp =
            GPUFiniteFieldMatrices.KaratsubaZeros(Float64, g_length, N1, N2, N1*N2, true)
        GPUFiniteFieldMatrices.initialize_plan!(g)
        GPUFiniteFieldMatrices.initialize_plan!(g_temp)

        return ControlledReductionContext{KaratsubaMatrix{Float64},KaratsubaVector{Float64}}(
            Ruv,
            A,
            B,
            temp,
            g,
            g_temp,
            params.linear_recurrence_algorithm,
            params.sqrt_linear_recurrence_min_length,
        )
    elseif params.use_gpu == true

        A = GPUFiniteFieldMatrices.zeros(Float64, g_length, g_length, M)
        B = GPUFiniteFieldMatrices.zeros(Float64, g_length, g_length, M)
        temp = GPUFiniteFieldMatrices.zeros(Float64, g_length, g_length, M)

        g = GPUFiniteFieldMatrices.zeros(Float64, g_length, M)
        g_temp = GPUFiniteFieldMatrices.zeros(Float64, g_length, M)

        return ControlledReductionContext{CuModMatrix{Float64},CuModVector{Float64}}(
            Ruv,
            A,
            B,
            temp,
            g,
            g_temp,
            params.linear_recurrence_algorithm,
            params.sqrt_linear_recurrence_min_length,
        )
    else # use the CPU
        A = matspace()
        B = matspace()
        temp = matspace()

        if params.always_use_bigints || ZZ(2)^64 <= ZZ(M)
            # Big modulus
            # NOTE: can't use `similar` for a pointer type
            g = [ZZ(0) for i = 1:g_length]
            g_temp = [ZZ(0) for i = 1:g_length]
        else
            # Small modulus
            g = zeros(UInt, g_length)
            g_temp = similar(g)
        end

        return ControlledReductionContext(
            Ruv,
            A,
            B,
            temp,
            g,
            g_temp,
            params.linear_recurrence_algorithm,
            params.sqrt_linear_recurrence_min_length,
        )
    end
end

function default_context_type(matspace, params)
    B = base_ring(matspace)
    m = modulus(B)

    if params.use_gpu == true
        ControlledReductionContext{CuModMatrix{Float64},CuModVector{Float64}}
    elseif params.always_use_bigints || ZZ(2)^64 <= ZZ(m)
        ControlledReductionContext{ZZModMatrix,Vector{ZZRingElem}}
    else
        ControlledReductionContext{zzModMatrix,Vector{UInt}}
    end
end
