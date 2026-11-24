"""
    sqrt_interval_product!(a, b, start, stop, g, temp, g_temp)

Experimental square-root-time interval product for the affine matrices
`F(x) = a*x + b`. It computes `F(start) * ... * F(stop) * g` using a product
tree and fast multipoint evaluation over the modulus ring.
"""
function sqrt_interval_product!(a, b, start, stop, g, temp, g_temp)
    start = Int(start)
    stop = Int(stop)
    interval_length = stop - start + 1
    interval_length <= 0 && return g

    rows = nrows(a)
    cols = ncols(a)
    rows == cols || throw(ArgumentError("interval product matrices must be square"))
    nrows(b) == rows && ncols(b) == cols ||
        throw(ArgumentError("interval product matrices must have matching sizes"))

    small_modulus = a isa zzModMatrix
    big_modulus = a isa ZZModMatrix
    (small_modulus || big_modulus) &&
    _supports_sqrt_interval_product(a, b, g, temp, g_temp) || throw(
        ArgumentError(
            "square-root interval products require compatible CPU modular arrays",
        ),
    )

    backend =
        small_modulus ? _NmodIntervalProductEvaluationBackend() :
        _FmpzModIntervalProductEvaluationBackend()
    ring = base_ring(parent(a))
    poly_ring, x = polynomial_ring(ring, "x"; cached = false)
    polynomial_matrix_space = matrix_space(poly_ring, rows, cols)

    block_length = isqrt(interval_length)
    blocks, _ = divrem(interval_length, block_length)
    factors = [polynomial_matrix_space() for _ = 1:block_length]
    for offset = 0:(block_length-1)
        factor = factors[offset+1]
        ring_offset = ring(offset)
        for j = 1:cols, i = 1:rows
            linear = ring(a[i, j])
            constant = ring(b[i, j]) + linear * ring_offset
            factor[i, j] = poly_ring(linear) * x + poly_ring(constant)
        end
    end

    while Base.length(factors) > 1
        next_level = typeof(factors)()
        for i = 1:2:Base.length(factors)
            if i == Base.length(factors)
                push!(next_level, factors[i])
            else
                push!(next_level, factors[i] * factors[i+1])
            end
        end
        factors = next_level
    end
    block_polynomial = factors[1]

    points = Int[start + i * block_length for i = 0:(blocks-1)]
    evaluation_points = _interval_product_evaluation_points(backend, ring, points)
    evaluated_blocks = [matrix_space(ring, rows, cols)() for _ = 1:blocks]
    try
        for j = 1:cols, i = 1:rows
            values = _evaluate_interval_product_polynomial(
                backend,
                block_polynomial[i, j],
                evaluation_points,
            )
            for k = 1:blocks
                evaluated_blocks[k][i, j] = ring(values[k])
            end
        end
    finally
        _release_interval_product_evaluation_points!(backend, evaluation_points)
    end

    product = evaluated_blocks[1]
    for k = 2:blocks
        product = product * evaluated_blocks[k]
    end

    tail_start = start + blocks * block_length
    for value = stop:-1:tail_start
        my_mul!(temp, a, value)
        my_add!(temp, temp, b)
        my_matvecmul!(g_temp, temp, g)
        my_copy!(g, g_temp)
    end
    my_matvecmul!(g_temp, product, g)
    my_copy!(g, g_temp)
    return g
end

"""
Backend boundary for modular polynomial evaluation. The interval algorithm uses
the backend methods below; concrete implementations own FLINT calls and raw
value storage.
"""
abstract type _IntervalProductEvaluationBackend end

struct _NmodIntervalProductEvaluationBackend <: _IntervalProductEvaluationBackend end

struct _FmpzModIntervalProductEvaluationBackend <: _IntervalProductEvaluationBackend end

struct _FmpzModEvaluationPoints
    values::Vector{Int}
end

_interval_product_evaluation_points(
    ::_NmodIntervalProductEvaluationBackend,
    ring,
    points::Vector{Int},
) = UInt[data(ring(point)) for point in points]

function _interval_product_evaluation_points(
    ::_FmpzModIntervalProductEvaluationBackend,
    ring,
    points::Vector{Int},
)
    reduced_points = ZZRingElem[lift(ring(point)) for point in points]
    raw_points = zeros(Int, length(points))
    initialized = 0
    try
        GC.@preserve ring reduced_points raw_points begin
            for i in eachindex(points)
                point = reduced_points[i]
                destination = Ptr{ZZRingElem}(pointer(raw_points, i))
                @ccall Oscar.Nemo.libflint.fmpz_set(
                    destination::Ptr{ZZRingElem},
                    point::Ref{ZZRingElem},
                )::Nothing
                initialized = i
            end
        end
    catch
        _clear_fmpz_values!(raw_points, initialized)
        rethrow()
    end
    return _FmpzModEvaluationPoints(raw_points)
end

_evaluate_interval_product_polynomial(
    ::_NmodIntervalProductEvaluationBackend,
    poly::Oscar.Nemo.zzModPolyRingElem,
    points::Vector{UInt},
) = _fast_nmod_evaluate(poly, points)

_evaluate_interval_product_polynomial(
    ::_FmpzModIntervalProductEvaluationBackend,
    poly::Oscar.Nemo.ZZModPolyRingElem,
    points::_FmpzModEvaluationPoints,
) = _fast_fmpz_mod_evaluate(poly, points)

_release_interval_product_evaluation_points!(
    ::_NmodIntervalProductEvaluationBackend,
    points,
) = nothing

function _release_interval_product_evaluation_points!(
    ::_FmpzModIntervalProductEvaluationBackend,
    points::_FmpzModEvaluationPoints,
)
    _clear_fmpz_values!(points.values, length(points.values))
    return nothing
end

function _fast_nmod_evaluate(poly::Oscar.Nemo.zzModPolyRingElem, points::Vector{UInt})
    values = Vector{UInt}(undef, length(points))
    GC.@preserve poly points values begin
        @ccall Oscar.Nemo.libflint.nmod_poly_evaluate_nmod_vec_fast(
            values::Ptr{UInt},
            poly::Ref{Oscar.Nemo.zzModPolyRingElem},
            points::Ptr{UInt},
            length(points)::Int,
        )::Nothing
    end
    return values
end

function _fast_fmpz_mod_evaluate(poly::Oscar.Nemo.ZZModPolyRingElem, points::Vector{Int})
    backend = _FmpzModIntervalProductEvaluationBackend()
    evaluation_points =
        _interval_product_evaluation_points(backend, base_ring(parent(poly)), points)
    try
        return _fast_fmpz_mod_evaluate(poly, evaluation_points)
    finally
        _release_interval_product_evaluation_points!(backend, evaluation_points)
    end
end

function _fast_fmpz_mod_evaluate(
    poly::Oscar.Nemo.ZZModPolyRingElem,
    points::_FmpzModEvaluationPoints,
)
    raw_points = points.values
    raw_values = zeros(Int, length(raw_points))
    values = Vector{ZZRingElem}(undef, length(raw_points))
    ring = base_ring(parent(poly))
    GC.@preserve poly points raw_points raw_values values ring begin
        try
            @ccall Oscar.Nemo.libflint.fmpz_mod_poly_evaluate_fmpz_vec_fast(
                raw_values::Ptr{Int},
                poly::Ref{Oscar.Nemo.ZZModPolyRingElem},
                raw_points::Ptr{Int},
                length(raw_points)::Int,
                ring.ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct},
            )::Nothing
            for i in eachindex(raw_values)
                value = ZZ(0)
                source = Ptr{ZZRingElem}(pointer(raw_values, i))
                @ccall Oscar.Nemo.libflint.fmpz_set(
                    value::Ref{ZZRingElem},
                    source::Ptr{ZZRingElem},
                )::Nothing
                values[i] = value
            end
        finally
            _clear_fmpz_values!(raw_values, length(raw_values))
        end
    end
    return values
end

function _clear_fmpz_values!(values::Vector{Int}, count::Int)
    GC.@preserve values begin
        for i = 1:count
            value = Ptr{ZZRingElem}(pointer(values, i))
            @ccall Oscar.Nemo.libflint.fmpz_clear(value::Ptr{ZZRingElem})::Nothing
        end
    end
    return nothing
end

function _supports_sqrt_interval_product(a, b, g, temp, g_temp)
    if a isa zzModMatrix
        return b isa zzModMatrix &&
               temp isa zzModMatrix &&
               g isa Vector{UInt} &&
               g_temp isa Vector{UInt}
    elseif a isa ZZModMatrix
        return b isa ZZModMatrix &&
               temp isa ZZModMatrix &&
               g isa Vector{ZZRingElem} &&
               g_temp isa Vector{ZZRingElem}
    end
    return false
end

function _sqrt_interval_product_min_length(
    matrix_width::Int,
    configured_min::Union{Nothing,Int},
)
    configured_min === nothing || return configured_min
    width_squared = if matrix_width > 0 && matrix_width > div(typemax(Int), matrix_width)
        typemax(Int)
    else
        matrix_width^2
    end
    return max(64, width_squared)
end

function _with_interval_product_test_hook(f::Function, hook::Function)
    storage = task_local_storage()
    key = :DeRham_interval_product_test_hook
    previous = get(storage, key, nothing)
    storage[key] = hook
    try
        return f()
    finally
        if previous === nothing
            delete!(storage, key)
        else
            storage[key] = previous
        end
    end
end

function _notify_interval_product_test_hook()
    hook = get(task_local_storage(), :DeRham_interval_product_test_hook, nothing)
    hook === nothing || hook()
    return nothing
end

function _use_sqrt_interval_product(a, b, g, temp, g_temp, params, interval_length)
    params.interval_product_algorithm === :sqrt || return false
    _supports_sqrt_interval_product(a, b, g, temp, g_temp) || return false
    min_length =
        _sqrt_interval_product_min_length(nrows(a), params.interval_product_min_length)
    return interval_length >= min_length
end

function prodeval_linear!(a, b, start, stop, g, temp, g_temp, params::ZetaFunctionParams)
    interval_length = stop - start + 1
    if _use_sqrt_interval_product(a, b, g, temp, g_temp, params, interval_length)
        _notify_interval_product_test_hook()
        return sqrt_interval_product!(a, b, start, stop, g, temp, g_temp)
    end
    return finitediff_prodeval_linear!(a, b, start, stop, g, temp, g_temp)
end
