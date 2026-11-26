"""
    sqrt_linear_recurrence!(a, b, start, stop, g, temp, g_temp)

Evaluate the linear recurrence `v[x] = F(x) * v[x+1]` with the experimental
square-root-time Bostan–Gaudry–Schost method. For `F(x) = a*x + b`, overwrite
`g` with `F(start) * ... * F(stop) * g` using a product tree and fast multipoint
evaluation over the modulus ring. The inclusive indices and the run length must
fit in an `Int`.

Use square `zzModMatrix` inputs with `Vector{UInt}` vectors, or square
`ZZModMatrix` inputs with `Vector{ZZRingElem}` vectors. The scratch matrix
`temp` must match the input dimensions and modulus. Both vectors must have
one entry per row. Supply separate scratch storage for `temp` and `g_temp`.
Each `UInt` input entry must be less than the modulus. With `ZZRingElem`
vectors, initialize each input and scratch entry as a distinct object.
"""
function sqrt_linear_recurrence!(a, b, start, stop, g, temp, g_temp)
    start = Int(start)
    stop = Int(stop)
    start > stop && return g
    interval_length = Base.Checked.checked_add(Base.Checked.checked_sub(stop, start), 1)

    rows = nrows(a)
    cols = ncols(a)
    rows == cols || throw(ArgumentError("linear recurrence matrices must be square"))
    nrows(b) == rows && ncols(b) == cols ||
        throw(ArgumentError("linear recurrence matrices must have matching sizes"))

    small_modulus = a isa zzModMatrix
    big_modulus = a isa ZZModMatrix
    (small_modulus || big_modulus) &&
    _supports_sqrt_linear_recurrence(a, b, g, temp, g_temp) || throw(
        ArgumentError(
            "square-root linear recurrence requires compatible CPU modular arrays",
        ),
    )
    nrows(temp) == rows && ncols(temp) == cols ||
        throw(DimensionMismatch("linear recurrence temporary matrix must match the inputs"))
    length(g) == cols ||
        throw(DimensionMismatch("linear recurrence vector must match the matrix columns"))
    length(g_temp) == rows || throw(
        DimensionMismatch("linear recurrence temporary vector must match the matrix rows"),
    )
    if big_modulus && !all(i -> isassigned(g, i) && isassigned(g_temp, i), eachindex(g))
        throw(ArgumentError("linear recurrence BigInt vectors require initialized entries"))
    end

    ring = base_ring(parent(a))
    m = modulus(ring)
    m == modulus(base_ring(b)) == modulus(base_ring(temp)) ||
        throw(ArgumentError("linear recurrence matrices must share a modulus"))
    small_modulus &&
        !all(value -> value < m, g) &&
        throw(
            ArgumentError("linear recurrence UInt input entries must be reduced modulo $m"),
        )
    rows == 0 && return g

    backend = small_modulus ? _NmodBGSEvaluationBackend() : _FmpzModBGSEvaluationBackend()
    poly_ring, _ = polynomial_ring(ring, "x"; cached = false)
    polynomial_matrix_space = matrix_space(poly_ring, rows, cols)

    block_length = isqrt(interval_length)
    blocks, tail_length = divrem(interval_length, block_length)
    factors = [polynomial_matrix_space() for _ = 1:block_length]
    for offset = 0:(block_length-1)
        factor = factors[offset+1]
        ring_offset = ring(offset)
        for j = 1:cols, i = 1:rows
            linear = ring(a[i, j])
            constant = ring(b[i, j]) + linear * ring_offset
            factor[i, j] = poly_ring([constant, linear])
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
    evaluated_blocks = [matrix_space(ring, rows, cols)() for _ = 1:blocks]
    evaluation_buffer = _bgs_evaluation_buffer(backend, blocks)
    evaluation_points =
        _bgs_evaluation_points(backend, ring, points, block_polynomial[1, 1])
    try
        for j = 1:cols, i = 1:rows
            values = _evaluate_bgs_polynomial!(
                evaluation_buffer,
                block_polynomial[i, j],
                evaluation_points,
            )
            for k = 1:blocks
                _set_bgs_entry!(backend, evaluated_blocks[k], ring, values[k], i, j)
            end
        end
    finally
        _release_bgs_evaluation_points!(backend, evaluation_points)
    end

    for offset = 0:(tail_length-1)
        my_mul!(temp, a, stop - offset)
        my_add!(temp, temp, b)
        my_matvecmul!(g_temp, temp, g)
        my_copy!(g, g_temp)
    end
    for k = blocks:-1:1
        my_matvecmul!(g_temp, evaluated_blocks[k], g)
        my_copy!(g, g_temp)
    end
    return g
end

"""
Modular polynomial evaluation backend for the BGS method. The concrete
implementations manage the FLINT calls and raw value storage.
"""
abstract type _BGSEvaluationBackend end

struct _NmodBGSEvaluationBackend <: _BGSEvaluationBackend end

struct _FmpzModBGSEvaluationBackend <: _BGSEvaluationBackend end

_set_bgs_entry!(::_NmodBGSEvaluationBackend, matrix, ring, value, i, j) =
    setindex!(matrix, ring(value), i, j)

_set_bgs_entry!(::_FmpzModBGSEvaluationBackend, matrix, _, value, i, j) =
    Oscar.Nemo.setindex_raw!(matrix, value, i, j)

struct _NmodBGSModulus
    n::UInt
    ninv::UInt
    norm::UInt
end

struct _NmodBGSEvaluationPoints
    tree::Ptr{Ptr{UInt}}
    count::Int
    modulus::_NmodBGSModulus
end

struct _FmpzModBGSEvaluationPoints
    tree::Ptr{Ptr{Cvoid}}
    count::Int
end

function _bgs_evaluation_points(
    ::_NmodBGSEvaluationBackend,
    ring,
    points::Vector{Int},
    poly::Oscar.Nemo.zzModPolyRingElem,
)
    values = UInt[data(ring(point)) for point in points]
    modulus = _NmodBGSModulus(poly.mod_n, poly.mod_ninv, poly.mod_norm)
    isempty(values) && return _NmodBGSEvaluationPoints(C_NULL, 0, modulus)

    tree = @ccall Oscar.Nemo.libflint._nmod_poly_tree_alloc(
        length(values)::Int,
    )::Ptr{Ptr{UInt}}
    tree == C_NULL && throw(OutOfMemoryError())
    try
        GC.@preserve values modulus begin
            @ccall Oscar.Nemo.libflint._nmod_poly_tree_build(
                tree::Ptr{Ptr{UInt}},
                values::Ptr{UInt},
                length(values)::Int,
                modulus::_NmodBGSModulus,
            )::Nothing
        end
    catch
        @ccall Oscar.Nemo.libflint._nmod_poly_tree_free(
            tree::Ptr{Ptr{UInt}},
            length(values)::Int,
        )::Nothing
        rethrow()
    end
    return _NmodBGSEvaluationPoints(tree, length(values), modulus)
end

function _bgs_evaluation_points(
    ::_FmpzModBGSEvaluationBackend,
    ring,
    points::Vector{Int},
    ::Oscar.Nemo.ZZModPolyRingElem,
)
    isempty(points) && return _FmpzModBGSEvaluationPoints(C_NULL, 0)
    reduced_points = ZZRingElem[lift(ring(point)) for point in points]
    raw_points = zeros(Int, length(points))
    tree = Ptr{Ptr{Cvoid}}(C_NULL)
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
        tree = @ccall Oscar.Nemo.libflint._fmpz_mod_poly_tree_alloc(
            length(points)::Int,
        )::Ptr{Ptr{Cvoid}}
        tree == C_NULL && throw(OutOfMemoryError())
        GC.@preserve ring raw_points tree begin
            @ccall Oscar.Nemo.libflint._fmpz_mod_poly_tree_build(
                tree::Ptr{Ptr{Cvoid}},
                raw_points::Ptr{Int},
                length(points)::Int,
                ring.ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct},
            )::Nothing
        end
    catch
        if tree != C_NULL
            @ccall Oscar.Nemo.libflint._fmpz_mod_poly_tree_free(
                tree::Ptr{Ptr{Cvoid}},
                length(points)::Int,
            )::Nothing
        end
        rethrow()
    finally
        _clear_fmpz_values!(raw_points, initialized)
    end
    return _FmpzModBGSEvaluationPoints(tree, length(points))
end

_bgs_evaluation_buffer(::_NmodBGSEvaluationBackend, count::Int) = Vector{UInt}(undef, count)

_bgs_evaluation_buffer(::_FmpzModBGSEvaluationBackend, count::Int) =
    (zeros(Int, count), ZZRingElem[ZZ(0) for _ = 1:count])

function _release_bgs_evaluation_points!(
    ::_NmodBGSEvaluationBackend,
    points::_NmodBGSEvaluationPoints,
)
    if points.tree != C_NULL
        @ccall Oscar.Nemo.libflint._nmod_poly_tree_free(
            points.tree::Ptr{Ptr{UInt}},
            points.count::Int,
        )::Nothing
    end
    return nothing
end

function _release_bgs_evaluation_points!(
    ::_FmpzModBGSEvaluationBackend,
    points::_FmpzModBGSEvaluationPoints,
)
    if points.tree != C_NULL
        @ccall Oscar.Nemo.libflint._fmpz_mod_poly_tree_free(
            points.tree::Ptr{Ptr{Cvoid}},
            points.count::Int,
        )::Nothing
    end
    return nothing
end

function _evaluate_bgs_polynomial!(
    values::Vector{UInt},
    poly::Oscar.Nemo.zzModPolyRingElem,
    points::_NmodBGSEvaluationPoints,
)
    length(values) == points.count ||
        throw(DimensionMismatch("evaluation buffer must match the evaluation points"))
    points.count == 0 && return values
    if poly.length == 0
        fill!(values, 0)
        return values
    end
    GC.@preserve poly points values begin
        @ccall Oscar.Nemo.libflint._nmod_poly_evaluate_nmod_vec_fast_precomp(
            values::Ptr{UInt},
            poly.coeffs::Ptr{Nothing},
            poly.length::Int,
            points.tree::Ptr{Ptr{UInt}},
            points.count::Int,
            points.modulus::_NmodBGSModulus,
        )::Nothing
    end
    return values
end

function _evaluate_bgs_polynomial!(
    buffers::Tuple{Vector{Int},Vector{ZZRingElem}},
    poly::Oscar.Nemo.ZZModPolyRingElem,
    points::_FmpzModBGSEvaluationPoints,
)
    raw_values, values = buffers
    length(values) == points.count && length(raw_values) == points.count ||
        throw(DimensionMismatch("evaluation buffers must match the evaluation points"))
    points.count == 0 && return values
    if poly.length == 0
        GC.@preserve values begin
            for value in values
                @ccall Oscar.Nemo.libflint.fmpz_set_si(
                    value::Ref{ZZRingElem},
                    0::Int,
                )::Nothing
            end
        end
        return values
    end
    ring = base_ring(parent(poly))
    GC.@preserve poly points raw_values values ring begin
        try
            @ccall Oscar.Nemo.libflint._fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(
                raw_values::Ptr{Int},
                poly.coeffs::Ptr{Nothing},
                poly.length::Int,
                points.tree::Ptr{Ptr{Cvoid}},
                points.count::Int,
                ring.ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct},
            )::Nothing
            for i in eachindex(values)
                value = values[i]
                source = Ptr{ZZRingElem}(pointer(raw_values, i))
                @ccall Oscar.Nemo.libflint.fmpz_set(
                    value::Ref{ZZRingElem},
                    source::Ptr{ZZRingElem},
                )::Nothing
            end
        finally
            _clear_fmpz_values!(raw_values, length(raw_values))
            fill!(raw_values, 0)
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

function _supports_sqrt_linear_recurrence(a, b, g, temp, g_temp)
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

function _use_sqrt_linear_recurrence(
    a,
    b,
    g,
    temp,
    g_temp,
    algorithm,
    configured_min,
    interval_length,
)
    algorithm === :sqrt && configured_min !== nothing || return false
    interval_length >= configured_min || return false
    return _supports_sqrt_linear_recurrence(a, b, g, temp, g_temp)
end

function evaluate_linear_recurrence!(
    a,
    b,
    start,
    stop,
    g,
    temp,
    g_temp,
    algorithm::Symbol,
    configured_min::Union{Nothing,Int},
)
    interval_length = stop - start + 1
    if _use_sqrt_linear_recurrence(
        a,
        b,
        g,
        temp,
        g_temp,
        algorithm,
        configured_min,
        interval_length,
    )
        @debug "Using square-root linear recurrence"
        return sqrt_linear_recurrence!(a, b, start, stop, g, temp, g_temp)
    end
    return finitediff_prodeval_linear!(a, b, start, stop, g, temp, g_temp)
end
