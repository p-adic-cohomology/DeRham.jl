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
        _supports_sqrt_interval_product(a, b, g, temp, g_temp) ||
        throw(ArgumentError("square-root interval products require compatible CPU modular arrays"))

    ring = base_ring(parent(a))
    poly_ring, x = polynomial_ring(ring, "x"; cached = false)
    polynomial_matrix_space = matrix_space(poly_ring, rows, cols)

    block_length = isqrt(interval_length)
    blocks, _ = divrem(interval_length, block_length)
    factors = [polynomial_matrix_space() for _ = 1:block_length]
    for offset = 0:(block_length - 1)
        factor = factors[offset + 1]
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
                push!(next_level, factors[i] * factors[i + 1])
            end
        end
        factors = next_level
    end
    block_polynomial = factors[1]

    points = Int[start + i * block_length for i = 0:(blocks - 1)]
    evaluated_blocks = [matrix_space(ring, rows, cols)() for _ = 1:blocks]
    for j = 1:cols, i = 1:rows
        if small_modulus
            point_data = UInt[data(ring(point)) for point in points]
            values = _fast_nmod_evaluate(block_polynomial[i, j], point_data)
        else
            values = _fast_fmpz_mod_evaluate(block_polynomial[i, j], points)
        end
        for k = 1:blocks
            evaluated_blocks[k][i, j] = ring(values[k])
        end
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
    raw_values = zeros(Int, length(points))
    GC.@preserve poly points raw_values begin
        @ccall Oscar.Nemo.libflint.fmpz_mod_poly_evaluate_fmpz_vec_fast(
            raw_values::Ptr{Int},
            poly::Ref{Oscar.Nemo.ZZModPolyRingElem},
            points::Ptr{Int},
            length(points)::Int,
            base_ring(parent(poly)).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct},
        )::Nothing
    end

    values = Vector{ZZRingElem}(undef, length(points))
    GC.@preserve raw_values values begin
        for i in eachindex(raw_values)
            value = ZZ(0)
            source = Ptr{ZZRingElem}(pointer(raw_values, i))
            @ccall Oscar.Nemo.libflint.fmpz_set(
                value::Ref{ZZRingElem},
                source::Ptr{ZZRingElem},
            )::Nothing
            values[i] = value
        end
        for i in eachindex(raw_values)
            source = Ptr{ZZRingElem}(pointer(raw_values, i))
            @ccall Oscar.Nemo.libflint.fmpz_clear(source::Ptr{ZZRingElem})::Nothing
        end
    end
    return values
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

function prodeval_linear!(a, b, start, stop, g, temp, g_temp, params::ZetaFunctionParams)
    interval_length = stop - start + 1
    if params.interval_product_algorithm === :sqrt &&
       interval_length >= params.interval_product_min_length &&
       _supports_sqrt_interval_product(a, b, g, temp, g_temp)
        return sqrt_interval_product!(a, b, start, stop, g, temp, g_temp)
    end
    return finitediff_prodeval_linear!(a, b, start, stop, g, temp, g_temp)
end
