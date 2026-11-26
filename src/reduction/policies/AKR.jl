struct _AKRReductionMap{T<:Integer}
    cofactor_length::Int
    reduced_length::Int
    offsets::Vector{Int}
    positions::Vector{Int}
    factors::Vector{T}
    derivative_positions::Vector{Int}
    derivative_factors::Vector{T}
end

const _AKR_REDUCTION_MAP_CACHE_BYTES = 64 * 1024^2

mutable struct _AKRRelationTable{T<:Integer,N}
    cutoff::Int
    modulus::BigInt
    relations::Dict{NTuple{N,Int},Vector{Tuple{Int,NTuple{N,Int},BigInt}}}
    reduction_maps::Dict{Tuple{Int,Symbol},_AKRReductionMap{T}}
    rejected_maps::Set{Tuple{Int,Symbol}}
    cached_map_bytes::Int
end

function _AKRRelationTable(
    inverse::_SparseJacobianRightInverse{N},
    always_use_bigints,
) where {N}
    T = always_use_bigints || inverse.modulus > typemax(UInt64) ? BigInt : UInt64
    return _AKRRelationTable{T,N}(
        inverse.cutoff,
        inverse.modulus,
        inverse.relations,
        Dict{Tuple{Int,Symbol},_AKRReductionMap{T}}(),
        Set{Tuple{Int,Symbol}}(),
        0,
    )
end

function _build_akr_relation_table(f, p, M; always_use_bigints = false)
    inverse = _build_sparse_jacobian_right_inverse(f, p, M)
    if inverse === nothing
        cutoff = nvars(parent(f)) * (total_degree(f) - 2) + 1
        throw(
            ArgumentError(
                "the Jacobian ideal does not reach every degree-$cutoff monomial; the input may be singular or $p may be bad",
            ),
        )
    end
    return _AKRRelationTable(inverse, always_use_bigints)
end

function _akr_cutoff_split(exponent, cutoff, ::Val{N}) where {N}
    left = cutoff
    h = ntuple(Val(N)) do i
        part = min(left, exponent[i])
        left -= part
        part
    end
    left == 0 || error("failed to split a monomial at the AKR Jacobian cutoff")
    shift = ntuple(i -> exponent[i] - h[i], Val(N))
    return h, shift
end

function _build_akr_reduction_map(
    table::_AKRRelationTable{T,N},
    m,
    f,
    termorder,
) where {T,N}
    ordering = _graded_monomial_order(termorder)
    d = total_degree(f)
    degree = d * m - N
    degree >= table.cutoff || error("AKR pole order is below the Jacobian cutoff")
    cofactor_degree = degree - d + 1
    source_length = binomial(degree + N - 1, N - 1)
    cofactor_length = binomial(cofactor_degree + N - 1, N - 1)
    offsets = Vector{Int}(undef, source_length + 1)
    positions = Int[]
    factors = T[]
    exponent = Vector{Int}(undef, N)

    for j = 1:source_length
        _graded_exponent_at!(exponent, degree, j, ordering)
        offsets[j] = length(positions) + 1
        h, shift = _akr_cutoff_split(exponent, table.cutoff, Val(N))
        for (i, cofactor, factor) in table.relations[h]
            shifted = ntuple(l -> shift[l] + cofactor[l], Val(N))
            index = _graded_monomial_index(shifted, cofactor_degree, ordering)
            push!(positions, (i - 1) * cofactor_length + index)
            push!(factors, convert(T, factor))
        end
    end
    offsets[end] = length(positions) + 1

    derivative_positions = Int[]
    derivative_factors = T[]
    reduced_length = 0
    if T === UInt64
        reduced_degree = cofactor_degree - 1
        reduced_length = binomial(reduced_degree + N - 1, N - 1)
        derivative_positions = zeros(Int, N * cofactor_length)
        derivative_factors = zeros(T, N * cofactor_length)
        for j = 1:cofactor_length
            _graded_exponent_at!(exponent, cofactor_degree, j, ordering)
            for i = 1:N
                exponent[i] == 0 && continue
                reduced_exponent = ntuple(l -> exponent[l] - (l == i), Val(N))
                position = (i - 1) * cofactor_length + j
                derivative_positions[position] =
                    _graded_monomial_index(reduced_exponent, reduced_degree, ordering)
                derivative_factors[position] = T(exponent[i])
            end
        end
    end

    return _AKRReductionMap{T}(
        cofactor_length,
        reduced_length,
        offsets,
        positions,
        factors,
        derivative_positions,
        derivative_factors,
    )
end

"""
Compute the payload size of a pole-order map within the available bytes.
The greedy cutoff split fixes every exponent before the last nonzero
coordinate of the cutoff monomial. Additional degree can go in that
coordinate or later.
"""
function _akr_reduction_map_bytes(
    table::_AKRRelationTable{T,N},
    degree,
    d,
    available,
) where {T,N}
    source_length = binomial(degree + N - 1, N - 1)
    source_length < div(available, sizeof(Int)) || return nothing
    bytes = (source_length + 1) * sizeof(Int)
    if T === UInt64
        cofactor_length = binomial(degree - d + N, N - 1)
        bytes_per_monomial = N * (sizeof(Int) + sizeof(T))
        cofactor_length <= div(available - bytes, bytes_per_monomial) || return nothing
        bytes += cofactor_length * bytes_per_monomial
    end

    bytes_per_entry = sizeof(Int) + Base.elsize(Vector{T})
    for (h, relation) in table.relations
        isempty(relation) && continue
        trailing = N - something(findlast(!iszero, h), 1)
        extensions = binomial(degree - table.cutoff + trailing, trailing)
        entry_count = length(relation)
        extensions <= div(div(available - bytes, bytes_per_entry), entry_count) ||
            return nothing
        bytes += extensions * entry_count * bytes_per_entry
    end
    return bytes
end

"""
Reuse pole-order maps whose combined array payload fits the cache budget.
Remember maps that exceed the remaining budget, which decreases as maps are added.
BigInt factors share the lifted coefficients in the cutoff relation table.
"""
function _akr_reduction_map(table::_AKRRelationTable{T,N}, m, f, termorder) where {T,N}
    key = (m, termorder)
    cached = get(table.reduction_maps, key, nothing)
    cached === nothing || return cached
    key in table.rejected_maps && return nothing

    d = total_degree(f)
    degree = d * m - N
    degree >= table.cutoff || error("AKR pole order is below the Jacobian cutoff")
    available = _AKR_REDUCTION_MAP_CACHE_BYTES - table.cached_map_bytes
    map_bytes = _akr_reduction_map_bytes(table, degree, d, available)
    if map_bytes === nothing
        push!(table.rejected_maps, key)
        return nothing
    end
    reduction_map = _build_akr_reduction_map(table, m, f, termorder)
    table.reduction_maps[key] = reduction_map
    table.cached_map_bytes += map_bytes
    return reduction_map
end

_akr_residue_uint64(value::UInt64) = value
_akr_residue_uint64(value) = UInt64(lift(ZZ, value))

_akr_residue(value, ::Type{UInt64}) = _akr_residue_uint64(value)
_akr_residue(value, ::Type{BigInt}) = BigInt(lift(ZZ, value))

function _akr_initial_residues(g)
    residues = zeros(UInt64, length(g))
    for i in eachindex(g)
        coefficient = g[i]
        iszero(coefficient) && continue
        residues[i] = _akr_residue_uint64(coefficient)
    end
    return residues
end

_akr_add_product(a::UInt64, b::UInt64, c::UInt64, modulus::UInt64) =
    UInt64(rem(UInt128(a) + UInt128(b) * UInt128(c), UInt128(modulus)))
_akr_add_product(a::BigInt, b::BigInt, c::BigInt, modulus::BigInt) = mod(a + b * c, modulus)

_akr_scale_factor(factor, power, modulus::UInt64) =
    UInt64(rem(UInt128(UInt64(factor)) * UInt128(power), UInt128(modulus)))
_akr_scale_factor(factor, power, ::BigInt) = factor * power

function _akr_apply_uncached(g, m, f, table::_AKRRelationTable{T,N}, termorder) where {T,N}
    ordering = _graded_monomial_order(termorder)
    d = total_degree(f)
    degree = d * m - N
    output_degree = degree - d
    result = zeros(T, binomial(output_degree + N - 1, N - 1))
    modulus = convert(T, table.modulus)
    exponent = Vector{Int}(undef, N)

    for j in eachindex(g)
        iszero(g[j]) && continue
        coefficient = _akr_residue(g[j], T)
        _graded_exponent_at!(exponent, degree, j, ordering)
        h, shift = _akr_cutoff_split(exponent, table.cutoff, Val(N))
        for (i, cofactor, factor) in table.relations[h]
            shifted = ntuple(k -> shift[k] + cofactor[k], Val(N))
            power = shifted[i]
            power == 0 && continue
            target = ntuple(k -> shifted[k] - (k == i), Val(N))
            position = _graded_monomial_index(target, output_degree, ordering)
            scaled_factor = _akr_scale_factor(factor, power, modulus)
            result[position] =
                _akr_add_product(result[position], coefficient, scaled_factor, modulus)
        end
    end
    return result
end

function _akr_cofactor_residues(g, reduction_map::_AKRReductionMap{UInt64}, modulus::UInt64)
    accumulators = zeros(UInt64, length(reduction_map.derivative_positions))
    for j in eachindex(g)
        value = g[j]
        iszero(value) && continue
        coefficient = _akr_residue_uint64(value)
        for k = reduction_map.offsets[j]:(reduction_map.offsets[j+1]-1)
            position = reduction_map.positions[k]
            accumulators[position] = _akr_add_product(
                accumulators[position],
                coefficient,
                reduction_map.factors[k],
                modulus,
            )
        end
    end
    return accumulators
end

function _akr_residues_to_ring(values, f)
    R = coefficient_ring(parent(f))
    result = fill(R(0), length(values))
    for i in eachindex(values)
        result[i] = R(values[i])
    end
    return result
end

function _akr_derivative_residues(
    cofactor_vector,
    reduction_map::_AKRReductionMap{UInt64},
    modulus::UInt64,
)
    reduced = zeros(UInt64, reduction_map.reduced_length)
    for i in eachindex(cofactor_vector)
        position = reduction_map.derivative_positions[i]
        position == 0 && continue
        coefficient = cofactor_vector[i]
        coefficient == 0 && continue
        reduced[position] = _akr_add_product(
            reduced[position],
            coefficient,
            reduction_map.derivative_factors[i],
            modulus,
        )
    end
    return reduced
end

function _akr_has_at_most_one_term(g)
    first_nonzero = findfirst(!iszero, g)
    return first_nonzero === nothing || findnext(!iszero, g, first_nonzero + 1) === nothing
end

function _akr_reduce_level_residues(
    g,
    m,
    f,
    table::_AKRRelationTable{UInt64},
    termorder,
    modulus,
)
    if _akr_has_at_most_one_term(g)
        return _akr_apply_uncached(g, m, f, table, termorder)
    end
    reduction_map = _akr_reduction_map(table, m, f, termorder)
    reduction_map === nothing && return _akr_apply_uncached(g, m, f, table, termorder)
    cofactors = _akr_cofactor_residues(g, reduction_map, modulus)
    return _akr_derivative_residues(cofactors, reduction_map, modulus)
end

function _akr_add_residues!(accumulator::Vector{UInt64}, values, modulus::UInt64)
    length(accumulator) == length(values) ||
        error("AKR residue vector has the wrong length")
    for i in eachindex(accumulator)
        value = values[i]
        iszero(value) && continue
        residue = _akr_residue_uint64(value)
        accumulator[i] =
            UInt64(rem(UInt128(accumulator[i]) + UInt128(residue), UInt128(modulus)))
    end
    return accumulator
end

function _akr_cofactor_vector(
    g,
    f,
    table::_AKRRelationTable{BigInt,N},
    reduction_map,
) where {N}
    R = coefficient_ring(parent(f))
    cofactor_vector = fill(R(0), N * reduction_map.cofactor_length)

    for j = 1:(length(reduction_map.offsets)-1)
        iszero(g[j]) && continue
        for k = reduction_map.offsets[j]:(reduction_map.offsets[j+1]-1)
            position = reduction_map.positions[k]
            cofactor_vector[position] += g[j] * R(reduction_map.factors[k])
        end
    end
    return cofactor_vector
end

_akr_cofactor_vector(g, m, f, picache::AbstractVector, params) = picache[m] * g

function _akr_pseudoinverse_derivative(gtemp, m, f, S, params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    l = d * m - n - 1
    blocks = _controlled_cofactor_blocks(n, d, l, S)
    last(last(blocks).columns) == length(gtemp) ||
        error("AKR cofactor blocks have the wrong length")
    vars = gens(PR)
    reduced = PR()

    for block in blocks
        component = vector_to_polynomial(
            gtemp[block.columns],
            n,
            block.degree,
            PR,
            params.termorder,
        )
        if !block.selected
            component = vars[block.variable] * component
        end
        reduced += derivative(component, block.variable)
    end

    reduced_vector = polynomial_to_vector(reduced, n + 1, params.termorder)
    if size(reduced_vector)[1] == 0
        reduced_vector = fill(R(0), binomial(d * (m - 1) - 1, n))
    end
    return reduced_vector
end

function _reducechain_akr_naive(g, m, f, p, picache, S, params)
    if picache isa _AKRRelationTable{UInt64}
        n = nvars(parent(f)) - 1
        current = _akr_initial_residues(g)
        modulus = UInt64(picache.modulus)
        while m > n
            current = _akr_reduce_level_residues(
                current,
                m,
                f,
                picache,
                params.termorder,
                modulus,
            )
            m -= 1
        end
        return _akr_residues_to_ring(current, f)
    end

    n = nvars(parent(f)) - 1
    while m > n
        g = _reduce_akr_level(g, m, f, picache, params, S)
        m = m - 1
    end

    return g
end

function _reducepoly_akr_naive(pol, S, f, p, picache, params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    result = PR()
    for term in pol
        terms = termsoforder(pol, term[2])
        for t in terms
            g = polynomial_to_vector(t[1], n+1, params.termorder)
            reduced = _reducechain_akr_naive(g, t[2], f, p, picache, S, params)
            reduced_poly = vector_to_polynomial(reduced, n, d*n-n-1, PR, params.termorder)
            result += reduced_poly
        end
    end

    return [[result, n]]
end

_akr_cutoff_eligible(n, d, S, p) = length(S) == n + 1 && all(in(S), 0:n) && mod(d, p) != 0

function _find_reduction_model(
    ::Union{Val{:akr},Val{:akr_naive},Val{:akr_levelbylevel}},
    f,
    P,
    p,
    M,
    S,
    params,
    changef,
    cache,
    precision_poly,
)
    if !_akr_cutoff_eligible(P.n, P.d, S, p)
        return _find_reduction_model(
            Val(:controlled),
            f,
            P,
            p,
            M,
            S,
            params,
            changef,
            cache,
            precision_poly,
        )
    end
    lifted = liftCoefficients(base_ring(precision_poly), precision_poly, f)
    right_inverse = _build_sparse_jacobian_right_inverse(lifted, p, M)
    if right_inverse === nothing
        changef && return false, false, nothing
        throw(ArgumentError("f is not smooth"))
    end
    return false, f, (; lifted, right_inverse)
end

function _build_akr_picache(
    highm,
    n,
    d,
    S,
    f,
    inverse::_SparseJacobianRightInverse,
    cache,
    params,
    p,
)
    return _AKRRelationTable(inverse, params.always_use_bigints)
end

function _build_akr_picache(highm, n, d, S, f, pseudoInverseMat, cache, params, p)
    M = factor(modulus(base_ring(parent(f))))[p]
    if _akr_cutoff_eligible(n, d, S, p)
        return _build_akr_relation_table(
            f,
            p,
            M;
            always_use_bigints = params.always_use_bigints,
        )
    end
    picache = Vector{typeof(pseudoInverseMat)}(undef, highm)
    picache[n] = pseudoInverseMat
    for i = (n+1):highm
        l = d*i - n - 1
        generate_degree_forward(cache, l)
        generate_degree_forward(cache, l-(d-1))
        generate_degree_forward(cache, l-d)
        pi_new = pseudo_inverse_controlled_lifted(f, S, l, M, params, cache)
        MS = matrix_space(base_ring(parent(f)), nrows(pi_new), ncols(pi_new))
        pseudo_inverse_mat = MS()
        for i = 1:nrows(pi_new)
            for j = 1:ncols(pi_new)
                pseudo_inverse_mat[i, j] = ZZ(pi_new[i, j])
            end
        end
        picache[i] = pseudo_inverse_mat
    end
    return picache
end

function _highest_akr_pole_order(FT)
    return FT[length(FT)][length(FT[length(FT)])][2]
end

function reducetransform_akr_naive(FT, N_m, S, f, pseudoInverseMat, p, cache, params)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    highm = _highest_akr_pole_order(FT)
    if 0 < params.verbose
        println("Building AKR relation table/cache")
        @time picache =
            _build_akr_picache(highm, n, d, S, f, pseudoInverseMat, cache, params, p)
    else
        picache = _build_akr_picache(highm, n, d, S, f, pseudoInverseMat, cache, params, p)
    end

    result = similar(FT)
    reduction_loop_start = time_ns()
    for i = 1:length(FT)
        pol = FT[i]
        if 0 < params.verbose
            println("Reducing vector $i")
            @time reduction = _reducepoly_akr_naive(pol, S, f, p, picache, params)
        else
            reduction = _reducepoly_akr_naive(pol, S, f, p, picache, params)
        end
        result[i] = reduction
    end
    0 < params.verbose &&
        println("AKR reduction loop: ", (time_ns() - reduction_loop_start) / 1e9, " s")

    return result
end

"""
    _reduce_akr_level(g, m, f, picache, params, S = 0:(nvars(parent(f)) - 1))

Apply one AKR reduction step to the combined coefficient vector at pole order
`m`, returning the coefficient vector for pole order `m - 1`.
The pseudo-inverse path uses `S` to recover the cofactor blocks. Pass the same
zero-based selection that built the inverse when it uses a proper subset.
"""
function _reduce_akr_level(g, m, f, picache, params, S = 0:(nvars(parent(f))-1))
    if picache isa _AKRRelationTable{UInt64}
        modulus = UInt64(picache.modulus)
        residues = _akr_initial_residues(g)
        reduced =
            _akr_reduce_level_residues(residues, m, f, picache, params.termorder, modulus)
        return _akr_residues_to_ring(reduced, f)
    end

    gtemp = if picache isa _AKRRelationTable{BigInt}
        reduction_map =
            _akr_has_at_most_one_term(g) ? nothing :
            _akr_reduction_map(picache, m, f, params.termorder)
        if reduction_map === nothing
            reduced = _akr_apply_uncached(g, m, f, picache, params.termorder)
            return _akr_residues_to_ring(reduced, f)
        end
        _akr_cofactor_vector(g, f, picache, reduction_map)
    else
        _akr_cofactor_vector(g, m, f, picache, params)
    end
    return _akr_pseudoinverse_derivative(gtemp, m, f, S, params)
end

"""
    _reducepoly_akr_levelbylevel(pol, S, f, p, picache, params)

Reduce all terms at each pole order together and add the input numerators
already at the terminal order. Both schedules apply cutoff relations directly
to sparse inputs or when there is insufficient room for a pole-order map.
The pseudo-inverse path differentiates polynomials.
"""
function _reducepoly_akr_levelbylevel(pol, S, f, p, picache, params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    isempty(pol) && return [[PR(), n]]
    highm = maximum(Int(term[2]) for term in pol)
    degree = d*highm - n - 1

    if picache isa _AKRRelationTable{UInt64}
        modulus = UInt64(picache.modulus)
        g = zeros(UInt64, binomial(degree + n, n))

        for m = highm:-1:(n+1)
            for term in termsoforder(pol, m)
                _akr_add_residues!(
                    g,
                    polynomial_to_vector(term[1], n + 1, params.termorder),
                    modulus,
                )
            end
            g = _akr_reduce_level_residues(g, m, f, picache, params.termorder, modulus)
        end

        g = _akr_residues_to_ring(g, f)
    else
        g = fill(R(0), binomial(degree + n, n))

        for m = highm:-1:(n+1)
            for term in termsoforder(pol, m)
                g .+= polynomial_to_vector(term[1], n + 1, params.termorder)
            end
            g = _reduce_akr_level(g, m, f, picache, params, S)
        end
    end

    result = vector_to_polynomial(g, n, d*n-n-1, PR, params.termorder)
    for (numerator, m) in pol
        m == n && (result += numerator)
    end
    return [[result, n]]
end

function reducetransform_akr_levelbylevel(FT, N_m, S, f, pseudoInverseMat, p, cache, params)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    highm = _highest_akr_pole_order(FT)
    if 0 < params.verbose
        println("Building AKR relation table/cache")
        @time picache =
            _build_akr_picache(highm, n, d, S, f, pseudoInverseMat, cache, params, p)
    else
        picache = _build_akr_picache(highm, n, d, S, f, pseudoInverseMat, cache, params, p)
    end

    result = similar(FT)
    reduction_loop_start = time_ns()
    for i = 1:length(FT)
        pol = FT[i]
        if 0 < params.verbose
            println("Reducing vector $i")
            @time reduction = _reducepoly_akr_levelbylevel(pol, S, f, p, picache, params)
        else
            reduction = _reducepoly_akr_levelbylevel(pol, S, f, p, picache, params)
        end
        result[i] = reduction
    end
    0 < params.verbose &&
        println("AKR reduction loop: ", (time_ns() - reduction_loop_start) / 1e9, " s")

    return result
end
