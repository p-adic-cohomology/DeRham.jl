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

function _build_akr_reduction_map(
    table::_AKRRelationTable{T,N},
    m,
    f,
    termorder,
) where {T,N}
    d = total_degree(f)
    degree = d * m - N
    degree >= table.cutoff || error("AKR pole order is below the Jacobian cutoff")
    cofactor_degree = degree - d + 1
    source_monomials = gen_exp_vec(N, degree, termorder)
    cofactor_monomials = gen_exp_vec(N, cofactor_degree, termorder)
    cofactor_indices = Dict(Tuple(e) => i for (i, e) in enumerate(cofactor_monomials))
    cofactor_length = length(cofactor_monomials)
    offsets = Vector{Int}(undef, length(source_monomials) + 1)
    positions = Int[]
    factors = T[]

    for (j, exponent) in enumerate(source_monomials)
        offsets[j] = length(positions) + 1
        left = table.cutoff
        h = ntuple(Val(N)) do i
            part = min(left, exponent[i])
            left -= part
            part
        end
        left == 0 || error("failed to split a monomial at the AKR Jacobian cutoff")
        shift = ntuple(i -> exponent[i] - h[i], Val(N))
        relation = table.relations[h]
        for (i, cofactor, factor) in relation
            shifted = ntuple(l -> shift[l] + cofactor[l], Val(N))
            index = get(cofactor_indices, shifted, 0)
            index != 0 || error("AKR cofactor has the wrong degree")
            push!(positions, (i - 1) * cofactor_length + index)
            push!(factors, convert(T, factor))
        end
    end
    offsets[end] = length(positions) + 1

    derivative_positions = Int[]
    derivative_factors = T[]
    reduced_length = 0
    if T === UInt64
        reduced_monomials = gen_exp_vec(N, cofactor_degree - 1, termorder)
        reduced_indices = Dict(Tuple(e) => i for (i, e) in enumerate(reduced_monomials))
        reduced_length = length(reduced_monomials)
        derivative_positions = zeros(Int, N * cofactor_length)
        derivative_factors = zeros(T, N * cofactor_length)
        for i = 1:N, (j, exponent) in enumerate(cofactor_monomials)
            exponent[i] == 0 && continue
            reduced_exponent = ntuple(l -> exponent[l] - (l == i), Val(N))
            position = (i - 1) * cofactor_length + j
            derivative_positions[position] = reduced_indices[reduced_exponent]
            derivative_factors[position] = T(exponent[i])
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
Reuse materialized pole-order maps while their combined array payload fits
the cache budget. Larger maps are used for one reduction step. BigInt factors
share the lifted coefficients in the cutoff relation table.
"""
function _akr_reduction_map(table::_AKRRelationTable, m, f, termorder)
    key = (m, termorder)
    cached = get(table.reduction_maps, key, nothing)
    cached === nothing || return cached

    reduction_map = _build_akr_reduction_map(table, m, f, termorder)
    map_bytes = sum(
        sizeof,
        (
            reduction_map.offsets,
            reduction_map.positions,
            reduction_map.factors,
            reduction_map.derivative_positions,
            reduction_map.derivative_factors,
        ),
    )
    if map_bytes <= _AKR_REDUCTION_MAP_CACHE_BYTES - table.cached_map_bytes
        table.reduction_maps[key] = reduction_map
        table.cached_map_bytes += map_bytes
    end
    return reduction_map
end

_akr_residue_uint64(value::UInt64) = value
_akr_residue_uint64(value) = UInt64(lift(ZZ, value))

function _akr_cofactor_residues(g, reduction_map::_AKRReductionMap{UInt64}, modulus::UInt64)
    accumulators = zeros(UInt64, length(reduction_map.derivative_positions))
    for j in eachindex(g)
        coefficient = _akr_residue_uint64(g[j])
        coefficient == 0 && continue
        for k = reduction_map.offsets[j]:(reduction_map.offsets[j+1]-1)
            position = reduction_map.positions[k]
            accumulators[position] = UInt64(
                rem(
                    UInt128(accumulators[position]) +
                    UInt128(coefficient) * UInt128(reduction_map.factors[k]),
                    UInt128(modulus),
                ),
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
        reduced[position] = UInt64(
            rem(
                UInt128(reduced[position]) +
                UInt128(coefficient) * UInt128(reduction_map.derivative_factors[i]),
                UInt128(modulus),
            ),
        )
    end
    return reduced
end

function _akr_reduce_level_residues(g, m, f, table::_AKRRelationTable, termorder, modulus)
    reduction_map = _akr_reduction_map(table, m, f, termorder)
    cofactors = _akr_cofactor_residues(g, reduction_map, modulus)
    return _akr_derivative_residues(cofactors, reduction_map, modulus)
end

function _akr_add_residues!(accumulator::Vector{UInt64}, values, modulus::UInt64)
    length(accumulator) == length(values) ||
        error("AKR residue vector has the wrong length")
    for i in eachindex(accumulator)
        value = _akr_residue_uint64(values[i])
        value == 0 && continue
        accumulator[i] =
            UInt64(rem(UInt128(accumulator[i]) + UInt128(value), UInt128(modulus)))
    end
    return accumulator
end

function _akr_cofactor_vector(g, m, f, table::_AKRRelationTable{UInt64}, params)
    reduction_map = _akr_reduction_map(table, m, f, params.termorder)
    residues = _akr_cofactor_residues(g, reduction_map, UInt64(table.modulus))
    return _akr_residues_to_ring(residues, f)
end

function _akr_cofactor_vector(g, m, f, table::_AKRRelationTable{BigInt,N}, params) where {N}
    reduction_map = _akr_reduction_map(table, m, f, params.termorder)
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
    in_degree = l - d + 1
    out_degree = l - d
    in_length = binomial(n + in_degree, n)
    out_length = binomial(n + out_degree, n)
    vars = gens(PR)
    reduced = PR()
    offset = 1

    for i = 0:n
        if i in S
            block = gtemp[offset:(offset+in_length-1)]
            component = vector_to_polynomial(block, n, in_degree, PR, params.termorder)
            offset += in_length
        else
            block = gtemp[offset:(offset+out_length-1)]
            component =
                vars[i+1] * vector_to_polynomial(block, n, out_degree, PR, params.termorder)
            offset += out_length
        end
        reduced += derivative(component, i + 1)
    end

    offset == length(gtemp) + 1 || error("AKR cofactor blocks have the wrong length")
    reduced_vector = polynomial_to_vector(reduced, n + 1, params.termorder)
    if size(reduced_vector)[1] == 0
        reduced_vector = fill(R(0), binomial(d * (m - 1) - 1, n))
    end
    return reduced_vector
end

function _reducechain_akr_naive(g, m, f, p, picache, S, params)
    if picache isa _AKRRelationTable{UInt64}
        n = nvars(parent(f)) - 1
        current = UInt64[_akr_residue_uint64(coefficient) for coefficient in g]
        modulus = UInt64(picache.modulus)
        while m > n
            current = _akr_reduce_level_residues(current, m, f, picache, params.termorder, modulus)
            m -= 1
        end
        return _akr_residues_to_ring(current, f)
    end

    n = nvars(parent(f)) - 1
    while m > n
        gtemp = _akr_cofactor_vector(g, m, f, picache, params)
        g = _akr_pseudoinverse_derivative(gtemp, m, f, S, params)
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

_akr_cutoff_eligible(n, d, S, p) = sort(S) == collect(0:n) && mod(d, p) != 0

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
        return _build_akr_relation_table(f, p, M; always_use_bigints = params.always_use_bigints)
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
    _reduce_akr_level(g, m, f, picache, params)

Apply one AKR reduction step to the combined coefficient vector at pole order
`m`, returning the coefficient vector for pole order `m - 1`.
"""
function _reduce_akr_level(g, m, f, picache, params, S = collect(0:(nvars(parent(f))-1)))
    if picache isa _AKRRelationTable{UInt64}
        modulus = UInt64(picache.modulus)
        residues = UInt64[_akr_residue_uint64(coefficient) for coefficient in g]
        reduced =
            _akr_reduce_level_residues(residues, m, f, picache, params.termorder, modulus)
        return _akr_residues_to_ring(reduced, f)
    end

    gtemp = _akr_cofactor_vector(g, m, f, picache, params)
    return _akr_pseudoinverse_derivative(gtemp, m, f, S, params)
end

"""
    _reducepoly_akr_levelbylevel(pol, S, f, p, picache, params)

Reduce all terms at each pole order together. Both schedules use coefficient
maps for derivatives on the UInt64 cutoff path. The BigInt and pseudo-inverse
paths differentiate polynomials.
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
