struct _AKRReductionMap{T<:Integer}
    cofactor_length::Int
    reduced_length::Int
    offsets::Vector{Int}
    positions::Vector{Int}
    factors::Vector{T}
    derivative_positions::Vector{Int}
    derivative_factors::Vector{T}
end

struct _AKRRelationTable{N}
    cutoff::Int
    modulus::BigInt
    relations::Dict{NTuple{N,Int},Vector{Tuple{Int,NTuple{N,Int},BigInt}}}
    reduction_maps::Dict{Tuple{Int,Symbol},_AKRReductionMap{UInt64}}
    bigint_reduction_maps::Dict{Tuple{Int,Symbol},_AKRReductionMap{BigInt}}
end

_akr_bigint(value) = parse(BigInt, string(value))

function _akr_relation_axpy!(row, multiple, pivot, modulus)
    for (j, value) in pivot
        a = mod(get(row, j, big(0)) - multiple * value, modulus)
        if iszero(a)
            delete!(row, j)
        else
            row[j] = a
        end
    end
end

function _build_akr_relation_table(f, p, M, ::Val{N}) where {N}
    d = total_degree(f)
    cutoff = N * (d - 2) + 1
    modulus = big(p)^M
    monomials = NTuple{N,Int}[Tuple(e) for e in gen_exp_vec(N, cutoff, :lex)]
    monomial_indices = Dict(u => i for (i, u) in enumerate(monomials))
    cofactor_monomials =
        NTuple{N,Int}[Tuple(e) for e in gen_exp_vec(N, cutoff - d + 1, :lex)]
    generators = [(i, g) for i = 1:N for g in cofactor_monomials]
    partials = [Dict{NTuple{N,Int},BigInt}() for _ = 1:N]
    for i = 1:N, monomial in terms(derivative(f, i))
        exponent = Tuple(Int.(exponent_vector(monomial, 1)))
        coefficient = mod(_akr_bigint(lift(ZZ, coeff(monomial, 1))), modulus)
        iszero(coefficient) || (partials[i][exponent] = coefficient)
    end

    rows = Dict{Int,BigInt}[]
    for (j, (i, g)) in enumerate(generators)
        row = Dict(-j => big(1))
        for (exponent, coefficient) in partials[i]
            product = ntuple(l -> g[l] + exponent[l], Val(N))
            row[monomial_indices[product]] = coefficient
        end
        push!(rows, row)
    end

    pivots = zeros(Int, length(monomials))
    rank = 0
    for j in eachindex(monomials)
        row_index = findnext(row -> !iszero(mod(get(row, j, big(0)), p)), rows, rank + 1)
        isnothing(row_index) && continue
        rank += 1
        rows[rank], rows[row_index] = rows[row_index], rows[rank]
        pivot = rows[rank]
        scale = invmod(pivot[j], modulus)
        for l in keys(pivot)
            pivot[l] = mod(pivot[l] * scale, modulus)
        end
        for (i, row) in enumerate(rows)
            if i != rank && haskey(row, j)
                _akr_relation_axpy!(row, row[j], pivot, modulus)
            end
        end
        pivots[j] = rank
    end

    rank == length(monomials) || throw(
        ArgumentError(
            "the Jacobian ideal does not reach every degree-$cutoff monomial; the input may be singular or $p may be bad",
        ),
    )

    relations = Dict{NTuple{N,Int},Vector{Tuple{Int,NTuple{N,Int},BigInt}}}()
    for (j, monomial) in enumerate(monomials)
        cofactors = Tuple{Int,NTuple{N,Int},BigInt}[]
        for (column, coefficient) in rows[pivots[j]]
            if column < 0
                i, exponent = generators[-column]
                push!(cofactors, (i, exponent, coefficient))
            end
        end
        relations[monomial] = cofactors
    end
    return _AKRRelationTable{N}(
        cutoff,
        modulus,
        relations,
        Dict{Tuple{Int,Symbol},_AKRReductionMap{UInt64}}(),
        Dict{Tuple{Int,Symbol},_AKRReductionMap{BigInt}}(),
    )
end

_build_akr_relation_table(f, p, M) =
    _build_akr_relation_table(f, p, M, Val(nvars(parent(f))))

function _build_akr_reduction_map(
    table::_AKRRelationTable{N},
    m,
    f,
    termorder,
    ::Type{T},
) where {N,T<:Integer}
    d = total_degree(f)
    degree = d * m - N
    degree >= table.cutoff || error("AKR pole order is below the Jacobian cutoff")
    cofactor_degree = degree - d + 1
    source_monomials = gen_exp_vec(N, degree, termorder)
    cofactor_monomials = gen_exp_vec(N, cofactor_degree, termorder)
    cofactor_indices = Dict(Tuple(e) => i for (i, e) in enumerate(cofactor_monomials))
    reduced_monomials = gen_exp_vec(N, cofactor_degree - 1, termorder)
    reduced_indices = Dict(Tuple(e) => i for (i, e) in enumerate(reduced_monomials))
    cofactor_length = length(cofactor_monomials)
    offsets = Vector{Int}(undef, length(source_monomials) + 1)
    positions = Int[]
    factors = T[]
    derivative_positions = zeros(Int, N * cofactor_length)
    derivative_factors = zeros(T, N * cofactor_length)

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
            push!(factors, T(mod(factor, table.modulus)))
        end
    end
    offsets[end] = length(positions) + 1

    for i = 1:N, (j, exponent) in enumerate(cofactor_monomials)
        exponent[i] == 0 && continue
        reduced_exponent = ntuple(l -> exponent[l] - (l == i), Val(N))
        position = (i - 1) * cofactor_length + j
        derivative_positions[position] = reduced_indices[reduced_exponent]
        derivative_factors[position] = T(exponent[i])
    end

    return _AKRReductionMap{T}(
        cofactor_length,
        length(reduced_monomials),
        offsets,
        positions,
        factors,
        derivative_positions,
        derivative_factors,
    )
end

function _akr_reduction_map(table::_AKRRelationTable, m, f, termorder)
    return get!(table.reduction_maps, (m, termorder)) do
        _build_akr_reduction_map(table, m, f, termorder, UInt64)
    end
end

function _akr_bigint_reduction_map(table::_AKRRelationTable, m, f, termorder)
    return get!(table.bigint_reduction_maps, (m, termorder)) do
        _build_akr_reduction_map(table, m, f, termorder, BigInt)
    end
end

_akr_residue_uint64(value::UInt64) = value
_akr_residue_uint64(value) = UInt64(lift(ZZ, value))

function _akr_cofactor_residues(
    g,
    reduction_map::_AKRReductionMap{UInt64},
    modulus::UInt64,
)
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

function _akr_cofactor_vector(g, m, f, table::_AKRRelationTable{N}, params) where {N}
    table.modulus > typemax(UInt64) &&
        return _akr_cofactor_vector_bigint(g, m, f, table, params)
    reduction_map = _akr_reduction_map(table, m, f, params.termorder)
    residues = _akr_cofactor_residues(g, reduction_map, UInt64(table.modulus))
    return _akr_residues_to_ring(residues, f)
end

function _akr_cofactor_vector_bigint(g, m, f, table::_AKRRelationTable{N}, params) where {N}
    reduction_map = _akr_bigint_reduction_map(table, m, f, params.termorder)
    R = coefficient_ring(parent(f))
    cofactor_vector = fill(R(0), N * reduction_map.cofactor_length)

    for j = 1:(length(reduction_map.offsets)-1)
        iszero(g[j]) && continue
        for k in reduction_map.offsets[j]:(reduction_map.offsets[j+1]-1)
            position = reduction_map.positions[k]
            cofactor_vector[position] += g[j] * R(reduction_map.factors[k])
        end
    end
    return cofactor_vector
end

_akr_cofactor_vector(g, m, f, picache::AbstractVector, params) = picache[m] * g

function _reducechain_akr_naive(g, m, f, p, picache, params)
    if picache isa _AKRRelationTable && picache.modulus <= typemax(UInt64)
        n = nvars(parent(f)) - 1
        current = UInt64[_akr_residue_uint64(coefficient) for coefficient in g]
        modulus = UInt64(picache.modulus)
        while m > n
            reduction_map = _akr_reduction_map(picache, m, f, params.termorder)
            cofactors = _akr_cofactor_residues(current, reduction_map, modulus)
            current = _akr_derivative_residues(cofactors, reduction_map, modulus)
            m -= 1
        end
        return _akr_residues_to_ring(current, f)
    end

    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))

    while m > n
        gtemp = _akr_cofactor_vector(g, m, f, picache, params)
        len = Int(length(gtemp)/(n+1))
        g = derivative(
            vector_to_polynomial(gtemp[1:len], n, d*(m-1)-n, PR, params.termorder),
            1,
        )
        for i = 1:n
            gpolytemp = derivative(
                vector_to_polynomial(
                    gtemp[(i*len+1):((i+1)*len)],
                    n,
                    d*(m-1)-n,
                    PR,
                    params.termorder,
                ),
                i+1,
            )
            g = g + gpolytemp
        end
        g = polynomial_to_vector(g, n+1, params.termorder)
        if size(g)[1] == 0
            g = fill(R(0), binomial(d*(m-1)-1, n))
        end
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
            reduced = _reducechain_akr_naive(g, t[2], f, p, picache, params)
            reduced_poly = vector_to_polynomial(reduced, n, d*n-n-1, PR, params.termorder)
            result += reduced_poly
        end
    end

    return [[result, n]]
end

function _build_akr_picache(highm, n, d, S, f, pseudoInverseMat, cache, params, p)
    M = factor(modulus(base_ring(parent(f))))[p]
    if sort(S) == collect(0:n) && mod(d, p) != 0
        return _build_akr_relation_table(f, p, M)
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
function _reduce_akr_level(g, m, f, picache, params)
    if picache isa _AKRRelationTable && picache.modulus <= typemax(UInt64)
        modulus = UInt64(picache.modulus)
        reduction_map = _akr_reduction_map(picache, m, f, params.termorder)
        cofactors = _akr_cofactor_residues(g, reduction_map, modulus)
        reduced = _akr_derivative_residues(cofactors, reduction_map, modulus)
        return _akr_residues_to_ring(reduced, f)
    end

    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    gtemp = _akr_cofactor_vector(g, m, f, picache, params)
    len = Int(length(gtemp)/(n+1))
    reduced = derivative(
        vector_to_polynomial(gtemp[1:len], n, d*(m-1)-n, PR, params.termorder),
        1,
    )
    for i = 1:n
        component = vector_to_polynomial(
            gtemp[(i*len+1):((i+1)*len)],
            n,
            d*(m-1)-n,
            PR,
            params.termorder,
        )
        reduced += derivative(component, i+1)
    end
    reduced_vector = polynomial_to_vector(reduced, n+1, params.termorder)
    if size(reduced_vector)[1] == 0
        reduced_vector = fill(R(0), binomial(d*(m-1)-1, n))
    end
    return reduced_vector
end

"""
    _reducepoly_akr_levelbylevel(pol, S, f, p, picache, params)

Reduce all terms at each pole order together. The derivative still uses the
existing polynomial-vector conversions used by the naive AKR implementation.
"""
function _reducepoly_akr_levelbylevel(pol, S, f, p, picache, params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    isempty(pol) && return [[PR(), n]]
    highm = maximum(Int(term[2]) for term in pol)
    degree = d*highm - n - 1
    g = fill(R(0), binomial(degree+n, n))

    for m = highm:-1:(n+1)
        for term in termsoforder(pol, m)
            g .+= polynomial_to_vector(term[1], n+1, params.termorder)
        end
        g = _reduce_akr_level(g, m, f, picache, params)
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
