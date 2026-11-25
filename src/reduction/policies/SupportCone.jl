"""
    SupportConePlanner(support)

Build exact support directions and a bounded, deterministic set of nonsingular
support subcones. Directions are exponent tuples sorted lexicographically. The
lexicographically first full-rank subcone is always included; at most the first
256 combinations are searched and at most 64 subcones are retained. Subcones
may overlap; they are not a triangulation. Determinants and Cramer coordinates
use integer arithmetic.
"""
struct SupportConeSubcone{N}
    indices::NTuple{N,Int}
    determinant::BigInt
    adjugate::Matrix{BigInt}
end

struct SupportConePlanner{N}
    support::Vector{NTuple{N,Int}}
    subcones::Vector{SupportConeSubcone{N}}
end

struct SupportConeChunk{N}
    direction::NTuple{N,Int}
    steps::Int
end

struct SupportConePlan{N}
    subcone::NTuple{N,Int}
    chunks::Vector{SupportConeChunk{N}}
end

const _SUPPORTCONE_MAX_SUBCONE_COMBINATIONS = 256
const _SUPPORTCONE_MAX_SUBCONES = 64

function _supportcone_determinant(matrix::AbstractMatrix{<:Integer})
    size(matrix, 1) == size(matrix, 2) ||
        throw(ArgumentError("a support-cone determinant requires a square matrix"))
    dimension = size(matrix, 1)
    dimension == 0 && return BigInt(1)

    work = BigInt.(matrix)
    sign = BigInt(1)
    previous_pivot = BigInt(1)

    for pivot_index = 1:(dimension-1)
        pivot_row = pivot_index
        while pivot_row <= dimension && iszero(work[pivot_row, pivot_index])
            pivot_row += 1
        end
        pivot_row > dimension && return BigInt(0)

        if pivot_row != pivot_index
            for column = 1:dimension
                work[pivot_index, column], work[pivot_row, column] =
                    work[pivot_row, column], work[pivot_index, column]
            end
            sign = -sign
        end

        pivot = work[pivot_index, pivot_index]
        for row = (pivot_index+1):dimension
            for column = (pivot_index+1):dimension
                numerator =
                    work[row, column] * pivot -
                    work[row, pivot_index] * work[pivot_index, column]
                quotient, remainder = divrem(numerator, previous_pivot)
                iszero(remainder) ||
                    throw(ArgumentError("support-cone determinant lost exact divisibility"))
                work[row, column] = quotient
            end
            work[row, pivot_index] = 0
        end
        previous_pivot = pivot
    end

    return sign * work[dimension, dimension]
end

function _supportcone_adjugate(matrix::AbstractMatrix{<:Integer}, determinant::BigInt)
    size(matrix, 1) == size(matrix, 2) ||
        throw(ArgumentError("a support-cone adjugate requires a square matrix"))
    dimension = size(matrix, 1)
    dimension == 0 && return Matrix{BigInt}(undef, 0, 0)

    augmented = Matrix{Rational{BigInt}}(undef, dimension, 2 * dimension)
    for row = 1:dimension
        for column = 1:dimension
            augmented[row, column] = BigInt(matrix[row, column]) // BigInt(1)
            augmented[row, dimension + column] =
                row == column ? BigInt(1) // BigInt(1) : BigInt(0) // BigInt(1)
        end
    end

    for pivot_index = 1:dimension
        pivot_row = pivot_index
        while pivot_row <= dimension && iszero(augmented[pivot_row, pivot_index])
            pivot_row += 1
        end
        pivot_row > dimension &&
            throw(ArgumentError("a support-cone adjugate requires a nonsingular matrix"))

        if pivot_row != pivot_index
            for column = 1:(2 * dimension)
                augmented[pivot_index, column], augmented[pivot_row, column] =
                    augmented[pivot_row, column], augmented[pivot_index, column]
            end
        end

        pivot = augmented[pivot_index, pivot_index]
        for column = 1:(2 * dimension)
            augmented[pivot_index, column] /= pivot
        end
        for row = 1:dimension
            row == pivot_index && continue
            factor = augmented[row, pivot_index]
            iszero(factor) && continue
            for column = 1:(2 * dimension)
                augmented[row, column] -= factor * augmented[pivot_index, column]
            end
        end
    end

    adjugate = Matrix{BigInt}(undef, dimension, dimension)
    for row = 1:dimension
        for column = 1:dimension
            adjugate[row, column] =
                BigInt(determinant * augmented[row, dimension + column])
        end
    end
    return adjugate
end

function _supportcone_lexicographic_basis(::Val{N}, support) where {N}
    pivots = Int[]
    echelon = Vector{Vector{Rational{BigInt}}}()
    indices = Int[]

    for (index, direction) in enumerate(support)
        row = [BigInt(value) // BigInt(1) for value in direction]
        for (pivot, basis_row) in zip(pivots, echelon)
            factor = row[pivot]
            iszero(factor) && continue
            for column = pivot:N
                row[column] -= factor * basis_row[column]
            end
        end

        pivot = findfirst(value -> !iszero(value), row)
        pivot === nothing && continue
        scale = row[pivot]
        for column = pivot:N
            row[column] /= scale
        end
        push!(pivots, pivot)
        push!(echelon, row)
        push!(indices, index)
        length(indices) == N && return Tuple(indices)
    end

    return nothing
end

function _supportcone_subcone(support, indices::NTuple{N,Int}) where {N}
    matrix = Matrix{BigInt}(undef, N, N)
    for row = 1:N
        for column = 1:N
            matrix[row, column] = support[indices[column]][row]
        end
    end
    determinant = _supportcone_determinant(matrix)
    iszero(determinant) && return nothing
    return SupportConeSubcone(indices, determinant, _supportcone_adjugate(matrix, determinant))
end

function _supportcone_candidate_subcones(::Val{N}, support) where {N}
    subcones = SupportConeSubcone{N}[]
    basis = _supportcone_lexicographic_basis(Val(N), support)
    if basis !== nothing
        subcone = _supportcone_subcone(support, basis)
        subcone === nothing || push!(subcones, subcone)
    end

    combinations_examined = 0
    for selection in Combinatorics.combinations(1:length(support), N)
        combinations_examined >= _SUPPORTCONE_MAX_SUBCONE_COMBINATIONS && break
        combinations_examined += 1

        indices = Tuple(selection)
        indices == basis && continue
        subcone = _supportcone_subcone(support, indices)
        subcone === nothing && continue
        push!(subcones, subcone)
        length(subcones) >= _SUPPORTCONE_MAX_SUBCONES && break
    end

    return (; subcones, combinations_examined)
end

function _supportcone_planner(::Val{N}, raw_support) where {N}
    support = NTuple{N,Int}[Tuple(direction) for direction in raw_support]
    sort!(support)
    unique!(support)

    candidates = _supportcone_candidate_subcones(Val(N), support)
    return SupportConePlanner{N}(support, candidates.subcones)
end

function SupportConePlanner(support)
    raw_support = [Tuple(Int.(direction)) for direction in support]
    isempty(raw_support) && throw(ArgumentError("support must contain a direction"))
    dimension = length(first(raw_support))
    dimension > 0 || throw(ArgumentError("support directions must be nonempty"))
    all(direction -> length(direction) == dimension, raw_support) ||
        throw(ArgumentError("support directions must have equal dimension"))
    return _supportcone_planner(Val(dimension), raw_support)
end

function _supportcone_floors(subcone::SupportConeSubcone{N}, u) where {N}
    length(u) == N || return nothing

    numerators = Vector{BigInt}(undef, N)
    for coordinate = 1:N
        numerator = BigInt(0)
        for row = 1:N
            numerator += subcone.adjugate[coordinate, row] * u[row]
        end
        numerators[coordinate] = numerator
    end

    denominator = subcone.determinant
    contains = if denominator > 0
        all(numerator -> numerator >= 0, numerators)
    else
        all(numerator -> numerator <= 0, numerators)
    end
    contains || return nothing

    return BigInt[fld(numerator, denominator) for numerator in numerators]
end

function _supportcone_plan(
    planner::SupportConePlanner{N},
    key::ReductionChainKey,
    terminal_order::Integer,
    directions::ReductionDirectionSet,
) where {N}
    required_steps = key.order - terminal_order
    required_steps > 0 || return nothing

    for subcone in planner.subcones
        chunk_lengths = _supportcone_floors(subcone, key.u)
        chunk_lengths === nothing && continue
        sum(chunk_lengths) >= required_steps || return nothing

        excess = sum(chunk_lengths) - required_steps
        while excess > 0
            largest = 1
            for index = 2:N
                chunk_lengths[index] > chunk_lengths[largest] && (largest = index)
            end
            chunk_lengths[largest] -= 1
            excess -= 1
        end

        chunks = SupportConeChunk{N}[]
        current_u = collect(key.u)
        for (coordinate, support_index) in enumerate(subcone.indices)
            steps = Int(chunk_lengths[coordinate])
            steps == 0 && continue
            direction = planner.support[support_index]
            directions(Tuple(current_u), steps, direction) || return nothing
            push!(chunks, SupportConeChunk(direction, steps))
            for index = 1:N
                current_u[index] -= steps * direction[index]
            end
        end

        return SupportConePlan(subcone.indices, chunks)
    end

    return nothing
end

"""
SupportConePolicy(planner, degree, S)

Experimental controlled-reduction policy using exact plans from the defining
polynomial's support, with depth-first fallback for keys without a support plan.
Support vectors are sorted lexicographically, and the bounded candidate
subcones are tried in lexicographic support-index order; overlapping subcones
are not triangulated.
"""
struct SupportConePolicy{N} <: AbstractReductionPolicy
    planner::SupportConePlanner{N}
    degree::Int
    S::Vector{Int}
end

function SupportConePolicy(planner::SupportConePlanner{N}, degree::Integer, S) where {N}
    degree >= 0 || throw(ArgumentError("the reduction degree must be nonnegative"))
    return SupportConePolicy{N}(planner, Int(degree), Int.(collect(S)))
end

function _supportcone_planner(f)
    return SupportConePlanner([
        Tuple(exponent_vector(monomial, 1)) for monomial in terms(f)
    ])
end

function _supportcone_fallback_chunk(policy, key, terminal_order, directions)
    maximum_chunk = key.order - terminal_order
    available_degree = sum(exponent for exponent in key.u if exponent > 0)
    available_degree >= policy.degree || throw(
        ArgumentError(
            "supportcone fallback cannot form a degree-$(policy.degree) direction from key $(key.u)",
        ),
    )
    direction = rev_chooseV(collect(key.u), policy.degree, policy.S)
    chunk = _depthfirst_chunk_length(key.u, direction, maximum_chunk)
    chunk >= 1 || throw(ArgumentError("supportcone fallback selected an empty chunk"))
    directions(key.u, chunk, direction) ||
        throw(ArgumentError("supportcone fallback selected an invalid reduction direction"))
    return direction, chunk
end

function run_reduction_policy!(
    policy::SupportConePolicy{N},
    engine::AbstractReductionEngine,
    keys_by_order,
    directions::ReductionDirectionSet,
    on_terminal,
) where {N}
    terminal_order = reduction_engine_geometry(engine).terminal_order
    pending_by_order = Dict{Int,Vector{Any}}()
    plan_suffixes = Dict{ReductionChainKey{N},Vector{SupportConeChunk{N}}}()

    for (order, inputs) in sort(collect(keys_by_order); by = first, rev = true)
        isempty(inputs) && continue
        chains = get!(pending_by_order, order, Any[])
        for input in inputs
            key = ReductionChainKey(input.key[1], input.key[2])
            chain = materialize(engine, key, input)
            push!(chains, chain)
        end
    end

    while !isempty(pending_by_order)
        order = maximum(keys(pending_by_order))
        chains = pop!(pending_by_order, order)
        key_counts = Dict{Tuple,Int}()
        for chain in chains
            key = reduction_chain_key(engine, chain)
            signature = (key.u, key.order)
            key_counts[signature] = get(key_counts, signature, 0) + 1
        end
        duplicate_keys = Set(signature for (signature, count) in key_counts if count > 1)
        survivors = combine!(engine, order)

        for chain in survivors
            key = reduction_chain_key(engine, chain)
            key.order < terminal_order &&
                throw(ArgumentError("supportcone chain is below the terminal order"))
            if key.order == terminal_order
                delete!(plan_suffixes, key)
                on_terminal(finish!(engine, chain))
                continue
            end

            signature = (key.u, key.order)
            signature in duplicate_keys && delete!(plan_suffixes, key)
            chunks = pop!(plan_suffixes, key, nothing)
            if chunks === nothing
                plan = _supportcone_plan(policy.planner, key, terminal_order, directions)
                if plan === nothing
                    direction, chunk =
                        _supportcone_fallback_chunk(policy, key, terminal_order, directions)
                    record_fallback!(engine, key, chunk)
                    reduce!(engine, chain, direction, chunk)
                    next_key = reduction_chain_key(engine, chain)
                    push!(get!(pending_by_order, next_key.order, Any[]), chain)
                    continue
                end
                chunks = plan.chunks
            end

            isempty(chunks) &&
                throw(ArgumentError("supportcone plan has no moves for a nonterminal key"))
            next_chunk = first(chunks)
            directions(key.u, next_chunk.steps, next_chunk.direction) ||
                throw(ArgumentError("supportcone selected an invalid reduction direction"))
            reduce!(engine, chain, next_chunk.direction, next_chunk.steps)
            next_key = reduction_chain_key(engine, chain)
            if length(chunks) > 1
                plan_suffixes[next_key] = chunks[2:end]
            end
            push!(get!(pending_by_order, next_key.order, Any[]), chain)
        end
    end

    return nothing
end
