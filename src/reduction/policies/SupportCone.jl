"""
    SupportConePlanner(support)

Build exact support directions and all nonsingular support subcones. Directions
are exponent tuples sorted lexicographically. Subcones are ordered
lexicographically by their support-index tuples, and may overlap; they are not
a triangulation. Determinants and Cramer coordinates use integer arithmetic.
"""
struct SupportConeSubcone{N}
    indices::NTuple{N,Int}
    determinant::BigInt
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

function _supportcone_determinant(matrix::AbstractMatrix{<:Integer})
    size(matrix, 1) == size(matrix, 2) ||
        throw(ArgumentError("a support-cone determinant requires a square matrix"))
    dimension = size(matrix, 1)
    dimension == 0 && return BigInt(1)

    work = BigInt.(matrix)
    sign = BigInt(1)
    previous_pivot = BigInt(1)

    for pivot_index in 1:(dimension-1)
        pivot_row = pivot_index
        while pivot_row <= dimension && iszero(work[pivot_row, pivot_index])
            pivot_row += 1
        end
        pivot_row > dimension && return BigInt(0)

        if pivot_row != pivot_index
            for column in 1:dimension
                work[pivot_index, column], work[pivot_row, column] =
                    work[pivot_row, column], work[pivot_index, column]
            end
            sign = -sign
        end

        pivot = work[pivot_index, pivot_index]
        for row in (pivot_index+1):dimension
            for column in (pivot_index+1):dimension
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

function _supportcone_planner(::Val{N}, raw_support) where {N}
    support = NTuple{N,Int}[Tuple(direction) for direction in raw_support]
    sort!(support)
    unique!(support)

    subcones = SupportConeSubcone{N}[]
    for selection in Combinatorics.combinations(1:length(support), N)
        indices = Tuple(selection)
        matrix = Matrix{BigInt}(undef, N, N)
        for row in 1:N
            for column in 1:N
                matrix[row, column] = support[indices[column]][row]
            end
        end
        determinant = _supportcone_determinant(matrix)
        iszero(determinant) || push!(subcones, SupportConeSubcone(indices, determinant))
    end

    return SupportConePlanner{N}(support, subcones)
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

function _supportcone_floors(
    planner::SupportConePlanner{N},
    subcone::SupportConeSubcone{N},
    u,
) where {N}
    length(u) == N || return nothing

    matrix = Matrix{BigInt}(undef, N, N)
    for row in 1:N
        for column in 1:N
            matrix[row, column] = planner.support[subcone.indices[column]][row]
        end
    end

    numerators = Vector{BigInt}(undef, N)
    for column in 1:N
        replaced = copy(matrix)
        for row in 1:N
            replaced[row, column] = u[row]
        end
        numerators[column] = _supportcone_determinant(replaced)
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
        chunk_lengths = _supportcone_floors(planner, subcone, key.u)
        chunk_lengths === nothing && continue
        sum(chunk_lengths) >= required_steps || return nothing

        excess = sum(chunk_lengths) - required_steps
        while excess > 0
            largest = 1
            for index in 2:N
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
            for index in 1:N
                current_u[index] -= steps * direction[index]
            end
        end

        return SupportConePlan(subcone.indices, chunks)
    end

    return nothing
end
