struct _SparseJacobianRightInverse{N}
    cutoff::Int
    modulus::BigInt
    relations::Dict{NTuple{N,Int},Vector{Tuple{Int,NTuple{N,Int},BigInt}}}
end

function _sparse_jacobian_relation_axpy!(
    row,
    row_index,
    multiple,
    pivot,
    modulus,
    column_rows,
)
    for (j, value) in pivot
        old = get(row, j, big(0))
        a = mod(old - multiple * value, modulus)
        if iszero(a)
            delete!(row, j)
            j > 0 && delete!(column_rows[j], row_index)
        else
            if iszero(old) && j > 0
                push!(column_rows[j], row_index)
            end
            row[j] = a
        end
    end
    return row
end

"""
Build sparse Jacobian relations at degree `N * (d - 2) + 1`. Return
`nothing` when the map lacks full row rank modulo `p`.

Pass `f` in the working residue ring of modulus `p^M` when `M > 1`.
The derivatives are taken in the coefficient ring of `f` before their
coefficients are lifted.

Each row uses positive keys for target monomials and negative keys for the
weights of the gradient generators. The target part equals the weighted sum
of the generators modulo `p^M`. Elimination transforms both parts together,
so a row whose target part is a single unit monomial gives its relation.
"""
function _build_sparse_jacobian_right_inverse(f, p, M, ::Val{N}) where {N}
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
        coefficient = mod(BigInt(lift(ZZ, coeff(monomial, 1))), modulus)
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

    column_rows = [Set{Int}() for _ in eachindex(monomials)]
    for (i, row) in enumerate(rows), j in keys(row)
        j > 0 && push!(column_rows[j], i)
    end

    pivots = zeros(Int, length(monomials))
    active = trues(length(rows))
    rank = 0
    for j in eachindex(monomials)
        row_index = 0
        for i in column_rows[j]
            if active[i] && !iszero(mod(rows[i][j], p)) && (row_index == 0 || i < row_index)
                row_index = i
            end
        end
        row_index == 0 && continue
        rank += 1
        active[row_index] = false
        pivot = rows[row_index]
        scale = invmod(pivot[j], modulus)
        for l in keys(pivot)
            pivot[l] = mod(pivot[l] * scale, modulus)
        end
        for i in collect(column_rows[j])
            if i != row_index && haskey(rows[i], j)
                row = rows[i]
                _sparse_jacobian_relation_axpy!(row, i, row[j], pivot, modulus, column_rows)
            end
        end
        pivots[j] = row_index
    end

    rank == length(monomials) || return nothing

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
    return _SparseJacobianRightInverse(cutoff, modulus, relations)
end

_build_sparse_jacobian_right_inverse(f, p, M) =
    _build_sparse_jacobian_right_inverse(f, p, M, Val(nvars(parent(f))))
