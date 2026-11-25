# A numerator g of pole order k is represented by W^k g, deg(g) = kd - n.
# Public-facing basis exponents are those of x_1...x_n g; thus they are positive.
struct _SparseDworkRelation{n}
    remainder::Vector{Tuple{NTuple{n,Int},BigInt}}
    cofactors::Vector{Tuple{Int,NTuple{n,Int},BigInt}}
end

struct _SparseDworkReduction{n}
    degree::Int
    cutoff::Int
    modulus::BigInt
    positions::Dict{NTuple{n,Int},Int}
    relations::Dict{NTuple{n,Int},_SparseDworkRelation{n}}
end

function _sparse_dwork_monomials(::Val{n}, degree::Int) where {n}
    result = NTuple{n,Int}[]
    degree < 0 && return result
    exponents = zeros(Int, n)
    function append_monomials(i, left)
        if i == n
            exponents[i] = left
            push!(result, Tuple(exponents))
        else
            for a = left:-1:0
                exponents[i] = a
                append_monomials(i + 1, left - a)
            end
        end
    end
    append_monomials(1, degree)
    return result
end

# Positive keys are monomial columns. Negative keys record the combination
# of original Jacobian rows. No division by a power of p takes place.
function _sparse_dwork_axpy!(row, multiple, pivot, modulus)
    for (j, value) in pivot
        a = mod(get(row, j, 0) - multiple * value, modulus)
        if iszero(a)
            delete!(row, j)
        else
            row[j] = a
        end
    end
end

function _sparse_dwork_degree!(
    relations,
    free,
    partials,
    degree,
    d,
    p,
    modulus,
    ::Val{n},
) where {n}
    cols = _sparse_dwork_monomials(Val(n), degree)
    colindex = Dict(u => j for (j, u) in enumerate(cols))
    cofactors = _sparse_dwork_monomials(Val(n), degree - d + 1)
    generators = [(i, g) for i = 1:n for g in cofactors]
    rows = Dict{Int,BigInt}[]
    for (j, (i, g)) in enumerate(generators)
        row = Dict(-j => big(1))
        for (u, a) in partials[i]
            e = ntuple(l -> g[l] + u[l], Val(n))
            row[colindex[e]] = a
        end
        push!(rows, row)
    end
    # Unit pivots modulo p lift the same monomial complement to Z/p^N.
    pivots = zeros(Int, length(cols))
    rank = 0
    for j in eachindex(cols)
        q = findnext(row -> !iszero(mod(get(row, j, 0), p)), rows, rank + 1)
        isnothing(q) && continue
        rank += 1
        rows[rank], rows[q] = rows[q], rows[rank]
        pivot = rows[rank]
        scale = invmod(pivot[j], modulus)
        for l in keys(pivot)
            pivot[l] = mod(pivot[l] * scale, modulus)
        end
        for (i, row) in enumerate(rows)
            if i != rank && haskey(row, j)
                _sparse_dwork_axpy!(row, row[j], pivot, modulus)
            end
        end
        pivots[j] = rank
    end
    for i = (rank+1):length(rows)
        any(>(0), keys(rows[i])) &&
            throw(ArgumentError("Jacobian quotient is not free at degree $degree"))
    end
    for (j, u) in enumerate(cols)
        rem = Tuple{NTuple{n,Int},BigInt}[]
        cof = Tuple{Int,NTuple{n,Int},BigInt}[]
        if iszero(pivots[j])
            push!(free, u)
            push!(rem, (u, big(1)))
        else
            for (l, a) in rows[pivots[j]]
                if l < 0
                    i, g = generators[-l]
                    push!(cof, (i, g, a))
                elseif l != j
                    push!(rem, (cols[l], mod(-a, modulus)))
                end
            end
        end
        relations[u] = _SparseDworkRelation(rem, cof)
    end
    return length(cols) - rank
end

function _sparse_dwork_basis(
    p::Int,
    N::Int,
    support::Vector{NTuple{n,Int}},
    coeffs::Vector{Int},
) where {n}
    isempty(support) && throw(ArgumentError("empty hypersurface equation"))
    d = sum(first(support))
    if n < 2 || d < 2 || p < 5 || iszero(d % p) || N < 1
        throw(ArgumentError("expected n,d ≥ 2, p ≥ 5 prime to d, and N ≥ 1"))
    end
    if length(support) != length(coeffs) ||
       !allunique(support) ||
       any(e -> sum(e) != d || minimum(e) < 0, support) ||
       any(a -> iszero(a % p), coeffs)
        throw(
            ArgumentError(
                "expected distinct homogeneous monomials with nonzero coefficients",
            ),
        )
    end
    modulus = big(p)^N
    ahat = [powermod(big(c), modulus ÷ p, modulus) for c in coeffs]
    partials = [Dict{NTuple{n,Int},BigInt}() for _ = 1:n]
    for (e, a) in zip(support, ahat), i = 1:n
        iszero(e[i]) && continue
        u = ntuple(l -> e[l] - (l == i), Val(n))
        b = mod(a * e[i], modulus)
        iszero(b) || (partials[i][u] = b)
    end
    # For a smooth hypersurface with p not dividing d, the partials form a
    # regular sequence of degree d-1. The Jacobian quotient has socle degree
    # n(d-2); conversely surjectivity above this degree excludes every
    # common projective zero of the partials over the algebraic closure.
    cutoff = n * (d - 2) + 1
    relations = Dict{NTuple{n,Int},_SparseDworkRelation{n}}()
    free = NTuple{n,Int}[]
    missing =
        _sparse_dwork_degree!(relations, free, partials, cutoff, d, p, modulus, Val(n))
    iszero(missing) || throw(ArgumentError("the hypersurface is geometrically singular"))
    for e = mod(-n, d):d:(cutoff-1)
        _sparse_dwork_degree!(relations, free, partials, e, d, p, modulus, Val(n))
    end
    expected = ((big(d - 1)^n + (-1)^n * (d - 1)) ÷ d)
    length(free) == expected || error("primitive rank disagrees with the smooth rank")
    sort!(free; by = u -> (sum(u), u))
    basis = [ntuple(i -> u[i] + 1, Val(n)) for u in free]
    positions = Dict(u => i for (i, u) in enumerate(free))
    return basis, _SparseDworkReduction(d, cutoff, modulus, positions, relations)
end

function _sparse_dwork_reduce(raw, data::_SparseDworkReduction{n}) where {n}
    modulus = data.modulus
    layers = Dict{Int,Dict{NTuple{n,Int},BigInt}}()
    for (v, a) in raw
        length(v) == n && minimum(v) >= 1 && sum(v) % data.degree == 0 ||
            throw(ArgumentError("invalid positive projective exponent"))
        u = ntuple(i -> v[i] - 1, Val(n))
        layer = get!(() -> Dict{NTuple{n,Int},BigInt}(), layers, sum(v) ÷ data.degree)
        layer[u] = mod(get(layer, u, 0) + a, modulus)
    end
    out = zeros(BigInt, length(data.positions))
    for k = maximum(keys(layers); init = 0):-1:1
        haskey(layers, k) || continue
        lower = get!(() -> Dict{NTuple{n,Int},BigInt}(), layers, k - 1)
        for (u, a) in pop!(layers, k)
            iszero(a) && continue
            # Factor u = x + h with |h| = cutoff in high degree. Below the
            # cutoff h=u. The stored coefficients express x^h in terms of
            # basis monomials and multiples x^g ∂_i f, modulo p^N.
            left = data.cutoff
            h = ntuple(Val(n)) do i
                z = min(left, u[i])
                left -= z
                z
            end
            x = ntuple(i -> u[i] - h[i], Val(n))
            relation = data.relations[h]
            for (b, c) in relation.remainder
                j = data.positions[b]
                out[j] = mod(out[j] + a * c, modulus)
            end
            # W^k x^(x+g) ∂_i f ≡ -W^(k-1) ∂_i x^(x+g).
            for (i, g, c) in relation.cofactors
                m = x[i] + g[i]
                iszero(m) && continue
                v = ntuple(l -> x[l] + g[l] - (l == i), Val(n))
                lower[v] = mod(get(lower, v, 0) - a * c * m, modulus)
            end
        end
    end
    return out
end
