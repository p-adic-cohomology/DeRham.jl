# Prime-field Dwork–Cartier operator in the W = πw normalization.
# Sperber–Voight, "Computing zeta functions of nondegenerate hypersurfaces
# with few monomials", (1.4), §2 and §5 (projective variant).

struct _SparseCartierTerm
    exponent::Int
    weight::Int
    order::Int
    unit::BigInt
end

function _cartier_index_limit(p, N)
    # For ell_e = [T^e] exp(T + T^p/p), (1.4) gives
    # e/p + ord_p(ell_e) >= e * (p²-3p+1)/(p²(p-1)).
    # Every factor has nonnegative weight, so any e at or past the
    # corresponding cutoff can be dropped even before cohomology reduction.
    q = big(p)
    return cld(N * q^2 * (q - 1), q^2 - 3q + 1) - 1
end

function _cartier_splitting_tables(p, N, coefficients, modulus, last)
    q = big(p)
    # F_e = ord_p(e!). The integral recurrence stores C_e = p^F_e ell_e.
    # Retention implies ord_p(C_e) < N + F_e - floor(e/p), with the last
    # difference bounded by ord_p(floor(last/p)!). Keep another N digits
    # for extracting the units of all surviving coefficients.
    guard = 0
    u = last ÷ p
    while u >= p
        u ÷= p
        guard += u
    end
    work = modulus^2 * q^guard
    powers = [powermod(big(a), modulus ÷ p, modulus) for a in coefficients]
    current = fill(big(1), length(powers))
    tables = [[_SparseCartierTerm[] for _ = 1:p] for _ in powers]
    factorial_orders = zeros(Int, last + 1)
    history = zeros(BigInt, min(last + 1, p))
    history[1] = 1
    scaled = big(1)
    for e = 0:last
        if e > 0
            unit = e
            multiplicity = 0
            while rem(unit, p) == 0
                unit ÷= p
                multiplicity += 1
            end
            factorial_orders[e+1] = factorial_orders[e] + multiplicity
            slot = mod(e, p) + 1
            if e >= p
                loss = factorial_orders[e] - factorial_orders[e-p+1]
                scaled += q^loss * history[slot]
            end
            # e ell_e = ell_(e-1) + ell_(e-p), followed only by division
            # by the unit part of e. The ring buffer still holds C_(e-p).
            scaled = mod(scaled * invmod(big(unit), work), work)
            history[slot] = scaled
        end
        if !iszero(scaled)
            unit = scaled
            order = -factorial_orders[e+1]
            while rem(unit, p) == 0
                unit ÷= p
                order += 1
            end
            weight = e + widen(p) * order
            if weight < p * N
                bucket = mod(e, p) + 1
                for j in eachindex(tables)
                    push!(
                        tables[j][bucket],
                        _SparseCartierTerm(
                            e,
                            Int(weight),
                            order,
                            mod(unit * current[j], modulus),
                        ),
                    )
                end
            end
        end
        for j in eachindex(current)
            current[j] = mod(current[j] * powers[j], modulus)
        end
    end
    return tables
end

# Solve S*r = -b simultaneously for all input monomials. Each resulting
# row gives a pivot variable in terms of the independent residue indices.
function _cartier_residue_system(support::Vector{NTuple{n,Int}}, basis, p) where {n}
    s = length(support)
    aug = zeros(Int, n, s + length(basis))
    for j in eachindex(support), i = 1:n
        aug[i, j] = mod(support[j][i], p)
    end
    for j in eachindex(basis), i = 1:n
        aug[i, s+j] = mod(-basis[j][i], p)
    end
    pivots = Int[]
    for j = 1:s
        i = length(pivots) + 1
        i > n && break
        k = findfirst(!iszero, @view aug[i:n, j])
        isnothing(k) && continue
        k += i - 1
        if k != i
            aug[k, :], aug[i, :] = aug[i, :], aug[k, :]
        end
        scale = invmod(aug[i, j], p)
        for t = j:size(aug, 2)
            aug[i, t] = mod(widen(scale) * aug[i, t], p)
        end
        for k = 1:n
            k == i && continue
            factor = aug[k, j]
            iszero(factor) && continue
            for t = j:size(aug, 2)
                aug[k, t] = mod(aug[k, t] - widen(factor) * aug[i, t], p)
            end
        end
        push!(pivots, j)
    end
    return aug, pivots, setdiff(1:s, pivots)
end

function _accumulate_cartier_terms!(
    image::Dict{NTuple{n,Int},BigInt},
    support,
    buckets,
    j,
    exponent::NTuple{n,Int},
    total_index,
    budget,
    order,
    unit,
    p,
    modulus,
) where {n}
    if j > length(support)
        pole = div(total_index, p)
        key = ntuple(i -> div(exponent[i], p), Val(n))
        value = unit * ((-1)^pole * big(p)^(pole + order))
        image[key] = mod(get(image, key, 0) + value, modulus)
        return
    end
    for term in buckets[j]
        term.weight >= budget && continue
        next = ntuple(i -> exponent[i] + term.exponent * support[j][i], Val(n))
        _accumulate_cartier_terms!(
            image,
            support,
            buckets,
            j + 1,
            next,
            total_index + term.exponent,
            budget - term.weight,
            order + term.order,
            mod(unit * term.unit, modulus),
            p,
            modulus,
        )
    end
    return
end

"""
    _sparse_cartier_images(p, N, support, coefficients, basis)

Compute all unreduced Frobenius images modulo `p^N`. Each input and output key
is a positive homogeneous exponent vector; its pole is its coordinate sum
divided by the degree. The input form and reduction relations use the same
Teichmüller lifts. This enumerates `p^(#support - rank(S mod p))` residue
classes for an input with a compatible right-hand side.
"""
function _sparse_cartier_images(
    p::Int,
    N::Int,
    support::Vector{NTuple{n,Int}},
    coefficients::Vector{Int},
    basis::Vector{NTuple{n,Int}},
) where {n}
    p >= 5 && N >= 1 || throw(ArgumentError("expected p >= 5 and positive precision"))
    length(support) == length(coefficients) && !isempty(support) ||
        throw(ArgumentError("expected a nonempty polynomial"))
    d = sum(first(support))
    d > 0 && !iszero(rem(d, p)) || throw(ArgumentError("degree must be prime to p"))
    all(e -> minimum(e) >= 0 && sum(e) == d, support) ||
        throw(ArgumentError("expected homogeneous exponents"))
    all(b -> minimum(b) > 0 && iszero(rem(sum(b), d)), basis) ||
        throw(ArgumentError("invalid positive projective basis exponents"))
    # Bound all coordinates and accumulated indices used by the Int loop.
    # Nonnegative individual weights give the same cutoff for sum_j e_j.
    max_index = _cartier_index_limit(p, N)
    max_input = maximum((sum(big, b) for b in basis); init = big(0))
    max(max_input + max_index * d, big(p) * N, max_index + 1) <= typemax(Int) ||
        throw(ArgumentError("the splitting indices exceed machine integer range"))
    modulus = big(p)^N
    tables = _cartier_splitting_tables(p, N, coefficients, modulus, Int(max_index))
    aug, pivots, free = _cartier_residue_system(support, basis, p)
    s = length(support)
    outputs = [Dict{NTuple{n,Int},BigInt}() for _ in basis]
    residues = zeros(Int, s)
    buckets = Vector{Vector{_SparseCartierTerm}}(undef, s)
    for (i, b) in enumerate(basis)
        any(!iszero, @view aug[(length(pivots)+1):n, s+i]) && continue
        h = div(sum(b), d)
        for parameters in Iterators.product((0:(p-1) for _ in free)...)
            for (j, r) in zip(free, parameters)
                residues[j] = r
            end
            for (row, j) in enumerate(pivots)
                r = aug[row, s+i]
                for k in free
                    r = Int(mod(r - widen(aug[row, k]) * residues[k], p))
                end
                residues[j] = r
            end
            for j = 1:s
                buckets[j] = tables[j][residues[j]+1]
            end
            # S*r+b = 0 mod p implies h+sum(r) = 0 mod p since p ∤ d.
            # The term has output pole k=(h+sum(e))/p and coefficient
            # (-p)^k prod(ell_e * Teich(a)^e). The budget is p times
            # the remaining absolute precision; terms at equality vanish.
            _accumulate_cartier_terms!(
                outputs[i],
                support,
                buckets,
                1,
                b,
                h,
                p * N - h,
                0,
                big(1),
                p,
                modulus,
            )
        end
        filter!(pair -> !iszero(last(pair)), outputs[i])
    end
    return outputs
end
