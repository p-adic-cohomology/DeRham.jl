# Lift the characteristic polynomial of A = p F. The supplied columns are
# integral residues modulo p^N. No filtration-adapted basis is required.
function _sparse_dwork_lift(columns, p, N, weight)
    r = length(columns)
    orders = sort([
        minimum((iszero(a) ? N : _sparse_dwork_ord(a, p)) for a in col) for col in columns
    ])
    M = matrix(ZZ, r, r, [ZZ(columns[j][i]) for i = 1:r for j = 1:r])
    R, _ = polynomial_ring(ZZ, "T"; cached = false)
    cp = charpoly(R, M)
    # c[1] = 1. At step k, sums[j] is the exact trace of F^j for j < k.
    c = BigInt[1]
    sums = BigInt[]
    digits = N
    for k = 1:r
        raw = BigInt(coeff(cp, r-k))
        scaled = big(p)^k
        iszero(mod(raw, gcd(big(p)^digits, scaled))) ||
            error("inconsistent scaled Frobenius coefficient at index $k")
        m = big(p)^max(0, digits-k)
        residue = mod(div(raw, scaled), m)
        prev = -sum((c[i+1]*sums[k-i] for i = 1:(k-1)); init = big(0))
        bound = isqrt(big(r)^2 * big(p)^(weight*k))
        # The trace prev - k*(residue + m*j) has absolute value <= bound.
        first = cld(prev-bound-k*residue, k*m)
        last = fld(prev+bound-k*residue, k*m)
        first <= last ||
            error("Frobenius coefficient contradicts the Weil bound at index $k")
        first == last || return nothing
        ck = residue + m*first
        push!(c, ck)
        push!(sums, prev-k*ck)
        digits += orders[k]
    end
    c[end]^2 == big(p)^(r*weight) ||
        error("inconsistent determinant of primitive Frobenius")
    return ZZ.(reverse(c))
end

function _sparse_dwork_ord(x, p)
    v = 0
    while iszero(rem(x, p))
        x = div(x, p)
        v += 1
    end
    return v
end

function _sparse_dwork_run(f, p, d, ::Val{n}) where {n}
    support = NTuple{n,Int}[Tuple(exponent_vector(f, i)) for i = 1:length(f)]
    # Pair indexed exponents with indexed coefficients, including for FqMPoly.
    coeffs = [Int(lift(ZZ, coeff(f, i))) for i = 1:length(f)]
    pure = [ntuple(j -> j == i ? d : 0, n) for i = 1:n]
    if length(support) > n+1 || !all(e -> e in support, pure)
        throw(ArgumentError("expected every pure power and at most one mixed monomial"))
    end
    rank = div((d-1)^n + (-1)^n*(d-1), d)
    N = n+2
    while true
        basis, data = _sparse_dwork_basis(p, N, support, coeffs)
        length(basis) == rank || error("unexpected rank of primitive cohomology")
        images = _sparse_cartier_images(p, N, support, coeffs, basis)
        columns = [_sparse_dwork_reduce(image, data) for image in images]
        answer = _sparse_dwork_lift(columns, p, N, n-2)
        isnothing(answer) || return answer
        N += 2
    end
end

"""
    sparse_dwork_coefficients(f)
    sparse_dwork_coefficients(X::ProjectiveHypersurface)

Experimental Sperber–Voight calculation of the full primitive L-polynomial.
Return descending integer coefficients, ending in 1. The input must be a
smooth cubic or quartic in three or four variables over a prime field of
characteristic at least five. Every pure power must occur; at most one other
monomial may occur. Unsupported equations and singular inputs raise
`ArgumentError`. Precision is increased until each integer coefficient is
uniquely determined by its residue and the Weil trace bound.
"""
function sparse_dwork_coefficients(f::MPolyRingElem)
    K = coefficient_ring(parent(f))
    K isa Oscar.FinField && Oscar.degree(K) == 1 ||
        throw(ArgumentError("expected a polynomial over a prime field"))
    p = BigInt(characteristic(K))
    p >= 5 || throw(ArgumentError("expected characteristic at least five"))
    p <= typemax(Int) || throw(ArgumentError("characteristic exceeds Int range"))
    n = nvars(parent(f))
    n in (3, 4) || throw(ArgumentError("expected three or four variables"))
    !iszero(f) && is_homogeneous(f) ||
        throw(ArgumentError("expected a nonzero homogeneous form"))
    d = Int(total_degree(f))
    d in (3, 4) || throw(ArgumentError("expected a cubic or quartic"))
    return _sparse_dwork_run(f, Int(p), d, Val(n))
end

sparse_dwork_coefficients(X::ProjectiveHypersurface) =
    sparse_dwork_coefficients(polynomial(X))

"""
    sparse_dwork_zeta(f)
    sparse_dwork_zeta(X::ProjectiveHypersurface)

Return the full projective zeta function using [`sparse_dwork_coefficients`](@ref).
"""
function sparse_dwork_zeta(f::MPolyRingElem)
    cs = sparse_dwork_coefficients(f)
    return zeta_function(cs, characteristic(parent(f)), nvars(parent(f))-2)
end

sparse_dwork_zeta(X::ProjectiveHypersurface) = sparse_dwork_zeta(polynomial(X))
