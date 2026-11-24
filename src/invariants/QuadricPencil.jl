function _quadric_gram_matrix(f, K, n)
    iszero(f) && throw(ArgumentError("expected a nonzero quadratic form"))
    A = Oscar.zero_matrix(K, n, n)
    half = inv(K(2))
    for (c, e) in Oscar.coefficients_and_exponents(f)
        sum(e) == 2 || throw(ArgumentError("expected a homogeneous quadratic form"))
        i, j = findfirst(!iszero, e), findlast(!iszero, e)
        A[i, j] = A[j, i] = i == j ? c : c * half
    end
    return A
end

function _quadric_pencil_data(f, g)
    parent(f) === parent(g) || throw(ArgumentError("expected forms in the same ring"))
    R = parent(f)
    K = Oscar.base_ring(R)
    K isa Oscar.FinField || throw(ArgumentError("expected a finite field"))
    isodd(Oscar.characteristic(K)) || throw(ArgumentError("expected odd characteristic"))
    n = Oscar.nvars(R)
    n >= 3 || throw(ArgumentError("expected at least three variables"))
    return _quadric_gram_matrix(f, K, n), _quadric_gram_matrix(g, K, n), K, n
end

# The pivot columns give coordinate vectors complementary to the radical.
function _quadric_pencil_correction(M, q, n)
    k, echelon = Oscar.rref(M)
    isodd(k) && return zero(q)
    k == 0 && return q^(n - 2)
    pivots = Int[findfirst(j -> !iszero(echelon[i, j]), 1:n) for i = 1:k]
    disc = (-1)^(k ÷ 2) * Oscar.det(M[pivots, pivots])
    return (Oscar.is_square(disc) ? 1 : -1) * q^(n - 2 - k ÷ 2)
end

"""
    quadric_pencil_point_count(f, g)

Count the common projective zeros of two nonzero quadratic forms over an odd
finite field. Counts singular intersections as well. Enumerates the pencil of
quadrics in the two forms (one more member than the size of the field).
"""
function quadric_pencil_point_count(f::MPolyRingElem, g::MPolyRingElem)
    A, B, K, n = _quadric_pencil_data(f, g)
    q = BigInt(Oscar.order(K))
    count = sum(q^i for i = 0:(n-3)) + _quadric_pencil_correction(A, q, n)
    for t in K
        count += _quadric_pencil_correction(t * A + B, q, n)
    end
    return count
end

"""
    quadric_pencil_zeta(f, g)

Return the full zeta function of the smooth projective intersection of two
quadrics in an odd number of variables over an odd finite field. A squarefree
binary pencil determinant of full homogeneous degree certifies smoothness.
Raises an `ArgumentError` when the input falls outside this domain.
"""
function quadric_pencil_zeta(f::MPolyRingElem, g::MPolyRingElem)
    A, B, K, n = _quadric_pencil_data(f, g)
    isodd(n) || throw(ArgumentError("the zeta formula requires an odd number of variables"))
    S, t = Oscar.polynomial_ring(K, "t")
    M = Oscar.matrix(S, n, n, [t * A[i, j] + B[i, j] for i = 1:n for j = 1:n])
    D = Oscar.det(M)
    if iszero(D) || Oscar.degree(D) < n - 1 || !isone(gcd(D, Oscar.derivative(D)))
        throw(ArgumentError("the binary pencil determinant must have only simple roots"))
    end

    R, T = Oscar.polynomial_ring(ZZ, "T")
    q = ZZ(Oscar.order(K))
    denom = prod(1 - q^i * T for i = 0:(n-3))
    scale = q^((n - 3) ÷ 2)
    sign = (-1)^((n - 1) ÷ 2)
    # At infinity the coefficient of t^(n-1) is the cofactor of A.
    # Corank one guarantees a nonzero principal cofactor in each residue field.
    minors = Dict{Int,typeof(D)}()
    function root_sign(residue_norm)
        for i = 1:n
            minor = get!(minors, i) do
                indices = [j for j = 1:n if j != i]
                sign * Oscar.det(M[indices, indices])
            end
            value = residue_norm(minor)
            if !iszero(value)
                return Oscar.is_square(value) ? 1 : -1
            end
        end
        error("no nonsingular principal minor at a simple determinant root")
    end
    for (factor, _) in Oscar.factor(D)
        h = Oscar.divexact(factor, Oscar.leading_coefficient(factor))
        epsilon = root_sign(m -> Oscar.resultant(h, m))
        denom *= 1 - epsilon * (scale * T)^Oscar.degree(h)
    end
    if Oscar.degree(D) == n - 1
        denom *= 1 - root_sign(m -> Oscar.coeff(m, n - 1)) * scale * T
    end
    return one(R) // denom
end
