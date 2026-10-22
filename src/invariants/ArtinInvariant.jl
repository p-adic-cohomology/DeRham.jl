struct _DividedFrobenius{T}
    numerator::T
    exponent::Int
    p::ZZRingElem
end

"""
    divided_frobenius(F, i)
    divided_frobenius(f, i; kwargs...)

Represent ``F^i / p^i`` to absolute ``p``-adic precision 1. The matrix form
accepts a square Frobenius matrix over ``\\mathbb{Z}/p^M`` and returns an
internal object containing ``F^i`` reduced modulo ``p^{i+1}``, the exponent
``i``, and ``p``. The exponent is kept equal to ``i`` even when the numerator
is divisible by a higher power of ``p``.

The polynomial form accepts a quartic in four variables over a prime field
``\\mathbb{F}_p``. It computes the Frobenius matrix with relative precision
``[i-1, i, i+1]`` and delegates to the matrix form. It returns `false` if the
quartic does not satisfy the Frobenius computation's smoothness condition.

Both forms require a positive exponent and odd characteristic. The matrix form
also requires a square matrix over ``\\mathbb{Z}/p^M`` with ``M \\ge i+1``,
so the input contains enough precision for the numerator.

```jldoctest
julia> f = DeRham.fermat_hypersurface(4, 4, 3);

julia> divided = DeRham.divided_frobenius(f, 1);

julia> (divided.exponent, size(divided.numerator))
(1, (21, 21))
```
"""
function divided_frobenius(F::MatElem, i::Integer)
    i = Int(i)
    @assert i > 0
    @assert nrows(F) == ncols(F)

    base_modulus = modulus(base_ring(F))
    factors = factor(ZZ(base_modulus))
    @assert length(factors) == 1
    p, _ = first(factors)
    @assert p > 2

    target_modulus = p^(i + 1)
    @assert mod(base_modulus, target_modulus) == 0

    R, _ = residue_ring(ZZ, target_modulus)
    reduced_F = map_entries(R, lift.(F))
    numerator = reduced_F^i

    return _DividedFrobenius(numerator, i, p)
end

function divided_frobenius(f::MPolyRingElem, i::Integer; kwargs...)
    i = Int(i)
    @assert i > 0
    @assert total_degree(f) == 4
    @assert nvars(parent(f)) == 4

    p = characteristic(parent(f))
    @assert p > 2
    @assert is_prime(p)
    @assert order(coefficient_ring(parent(f))) == p

    F = frobenius_matrix_with_precision(f, [i - 1, i, i + 1]; kwargs...)
    F === false && return false
    return divided_frobenius(F, i)
end
