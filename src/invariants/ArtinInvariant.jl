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
Keyword arguments are forwarded to `frobenius_matrix_with_precision` and then
to `zeta_coefficients`.

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

    numerator = _frobenius_power_mod(F, i, target_modulus)

    return _DividedFrobenius(numerator, i, p)
end

function _frobenius_power_mod(F::MatElem, exponent::Int, target_modulus::ZZRingElem)
    R, _ = residue_ring(ZZ, target_modulus)
    return map_entries(R, lift.(F))^exponent
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

"""
    ogus_nygaard_filtration(F, n)
    ogus_nygaard_filtration(f, n; kwargs...)

Compute the level-``n`` Ogus–Nygaard filtration for a Frobenius matrix or a
quartic over a prime field of odd characteristic. The result is a named tuple
with `bases` and `colengths`; each basis matrix has lattice generators as its
columns, with `bases[1]` representing ``N^0 = H``. The lists include every
computed level and the repeated basis at the first stable level. The colengths
are the corresponding ``p``-adic valuations of the basis determinants.
The matrix method requires a square matrix over ``\\mathbb{Z}/p^M`` with
``p^n`` dividing the modulus.

The polynomial method computes the Frobenius matrix with relative precision
``[max(n-2, 0), max(n-1, 0), n]`` before delegating to the matrix method. It
returns `false` if the Frobenius computation rejects the quartic as non-smooth.
Keyword arguments are forwarded to `frobenius_matrix_with_precision` and then
to `zeta_coefficients`.

The filtration is the one described in Ogus, “Supersingular K3 crystals”
(1979).
"""
function ogus_nygaard_filtration(F::MatElem, n::Integer)
    n = Int(n)
    @assert n >= 0
    @assert nrows(F) == ncols(F)
    @assert nrows(F) > 0

    base_modulus = modulus(base_ring(F))
    factors = factor(ZZ(base_modulus))
    @assert length(factors) == 1
    p, _ = first(factors)
    @assert p > 2
    @assert mod(base_modulus, p^n) == 0

    dimension = nrows(F)
    basis = identity_matrix(ZZ, dimension)
    bases = [basis]
    colengths = [0]
    precision_modulus = p^n

    for level = 1:n
        congruence_modulus = p^level
        F_power = _frobenius_power_mod(F, level, congruence_modulus)
        F_power_ZZ = map_entries(ZZ, lift.(F_power))

        for row = 1:dimension
            coefficients = [
                sum(F_power_ZZ[row, j] * basis[j, column] for j = 1:dimension) for
                column = 1:dimension
            ]
            kernel_basis = _congruence_kernel_basis(coefficients, congruence_modulus)
            basis = _column_lattice_hnf(basis * kernel_basis)
        end

        basis = _column_lattice_hnf(
            hcat(basis, precision_modulus * identity_matrix(ZZ, dimension)),
        )
        colength = _p_adic_valuation(abs(LinearAlgebra.det(basis)), p)
        previous_colength = colengths[end]

        if level == 1 && colength != 1
            error(
                "ogus_nygaard_filtration: unexpected colength at level 1: " *
                "got $colength; expected 1",
            )
        end

        if colength == previous_colength
            push!(bases, basis)
            push!(colengths, colength)
            break
        elseif colength == previous_colength + 1
            push!(bases, basis)
            push!(colengths, colength)
        else
            sequence = vcat(colengths, colength)
            error(
                "ogus_nygaard_filtration: unexpected colength at level $level: " *
                "sequence $sequence must increase by one until it stabilizes",
            )
        end
    end

    return (bases = bases, colengths = colengths)
end

function ogus_nygaard_filtration(f::MPolyRingElem, n::Integer; kwargs...)
    n = Int(n)
    @assert n >= 0
    @assert total_degree(f) == 4
    @assert nvars(parent(f)) == 4

    p = characteristic(parent(f))
    @assert p > 2
    @assert is_prime(p)
    @assert order(coefficient_ring(parent(f))) == p

    relative_precision = [max(n - 2, 0), max(n - 1, 0), n]
    F = frobenius_matrix_with_precision(f, relative_precision; kwargs...)
    F === false && return false
    return ogus_nygaard_filtration(F, n)
end

struct ArtinInvariantResult
    sigma::Int
    exact::Bool
end

function Base.show(io::IO, result::ArtinInvariantResult)
    if result.exact
        print(io, "σ = ", result.sigma)
    else
        print(io, "σ ≥ ", result.sigma)
    end
end

"""
    artin_invariant(f, bound; kwargs...)
    artin_invariant(F, bound)

Compute the Artin invariant ``σ`` of a supersingular quartic K3 surface via
the Ogus–Nygaard filtration. Supersingularity is a precondition and is not
checked; for a non-supersingular surface the result is meaningless and the
computation may fail.

The required positive integer `bound` sets the absolute precision ``p^B``.
The result is an [`ArtinInvariantResult`](@ref) with `sigma` and `exact` fields.
If the filtration stabilizes or its final colength is less than `bound`, then
`exact` is true and `sigma` is the Artin invariant. Otherwise `sigma` is a
lower bound and `exact` is false.

The polynomial method accepts a smooth quartic in four variables over the
prime field ``\\mathbb{F}_p`` for odd ``p``. It computes the Frobenius matrix
with relative precision `[max(bound - 2, 0), max(bound - 1, 0), bound]` and
returns `false` if the Frobenius computation rejects the quartic as
non-smooth. Keyword arguments are forwarded to `frobenius_matrix_with_precision`
and then to `zeta_coefficients`.

The matrix method trusts the caller to provide ``F`` computed with that same
relative precision. It requires a square matrix over ``\\mathbb{Z}/p^M`` with
``p^B`` dividing the modulus, where ``B = bound``.

If all of ``NS(X)`` is defined over ``\\mathbb{F}_{p^m}``, then ``Φ = F/p``
has period ``m`` and Ogus' theorem gives ``σ ≤ m - 1``. Thus `bound = m - 1`
combined with this upper bound determines ``σ``; use `bound = m` when the
filtration itself must certify an exact result.

The filtration theorem is from Ogus, “Supersingular K3 crystals” (1979).
"""
function artin_invariant(F::MatElem, bound::Integer)
    bound = Int(bound)
    @assert bound >= 1
    @assert nrows(F) == ncols(F)
    @assert nrows(F) > 0

    base_modulus = modulus(base_ring(F))
    factors = factor(ZZ(base_modulus))
    @assert length(factors) == 1
    p, _ = first(factors)
    @assert p > 2
    @assert mod(base_modulus, p^bound) == 0

    filtration = ogus_nygaard_filtration(F, bound)
    sigma = filtration.colengths[end]
    stabilized =
        length(filtration.bases) > 1 && filtration.bases[end] == filtration.bases[end-1]
    return ArtinInvariantResult(sigma, stabilized || sigma < bound)
end

function artin_invariant(f::MPolyRingElem, bound::Integer; kwargs...)
    bound = Int(bound)
    @assert bound >= 1
    @assert total_degree(f) == 4
    @assert nvars(parent(f)) == 4

    p = characteristic(parent(f))
    @assert p > 2
    @assert is_prime(p)
    @assert order(coefficient_ring(parent(f))) == p

    relative_precision = [max(bound - 2, 0), max(bound - 1, 0), bound]
    F = frobenius_matrix_with_precision(f, relative_precision; kwargs...)
    F === false && return false
    return artin_invariant(F, bound)
end

function _column_lattice_hnf(generators::ZZMatrix)
    row_hnf = hnf(transpose(generators))
    nonzero_rows = filter(
        row -> any(!iszero(row_hnf[row, column]) for column = 1:ncols(row_hnf)),
        1:nrows(row_hnf),
    )
    @assert length(nonzero_rows) == nrows(generators)
    return transpose(row_hnf[nonzero_rows, :])
end

function _congruence_kernel_basis(
    coefficients::Vector{ZZRingElem},
    congruence_modulus::ZZRingElem,
)
    dimension = length(coefficients)
    transformation = identity_matrix(ZZ, dimension)
    reduced_coefficients = copy(coefficients)

    for column = 2:dimension
        first = reduced_coefficients[1]
        next = reduced_coefficients[column]
        if iszero(first) && iszero(next)
            continue
        end

        common, first_factor, next_factor = gcdx(first, next)
        next_over_common = div(next, common)
        first_over_common = div(first, common)

        for row = 1:dimension
            old_first = transformation[row, 1]
            old_next = transformation[row, column]
            transformation[row, 1] = first_factor * old_first + next_factor * old_next
            transformation[row, column] =
                -next_over_common * old_first + first_over_common * old_next
        end

        reduced_coefficients[1] = common
        reduced_coefficients[column] = ZZ(0)
    end

    first_step = div(congruence_modulus, gcd(reduced_coefficients[1], congruence_modulus))
    diagonal = identity_matrix(ZZ, dimension)
    diagonal[1, 1] = first_step
    return transformation * diagonal
end

function _p_adic_valuation(value::ZZRingElem, p::ZZRingElem)
    remaining = abs(value)
    @assert !iszero(remaining)
    exponent = 0
    while mod(remaining, p) == 0
        remaining = divexact(remaining, p)
        exponent += 1
    end
    return exponent
end
