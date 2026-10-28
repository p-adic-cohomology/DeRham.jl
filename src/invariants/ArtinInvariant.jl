"""
    artin_invariant(f::MPolyRingElem; basis=nothing, kwargs...)
    artin_invariant(F)

Compute the Artin invariant ``σ`` of a supersingular quartic K3 surface from
its Frobenius matrix modulo ``p²``. The computation uses linear algebra over
``𝔽_p`` on primitive ``H²``. The matrix method returns a plain `Int`; the
polynomial method returns a plain `Int` on success, or `false` if its Frobenius
computation rejects the polynomial as non-smooth.

Supersingularity is a precondition and is not checked. A non-supersingular
input, such as an ordinary quartic, can produce an incorrect value or an
assertion error. The polynomial method requires a homogeneous quartic in four
variables over the prime field ``𝔽_p`` for an odd prime. It computes the
Frobenius matrix with relative precision `[0, 1, 2]`, using the optional `basis`
for its cohomology basis and forwarding the remaining keyword arguments to
`frobenius_matrix_with_precision`. It returns `false` if that computation
rejects the polynomial as non-smooth.

The matrix method requires a 21×21 Frobenius matrix computed with relative
precision request `[0, 1, 2]` in the chosen pole-order basis. That request
provides the precision needed by the columns used in this computation. The
matrix-level method assumes this contract and does not carry or check precision
metadata. Its first column must be divisible by ``p²`` in the chosen pole-order
basis.

The algorithm is from forthcoming work of Bhamidipati, Broe, Garzella, Nair,
Pagano, and Zhang.

# Examples

A supersingular quartic over ``𝔽_3`` with Artin invariant 2:

```jldoctest
julia> R, (x, y, z, w) = polynomial_ring(GF(3), ["x", "y", "z", "w"]);

julia> f = x^4 + y^3*z + y*z^3 + x*y^2*w + w^4;

julia> DeRham.artin_invariant(f)
2
```
"""
function artin_invariant(F::Union{MatElem,AbstractMatrix})
    @assert size(F) == (21, 21)
    F = _as_frobenius_matrix(F)

    base_modulus = modulus(base_ring(F))
    factors = factor(ZZ(base_modulus))
    @assert length(factors) == 1
    p, _ = first(factors)
    @assert p > 2

    precision_modulus = p^2
    R, _ = residue_ring(ZZ, precision_modulus)
    F_mod_p2 = map_entries(entry -> R(lift(ZZ, entry)), F)
    A = map_entries(entry -> lift(ZZ, entry), F_mod_p2)
    @assert all(iszero(A[row, 1]) for row = 1:21)

    Fp = GF(Int(p))
    Abar = map_entries(Fp, A)
    @assert rank(Abar) == 1
    @assert iszero(Abar * Abar)

    _, K = nullspace(Abar)
    @assert ncols(K) == 20
    @assert 21 - ncols(K) == 1

    K_integer = map_entries(ZZ, lift.(Ref(ZZ), K))
    B_integer = A * K_integer
    @assert all(iszero(mod(B_integer[row, column], p)) for row = 1:21, column = 1:20)
    B = map_entries(Fp, divexact.(B_integer, p))

    previous_basis = K
    previous_colength = 1
    for level = 2:22
        @assert rank(hcat(previous_basis, Abar)) == ncols(previous_basis)

        current_basis = _artin_invariant_preimage(K, B, previous_basis)
        @assert rank(hcat(previous_basis, current_basis)) == ncols(previous_basis)

        current_colength = 21 - ncols(current_basis)
        if current_colength == previous_colength
            return previous_colength
        elseif current_colength == previous_colength + 1
            previous_basis = current_basis
            previous_colength = current_colength
        else
            error(
                "artin_invariant: unexpected colength at level $level: " *
                "got $current_colength after $previous_colength; expected " *
                "stabilization or an increase by one",
            )
        end
    end

    error("artin_invariant: filtration did not stabilize")
end

function _artin_invariant_preimage(K::MatElem, B::MatElem, target::MatElem)
    _, annihilator = nullspace(transpose(target))
    constraints = transpose(annihilator) * B
    _, coefficients = nullspace(constraints)
    return K * coefficients
end

function artin_invariant(f::MPolyRingElem; basis = nothing, kwargs...)
    @assert total_degree(f) == 4
    @assert is_homogeneous(f)
    @assert nvars(parent(f)) == 4

    p = characteristic(parent(f))
    @assert p > 2
    @assert is_prime(p)
    @assert Oscar.order(coefficient_ring(parent(f))) == p

    cohomology_basis = basis === nothing ? get_basis_of_cohomology(f) : basis
    relative_precision = [0, 1, 2]
    F = frobenius_matrix_with_precision(
        f,
        relative_precision;
        basis = cohomology_basis,
        kwargs...,
    )
    F === false && return false
    return artin_invariant(F)
end
