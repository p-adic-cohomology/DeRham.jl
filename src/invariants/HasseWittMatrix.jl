"""
    hasse_witt_matrix(FM, h)

Given a Frobenius matrix `FM` (computed at precision 1, i.e. over a residue
ring `Z/p^M Z`) and the curve's genus/Hodge number `h`, returns the
Hasse-Witt matrix: the last `h x h` block of `FM`, reduced mod `p` (the
unique prime factor of the modulus of `base_ring(FM)`).
"""
function hasse_witt_matrix(FM::MatElem, h::Integer)
    n = modulus(base_ring(FM))
    (p, _) = first(factor(ZZ(n)))

    highprec = FM[(end-h+1):end, (end-h+1):end]

    K = GF(Int(p))

    map_entries(K, lift.(highprec))
end

"""
    hasse_witt_matrix(f; verbose = 0, kwargs...)

[`hasse_witt_matrix`](@ref) for the smooth curve defined by the homogeneous
polynomial `f`, computed via `frobenius_matrix_with_precision(f, r_m;
kwargs...)` (precision `r_m = [1, 0, ..., 0]`) and delegating to the matrix
method.
"""
function hasse_witt_matrix(f::MPolyRingElem; verbose = 0, kwargs...)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    if d < n
        throw(
            ArgumentError(
                "Hasse-Witt matrices are no supported for varieties of negative Kodaira dimensin (i.e. degree < number of variables",
            ),
        )
    end

    r_m = fill(0, n)
    r_m[1] = 1

    basis = get_basis_of_cohomology(f)

    hodge_polygon = _hodge_polygon(basis, n)
    hodge_numbers = hodge_polygon.slopelengths
    fh = hodge_numbers[1]

    F = frobenius_matrix_with_precision(f, r_m, basis = basis, verbose = verbose, kwargs...)

    hasse_witt_matrix(F, fh)
end
