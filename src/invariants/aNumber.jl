"""
    a_number(FM, h)

`a_number(FM, h) = h - rank(hasse_witt_matrix(FM, h))`.
"""
function a_number(FM::MatElem, h::Integer)
    return h - rank(hasse_witt_matrix(FM, h))
end

"""
    a_number(f; verbose = 0, kwargs...)

[`a_number`](@ref) for the smooth curve defined by the homogeneous
polynomial `f`, computed via `frobenius_matrix_with_precision(f, r_m;
kwargs...)` (precision `r_m = [1, 0, ..., 0]`) and delegating to the matrix
method. Only supported for curves.
"""
function a_number(f::MPolyRingElem; verbose = 0, kwargs...)
    n = nvars(parent(f)) - 1

    if n != 2
        throw(ArgumentError("a-numbers are only supported for curves so far"))
    end

    if 0 < verbose
        println("Calculating Hasse-Witt matrix...")
    end

    r_m = fill(0, n)
    r_m[1] = 1

    basis = get_basis_of_cohomology(f)

    hodge_polygon = _hodge_polygon(basis, n)
    hodge_numbers = hodge_polygon.slopelengths
    g = hodge_numbers[1]

    FM = frobenius_matrix_with_precision(f, r_m, basis = basis, verbose = verbose, kwargs...)

    a_number(FM, g)
end
