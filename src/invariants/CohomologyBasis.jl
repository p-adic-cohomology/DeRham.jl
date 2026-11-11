function get_basis_of_cohomology(f)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    params = default_params()
    S = collect(0:n)
    cache = controlled_reduction_cache(n, d, S, params)

    get_basis_of_cohomology(f, S, params, cache)
end

"""
    get_basis_of_cohomology(fe::ZetaFrontend)

Return the cohomology basis for `fe` as a vector of `[monomial, pole_order]`
pairs. Frontend-specific input requirements are checked before the basis is
constructed.
"""
function get_basis_of_cohomology(fe::ZetaFrontend)
    check_input(fe)
    f = polynomial(fe)
    P = presentation(fe)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    params = default_params()
    S = collect(0:terminal_pole_order(P))
    cache = controlled_reduction_cache(n, d, S, params, P)

    get_basis_of_cohomology(f, P, S, params, cache)
end
