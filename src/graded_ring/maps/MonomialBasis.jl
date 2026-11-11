function compute_basis_matrix(f, l, params, cache)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    vars = gens(PR)
    R = base_ring(PR)

    section = binomial(n + l - (d-1), n)
    domain_mons = compute_monomials(n+1, l - (d - 1), PR, params.termorder, cache)#, vars_reversed=cache.vars_reversed)

    if length(domain_mons) <= 0
        return []
    end

    M = zero_matrix(R, binomial(n + l, n), (n+1) * section)

    partials = [derivative(f, i) for i = 1:(n+1)]

    #TODO: if the user makes a mistake and puts in a surface that
    #has two few variables, the error happens here.
    #We should probably double check for this in ZetaFunction.jl
    for i = 1:(n+1)
        for monomial in eachindex(domain_mons)
            vec = polynomial_to_vector(
                domain_mons[monomial] * partials[i],
                n+1,
                params.termorder,
                cache,
            )#, vars_reversed=cache.vars_reversed)
            if length(vec) == 0 && size(M, 2) != 0
                return nothing # this f is not smooth (or something)
            end
            M[:, section*(i-1)+monomial] = vec
        end
    end

    return M
end

# Computes the monomial basis of $H_{dR}^n(U_{\QQ_p})$. In particular, we find the monomials
# of degree $l = dm - n - 1$ in $F_p[x_0, \dots, x_n]$ that project onto a basis of the cokernel
# of the map computed in `compute_classical_mat()`.
function monomial_basis_in_degree(f, l, params, cache)
    n = nvars(parent(f)) - 1
    PR = parent(f)

    ev = cache[l]

    if cache.vars_reversed
        reverse!.(ev)
    end

    row_monomials = gen_mon(ev, PR)

    if cache.vars_reversed
        reverse!.(ev)
    end

    M = compute_basis_matrix(f, l, params, cache)
    if M == nothing
        return nothing
    end
    if isempty(M)
        return row_monomials
    end

    rows = size(M)[1]

    temp = size(M)

    if (0 < params.verbose)
        println("Finding pivot rows of matrix of size $(size(M))")
        @time pivot_rows = pivot_columns(transpose(M))
    else
        pivot_rows = pivot_columns(transpose(M))
    end
    non_pivot_rows = setdiff([1:rows;], pivot_rows)
    return map((i) -> row_monomials[i], non_pivot_rows)
end

function compute_monomial_basis(f, P::CohomologyPresentation, m, params, cache)
    monomial_basis_in_degree(f, numerator_degree(P, m), params, cache)
end

# Computes the the monomial bases for different `m`. That is,
# `compute_monomial_bases(f, R, PR)[m]` will give the `m`-th case.
function compute_monomial_bases(f, params, cache)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    P = HypersurfacePresentation(n, d)
    compute_monomial_bases(f, P, params, cache)
end

function compute_monomial_bases(f, P::CohomologyPresentation, params, cache)
    res = []

    for m in pole_orders(P)
        b = compute_monomial_basis(f, P, m, params, cache)
        if b == nothing
            return nothing
        end
        push!(res, b)
    end
    return res
end

## Computes the pseudo_inverse for the classical case.
##TODO: update this to reflex changes to pseudo_inverse_controlled
##it's used in standard reduction, I'll plan to take care of it then

function get_basis_of_cohomology_twoflavors(f, S, params, cache)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    P = HypersurfacePresentation(n, d)
    get_basis_of_cohomology_twoflavors(f, P, S, params, cache)
end

function get_basis_of_cohomology_twoflavors(f, P::CohomologyPresentation, S, params, cache)
    basis = compute_monomial_bases(f, P, params, cache) # basis of cohomology

    if basis == nothing
        (0 < verbose) &&
            println("Cannont compute monomial basis, this f appears to be non-smooth")
        return false
    end

    Basis = []
    for (i, m) in enumerate(pole_orders(P))
        for j in basis[i]
            push!(Basis, [j, m])
        end
    end

    (basis, Basis)
end

function get_basis_of_cohomology(f, S, params, cache)
    get_basis_of_cohomology_twoflavors(f, S, params, cache)[2]
end

function get_basis_of_cohomology(f, P::CohomologyPresentation, S, params, cache)
    get_basis_of_cohomology_twoflavors(f, P, S, params, cache)[2]
end
