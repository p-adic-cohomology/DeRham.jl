function compute_basis_matrix(f, l, m, params, cache)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    vars = gens(PR)
    R = base_ring(PR)

    @assert(0 <= m && m <= n)

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
function compute_monomial_basis(f, m, params, cache)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    vars = gens(PR)

    ev = cache[m*d-n-1]#gen_exp_vec(n + 1, m*d - n - 1, params.termorder)

    if cache.vars_reversed
        reverse!.(ev)
    end

    row_monomials = gen_mon(ev, PR)

    if cache.vars_reversed
        reverse!.(ev)
    end

    M = compute_basis_matrix(f, d*m - n - 1, m, params, cache)
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

# Computes the the monomial bases for different `m`. That is,
# `compute_monomial_bases(f, R, PR)[m]` will give the `m`-th case.
function compute_monomial_bases(f, params, cache)
    n = nvars(parent(f)) - 1

    res = []

    for m = 1:n
        b = compute_monomial_basis(f, m, params, cache)
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

    basis = compute_monomial_bases(f, params, cache) # basis of cohomology

    if basis == nothing
        (0 < verbose) &&
            println("Cannont compute monomial basis, this f appears to be non-smooth")
        return false
    end

    Basis = []
    for i = 1:n
        for j in basis[i]
            push!(Basis, [j, i])
        end
    end

    (basis, Basis)
end

function get_basis_of_cohomology(f, S, params, cache)
    get_basis_of_cohomology_twoflavors(f, S, params, cache)[2]
end
