

"""
    compute_controlled_matrix(f, l, S, R, PR, params)

Given f, m, l, over the ring PR = R[x_0, ..., x_n], and the set S subseteq {0, ..., n}, computes
the matrix for the map
(mu_0, ..., mu_n) |--> sum_{i in S} mu_i partial_i f + sum_{i notin S} mu_i x_i partial_i f
 where mu_i for 0<= i < |S| are of degree l - (d - 1) and mu_i for
 |S| <= i <= n are of degree l - d.

 fields
 ------
 f - polynomial
 l - int
 S - vector of ints
 R - ring
 PR - polynomial ring
 params - the ControlledReductionParamaters
 cache - the GradedExpCache used for this controlled reduction
"""
function compute_controlled_matrix(f, l, S, R, PR, params, cache)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    vars = gens(PR)

    len_S = length(S)
    notS = setdiff(collect(0:n), S)
    # notS = [2]
    # Stilda = [1,1,0]
    len_notS = length(notS)

    @assert(len_S >= 0 && len_S <= n+1)

    Stilda = zeros(Int, n+1)
    for i in S
        Stilda[i+1] = 1
    end

    # The code that was formerly here looked as follows:
    #
    #    in_S_mons_vec = gen_exp_vec(n+1, l-(d-1), params.termorder)
    #    not_in_S_mons_vec = gen_exp_vec(n+1, l-d, params.termorder)
    #
    # Thus, I figured that when we get them from the cache, which
    # is calculated using vars_reversed=true in gen_exp_vec,
    # I figured we'd need to reverse here:

    in_S_mons_vec = cache[l-(d-1)]
    not_in_S_mons_vec = cache[l-d]


    # if cache.vars_reversed
    #     reverse!.(in_S_mons_vec)
    #     reverse!.(not_in_S_mons_vec)
    # end

    in_S_mons = gen_mon(in_S_mons_vec, PR)
    not_in_S_mons = gen_mon(not_in_S_mons_vec, PR)

    # if cache.vars_reversed
    #     reverse!.(in_S_mons_vec)
    #     reverse!.(not_in_S_mons_vec)
    # end

    in_set_section = binomial(n + l - (d-1), n)
    not_in_set_section = binomial(n + l - d, n)
    cols = len_S * in_set_section + (n + 1 - len_S) * not_in_set_section

    if len_S > 0
        @assert(length(in_S_mons_vec) > 0)
    end
    if len_S < n+1
        @assert(length(not_in_S_mons_vec) > 0)
    end

    U = matrix_space(R, binomial(n + l, n), cols)
    M = U()

    partials = zeros(parent(f), n+1)
    for i = 1:(n+1)
        partials[i] = derivative(f, i)
    end
    #partials = [ derivative(f, i) for i in 1:n+1 ] # ∂_0, ∂_1, ∂_2

    #println("partials before possibly reversing = $partials")
    #if params.vars_reversed == true
    #    # one needs to be quite careful with the ordering of partials
    #    reverse!(partials)
    #    #partials = reverse(partials)
    #    reverse!(Stilda)
    #    #Stilda = reverse(Stilda)
    #    reverse!(vars)
    #    #vars = reverse(vars)
    #end
    #println("partials after possibly reversing = $partials")
    #println("Stilda = $Stilda")
    #println("vars = $vars")
    #println("in_S_mons = $in_S_mons")
    #println("in_S_mons_vec = $in_S_mons_vec")

    v = fill(R(0), length(cache[l]))

    col_idx = 1
    for i = 1:(n+1)
        if Stilda[i] == 1
            for monomial in eachindex(in_S_mons)
                #M[:, col_idx] = polynomial_to_vector(in_S_mons[monomial] * partials[i], n+1, params.termorder, cache, vars_reversed=cache.vars_reversed)
                polynomial_to_vector!(
                    v,
                    in_S_mons[monomial] * partials[i],
                    n+1,
                    params.termorder,
                    cache,
                )#), vars_reversed=cache.vars_reversed)
                M[:, col_idx] = v
                col_idx = col_idx + 1
            end
        else
            for monomial in eachindex(not_in_S_mons)
                #M[:, col_idx] = polynomial_to_vector(not_in_S_mons[monomial] * vars[i] * partials[i], n+1, params.termorder, cache, vars_reversed=cache.vars_reversed)
                polynomial_to_vector!(
                    v,
                    not_in_S_mons[monomial] * vars[i] * partials[i],
                    n+1,
                    params.termorder,
                    cache,
                )#, vars_reversed=cache.vars_reversed)
                M[:, col_idx] = v
                col_idx = col_idx + 1
            end
        end
    end

    return M
end


# Computes the pseudo_inverse for the controlled case.
"""
Solves the linear algebra problem in section 1.5.2 of Costa's thesis,
top of page 23.
In other words, finds a pseudo-inverse to the linear map described
there.
The map is constructed as a matrix from the polynomial f and the set S.

f - the polynomial defining the hypersurface
S - the set in [0..n] to be used for the linear algebra problem
termorder - the order of the monomials used for vectors

I think these are correct: (TODO)
R - coefficient_ring(parent(f))
PR- parent(f)
params - the ControlledReductionParamaters
cache - the GradedExpCache used for this controlled reduction
"""
function pseudo_inverse_controlled(f, S, l, R, PR, params, cache)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i = 0:n])
    fLift = liftCoefficients(ZZ, PRZZ, f)
    #println("fLift=$fLift")

    if (0 < params.verbose)
        println("Computing controlled matrix...")
        @time U = compute_controlled_matrix(fLift, l, S, ZZ, PRZZ, params, cache)
    else
        U = compute_controlled_matrix(fLift, l, S, ZZ, PRZZ, params, cache)
    end
    #U = compute_controlled_matrix(f, l, S, R, PR, params)

    #println("controlled matrix: \n$U")
    (6 < params.verbose) && println("controlled matrix: \n$U")

    #temp = size(U)
    if (0 < params.verbose)
        println("Computing pseudoinverse of matrix of size $(size(U))")
        @time flag, B = Oscar.is_invertible_with_inverse(
            matrix(R, [R(x) for x in Array(U)]),
            side = :right,
        )
    else
        flag, B = Oscar.is_invertible_with_inverse(
            matrix(R, [R(x) for x in Array(U)]),
            side = :right,
        )
    end

    (6 < params.verbose) && println("pinv mod p: \n$B")

    if flag
        return (U, B)
    else
        if S == collect(0:n) && l == ((d-2)*(n+1)+1)
            throw(ArgumentError("f is not smooth"))
        else
            throw(ArgumentError("matrix from f is not right invertible"))
        end
    end
end

"""
    pseudo_inverse_controlled_lifted(f,S,l,M,params,cache)

Solves the linear algebra problem as in
`pseudo_inverse_controlled`, but then
hensel lifts to Z/p^MZ

f - the polynomial definitng the hypersurface
S - the set in [0..n] to be used for the linear algebra problem
l - ???? we need to document this, it's something used by compute_contolled_matrix
M - the absolute precision to lift to.
params - the ControlledReductionParamaters
cache - the GradedExpCache used for this controlled reduction
"""
function pseudo_inverse_controlled_lifted(f, S, l, M, params, cache)
    PR = parent(f)
    R = coefficient_ring(PR)
    (U, Sol_fp) = pseudo_inverse_controlled(f, S, l, R, PR, params, cache)
    lift_to_int(s) = map(x -> lift(ZZ, x), s)

    if (0 < params.verbose)
        println("Lifting mod p solution to the integers")
        @time Sol_mod_p_int = lift_to_int(Sol_fp)
    else
        Sol_mod_p_int = lift_to_int(Sol_fp)
    end

    p = characteristic(PR)
    GC.gc()
    if (0 < params.verbose)
        println("Hensel lifting matrix of size $(size(Sol_fp))")
        @time henselLift(p, M, U, Sol_mod_p_int)
    else
        henselLift(p, M, U, Sol_mod_p_int)
    end
end
