
"""
    chooseV(I, d, S)
Choose the direction of reduction V following Edgar's method

INPUTS:
* "I" -- list/tuple, exponents of monomials
* "d" -- integer, degree of f
* "S" -- list
"""
function chooseV(I, d, S)
    v = zeros(Int, length(I))
    n = length(I)

    sum = 0

    for i in S
        if I[n-i] > 0
            v[n-i] = v[n-i] + 1
            sum = sum + 1
            if sum >= d
                break
            end
        end
    end

    while sum < d
        for i = 0:(n-1)
            if I[n-i] > v[n-i]
                v[n-i] = v[n-i] + 1
                sum = sum + 1
                break
            end
        end
    end
    return v
end

"""
    skim_chooseV(I, d)

choose direction of reduction V by scanning from right to left
Note: this differs from Edgar's method, see chooseV

INPUTS:
* "I" -- list/tuple, exponents of monomials
* "d" -- integer, degree of f
"""
function skim_chooseV(I, d)
    V = zeros(Int, length(I))
    i = 0
    s = length(I)
    foundNonZero = false
    while i < d
        #FIXME reversed to match Costa's
        if s == 0 && foundNonZero == false
            return V
        elseif s == 0
            s = length(I)
            foundNonZero = false
        end
        if (I-V)[s] > 0
            V[s] = V[s] + 1
            i = i + 1
            foundNonZero = true
        end
        s = s-1
    end
    return V
end

"""
    rev_chooseV(I, d)

choose direction of reduction in the same way as Costa's code

INPUTS:
* "I" -- list/tuple, exponents of monomials
* "d" -- integer, degree of f
* "S" -- Vector of ints
"""
function rev_chooseV(I, d, S)
    reverse!(I)

    V = chooseV(I, d, S)

    reverse!(I)

    reverse!(V)

    return V
end

"""
    tweak(I,m)

If for a vectors of ints, I, we let |I| = sum(I[i]). This function returns I after removing an integer vector, J, with |J| = m
from it

Removes from the "front" of the vector"

INPUTS
* "I" -- vector of nonnegative integers
* "m" -- nonnegative integer
"""
function tweak(J, m)
    count = 0
    o = m
    I = copy(J)
    while m > 0
        for i in axes(I, 1)
            if (I[i] > 0)
                I[i] = I[i] - 1
                m = m - 1
                break
            end
        end
        if count > length(I)*o
            return I
        end
        count = count+1
    end
    return I
end


"""
    rev_tweak(J,m)

If for a vectors of ints, I, we let |I| = sum(I[i]). This function returns I after removing an integer vector, J, with |J| = m
from it

unlike tweak, this, removes from the "back" of the vector instead of the front.

INPUTS
* "I" -- vector of nonnegative integers
* "m" -- nonnegative integer
"""
function rev_tweak(J, m)
    count = 0
    o = m
    I = copy(J)
    while m > 0
        for i in reverse(axes(I, 1))
            if (I[i] > 0)
                I[i] = I[i] - 1
                m = m - 1
                break
            end
        end
        if count > length(I)*o
            return I
        end
        count = count+1
    end
    return I
end

"""
    undo_rev_tweak(I,p)

gives the unique vector J containing only
multiples of p such that tweak(J, |J-I|) = I

We're assuming that |J-I| < p here.

INPUTS
* "I" -- vector of nonnegative integers
* "p" -- prime number
"""
function undo_rev_tweak(I, p)

    # is this correct for large d and n?
    J = copy(I)

    k = 1
    #accumulator = 0
    while k ≤ length(J)
        while (J[k] % p) != 0
            J[k] = J[k] + 1
            #accumulator += 1
        end
        k = k + 1
    end

    # @assert accumulator == m

    J
end

"""
    costadata_of_intial_term
Returns the data used by costa's code given a polynomial term.
More specifically, given a polynomial h, it writes h in the form
h = x^{rev_tweak(U)}g

Note: this only works with a single term, so it should
only be used at the beginning of reduction

fields
------
term - monomial
g - vector
n - # of vars - 1
d - deg of f
p - prime number
S - vector
cache - the GradedExpCache used for this controlled reduction
params - the ControlledReductionParamaters
"""
function costadata_of_initial_term!(term, g, n, d, P, p, S, cache, params)
    termorder = params.termorder
    vars_reversed = params.vars_reversed

    R = base_ring(parent(term[1]))
    i = term

    Stilda = zeros(Int, n+1)
    for j in S
        Stilda[j+1] = 1
    end

    U = exponent_vector(i[1], 1) + Stilda

    gCoeff = coeff(i[1], 1)

    reduction_degree = _reduction_g_degree(P)
    g_exps = U - tweak(U, reduction_degree)  # TODO: make our code work with tweak

    ev = cache[reduction_degree]

    my_zero!(g)

    for j in axes(g, 1)
        if g_exps == ev[j]
            CUDA.@allowscalar g[j] = lift(gCoeff)
            break
        end
    end

    return (U, g)
end

function costadata_of_initial_term!(term, g, n, d, p, S, cache, params)
    P = HypersurfacePresentation(n, d)
    return costadata_of_initial_term!(term, g, n, d, P, p, S, cache, params)
end

"""
    incorporate_initial_term

Takes an array of Costa's data tuples, and a new term
to be added. If the new term already exists, it is
added, if not then it concatenates.

Note: perhaps this would be better served by a custom struct?

Note: right now, we don't have a method that just
takes an array of Costa's data tuples and consolodates it.
Perhaps this also would be well suited by a custom struct?

^^ such concerns have speed implications, so better to wait
until we're really trying to optimize this.

fields
------
costadata_arr - vector of costadata
costadata - output of costadata_of_initial_term
"""
function incorporate_initial_term!(costadata_arr, costadata)
    ind_already = false

    # if the u is already there, just add the vectors
    for i in eachindex(costadata_arr)
        (u, g) = costadata_arr[i]
        if all(u .== costadata[1])
            costadata_arr[i] = (u, g .+ costadata[2])
            ind_already = true
        end
    end

    # otherwise, push on a new term
    if !ind_already
        #(9 < verbose) && println("incorporation: u has type $(typeof(costadata[1]))")
        push!(costadata_arr, costadata)
    end
end

function remove_duplicates!(costadata_arr)
    i = 1
    while i <= (length(costadata_arr)-1)
        j = i+1
        while j <= length(costadata_arr)
            if all(costadata_arr[i][1][1] .== costadata_arr[j][1][1])
                costadata_arr[i] = (
                    (
                        costadata_arr[i][1][1],
                        costadata_arr[i][1][2] + costadata_arr[j][1][2],
                    ),
                    costadata_arr[i][2][1],
                )
                deleteat!(costadata_arr, j)
            else
                j = j + 1
            end
        end
        i = i + 1
    end
end


"""
    poly_of_end_costadata

Converts a Costa's data tuple to a polynomial with pole
after reduction.

Thus, the pole order is n.

fields
------
costadata - tuple of two vectors
PR - polynomial ring
p - prime number
d - degree of f
n - # of vars of f
params - the ControlledReductionParamaters
"""
function poly_of_end_costadata(costadata, PR, P, p, d, n, params)
    (u, g_vec) = costadata
    vars = gens(PR)
    cpu_g = Array(g_vec)
    reduction_degree = _reduction_g_degree(P)
    CUDA.@allowscalar g =
        vector_to_polynomial(cpu_g, n, reduction_degree, PR, params.termorder)

    return [prod(vars .^ u) * g, terminal_pole_order(P)]
end

function poly_of_end_costadata(costadata, PR, p, d, n, params)
    P = HypersurfacePresentation(n, d)
    return poly_of_end_costadata(costadata, PR, P, p, d, n, params)
end


"""
    poly_of_end_costadatas

Converts an array of Costa's data tuples to an array of polynomials with pole
after reduction.

Thus, the pole order is always n.

fields
------
costadata - tuple of two vectors
PR - polynomial ring
p - prime number
d - degree of f
n - # of vars of f - 1
S - vector of ints
params - the ControlledReductionParamaters
"""
function poly_of_end_costadatas(costadatas, PR, P, p, d, n, S, params)
    res = PR(0)
    vars = gens(PR)
    for costadata in costadatas
        res += poly_of_end_costadata(costadata, PR, P, p, d, n, params)[1]
    end
    XS = prod(PR(vars[i+1]) for i in S; init = PR(1))
    [[div(res, XS), terminal_pole_order(P)]]
end

function poly_of_end_costadatas(costadatas, PR, p, d, n, S, params)
    P = HypersurfacePresentation(n, d)
    return poly_of_end_costadatas(costadatas, PR, P, p, d, n, S, params)
end
