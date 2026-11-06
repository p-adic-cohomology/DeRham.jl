
# TODO: Compute Omega in terms of differential symbol and exterior
# products (perhaps as symbols).

# TODO: For each linearly independent vector m, compute $$m * Omega / P^h$$.
# Add these to the list of basis vectors.

# Future TODO: Parallelize these algorithms either by hand or by trying to use
# more library functions (for example the Combinatorics.jl library seems promising
# and likely is automatically parallel).

# MARK - stuff copied from the old utils file

# Given an Oscar Matrix M, return a the indices of those columns whcih are pivots after being particular
# into reduced row echelon form.
function pivot_columns(M)
    rank, M = rref(M)
    res = fill(0, rank)
    ncols = size(M, 2)
    j = 1
    for i = 1:rank
        while j <= ncols && M[i, j] == 0
            j += 1
        end
        res[i] = j
        j+=1
    end
    return res
end

"""
    henselLift(p, precision, A, T)
Hensel lifts mod p solution T to the linear system AX-I=0 to mod p^precision

INPUTS:
* "p" -- integer, a prime number
* "precision" -- integer
* "A" -- matrix, integer coefficients
* "T" -- matrix, integer coefficients, satisfies AT-I=0 mod p
"""
function henselLift(p, precision, A, T)
    i = 1

    temp_T = zero_matrix(base_ring(T), size(T)...)
    copy!(temp_T, T)
    #add!(temp_T,temp_T,T)
    temp_2T = zero_matrix(base_ring(T), size(T)...)
    temp_AT = zero_matrix(base_ring(T), size(A, 1), size(T, 2))
    temp_TAT = zero_matrix(base_ring(T), size(T)...)
    while i < precision
        Oscar.Nemo.mul!(temp_2T, 2, temp_T)
        #println("Computed: $temp_2T")
        #actual = 2*T
        #println("Actual: $(2*T)")
        #println("Correct? $(temp_2T == actual)")

        Oscar.Nemo.mul!(temp_AT, A, temp_T)
        #println("Computed: $temp_AT")
        #actual = A*T
        #println("Actual: $(A*T)")
        #println("Correct? $(temp_AT == actual)")

        Oscar.Nemo.mul!(temp_TAT, temp_T, temp_AT)
        #println("Computed: $temp_TAT")
        #actual = T*(A*T)
        #println("Actual: $(T*(A*T))")
        #println("Correct? $(temp_TAT == actual)")

        Oscar.Nemo.sub!(temp_T, temp_2T, temp_TAT)
        #println("Computed: $temp_T")
        #actual = 2*T - T * (A*T)
        #println("Actual: $(2*T - T * (A*T))")
        #println("Correct? $(temp_T == actual)")

        #T = 2*T - T * (A*T)
        #println("After step $i: $(julia_signed_mod.(T,p^(i+1)))")
        #println("After step $i: $(julia_signed_mod.(temp_T,p^(i+1)))")

        i *= 2
    end
    R, pi = residue_ring(ZZ, p^precision)
    #stuff = [R(x) for x in Array(T)]
    return matrix(R, temp_T)
end

function findpivotcolumnsU(U)
    n = ncols(U)
    pivotcols = []
    for i = 1:nrows(U)
        for j = 1:ncols(U)
            if U[i, j] != 0
                push!(pivotcols, j)
                break
            end
        end
    end
    return pivotcols
end
