
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
