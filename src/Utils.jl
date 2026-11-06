
function Factorial(x, y)
    fact = 1
    if x == 0
        return 1
    end
    #BEWARE could be infinite
    while x != y
        fact = fact*x
        x = x - 1
    end
    return fact
end

"""
Performs a mathematical mod (in julia if you mod a negative number, it stays negative).
Was using for debugging.

"""
function julia_signed_mod(x, n)
    ((x % n) + n) % n
end
