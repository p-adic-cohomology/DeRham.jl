"""
    prec_vec(polygon, relative_precision)

returns the vector by adding the height of the Hodge polygon and appropriate relative precision for each position

INPUTS:
* "polygon" -- A SlopesPolygon struct describing a polygon whose values describe the
divisibility of the roots.
* "relative_precision" -- list, list of relative precisions, can be computed by "calculate_relative_precision"
"""
function prec_vec(polygon, relative_precision, verbose = 0)
    vals = Int.(polygon.values)
    hodge_numbers = polygon.slopelengths
    (0 < verbose) && println(hodge_numbers)
    (0 < verbose) && println(relative_precision)
    prec_vec = reverse(vals)
    i = 1
    num_term = 1
    for j = 1:length(hodge_numbers)
        num_term = num_term + hodge_numbers[j]
        while i < num_term
            prec_vec[i] = prec_vec[i] + relative_precision[j]
            i = i + 1
        end
    end
    prec_vec[length(prec_vec)] =
        prec_vec[length(prec_vec)] + relative_precision[length(hodge_numbers)]

    return [ZZ(x) for x in prec_vec]
end

"""
    sign_fe(weight, d, p, prec_vec, cp_coeffs)

Returns the sign of the functional equation. For a degree-one characteristic polynomial, the constant coefficient determines the sign when it is known to sufficient p-adic precision.

INPUTS:
* "weight" -- integer, motivic weight of the cohomology group
* "d" -- integer, degree of characteristic polygon = dimension of cohomology group
* "p" -- integer, prime number
* "prec_vec" -- list, prec_vec[i] = maximum p-adic valuation of the i-th coefficient of the Frobenius char poly?
* "cp_coeffs" -- list, cp_coeffs[i] = i-th coefficient of the Frobenius char poly, which a priori may be incorrect
"""
function sign_fe(weight, d, p, prec_vec, cp_coeffs)
    sign = 1
    dimension = weight
    cp = cp_coeffs #TODO: why the naming discrepancy?
    if d == 1
        k = div(dimension, 2)
        required_precision = k + (p == 2 ? 2 : 1)
        if prec_vec[1] < required_precision
            throw(ArgumentError("not enough p-adic precision to determine the sign of a degree-one L-polynomial (need c_0 mod p^$required_precision)"))
        end

        p = ZZ(p)
        p_power = p^k
        c0 = ZZ(cp[1])
        if !iszero(mod(c0, p_power))
            throw(ArgumentError("constant coefficient is not divisible by p^k for a degree-one L-polynomial"))
        end

        unit_modulus = p == 2 ? ZZ(4) : p
        unit = mod(divexact(c0, p_power), unit_modulus)
        if unit == 1
            return 1
        elseif unit == unit_modulus - 1
            return -1
        end
        throw(ArgumentError("constant coefficient does not determine a sign of ±1 for a degree-one L-polynomial"))
    end

    if dimension % 2 == 0
        for i = 1:Int(floor(d/2))
            # NOTE-TO-SELF: check if the symmetric index is d-i or d+1-i
            k = div((d-2*i)*dimension, 2)
            p_power = p^(min(prec_vec[i+1], prec_vec[d+1-i]+k))

            if (mod(cp[i+1], p_power) != 0) && (mod(cp[d+1-i], p_power) != 0)
                if 0 == mod(cp[i+1] + cp[d+1-i] * p^k, p_power)
                    sign = -1
                else
                    sign = 1
                    #println(mod(cp[i+1] - cp[d+1-i] * ZZ(p^k), p_power))
                    # commented out for K3 over F3
                    # @assert 0 == mod(cp[i+1] - cp[d+1-i] * ZZ(p^k), p_power)
                end
                break
            end
        end
    end

    return sign
end

"""
    apply_symmetry(weight, d, p, prec_vec, cp_coeffs, modulus, sign)

Update the coefficients in the Frobenius characteristic polynomial using symmetry in the coefficients coming from
the Weil conjecture/Poincare duality

INPUTS:
* "weight" -- integer, motivic weight of the cohomology group
* "d" -- integer, degree of characteristic polygon = dimension of cohomology group
* "p" -- integer, prime number
* "prec_vec" -- list, prec_vec[i] = maximum p-adic valuation of the i-th coefficient of the Frobenius char poly?
* "cp_coeffs" -- list, cp_coeffs[i] = i-th coefficient of the Frobenius char poly, which a priori may be incorrect
* "modulus" -- list, modulus[i] = p^(prec_vec[i])
* "sign" -- integer, sign of the functional equation
"""
function apply_symmetry(weight, d, p, prec_vec, cp_coeffs, modulus, sign)
    dimension = weight
    for i = 0:div(d, 2)
        k = div((d-2*i) * dimension, 2)
        if prec_vec[i+1] >= (prec_vec[d+1-i] + k)
            prec_vec[d+1-i] = prec_vec[i+1] - k
            modulus[d+1-i] = p^(prec_vec[d+1-i])
            cp_coeffs[d+1-i] = div(sign * cp_coeffs[i+1], p^k)
            cp_coeffs[d+1-i] = mod(cp_coeffs[d+1-i], modulus[d+1-i])
        else
            prec_vec[i+1] = prec_vec[d+1-i] + k
            modulus[i+1] = p^(prec_vec[i+1])
            cp_coeffs[i+1] = sign * cp_coeffs[d+1-i] * p^k
            cp_coeffs[i+1] = mod(cp_coeffs[i+1], modulus[i+1])
        end
    end
    #println("updated modulus=$modulus")
    #println("updated coeffs = $cp_coeffs")
    return cp_coeffs
end

"""
    apply_newton_identity(weight, d, p, prec_vec, cp_coeffs, modulus)
calculate the i-th power sum of the roots and correct the coefficients along the way

INPUTS:
* "weight" -- integer, motivic weight of the cohomology group
* "d" -- integer, degree of characteristic polygon = dimension of cohomology group
* "p" -- integer, prime number
* "prec_vec" -- list, prec_vec[i] = maximum p-adic valuation of the i-th coefficient of the Frobenius char poly?
* "cp_coeffs" -- list, cp_coeffs[i] = i-th coefficient of the Frobenius char poly, which a priori may be incorrect
* "modulus" -- list, modulus[i] = p^(prec_vec[i])
"""
function apply_newton_identity(weight, d, p, prec_vec, cp_coeffs, modulus)
    e = [ZZ(0) for i = 1:(d+1)]  # e[k] = (k-1)-th elementary symmetric polynomial evaluated at the roots
    s = [ZZ(0) for i = 1:(d+1)]  # s[k] = k-th power sum of the roots
    for k = 0:d
        if k%2 == 0
            e[k+1] = cp_coeffs[d+1-k]  # FIXME: check the index
        else
            e[k+1] = -cp_coeffs[d+1-k]
        end
    end

    for k = 1:d
        sum = 0
        for i = 1:(k-1)
            if i%2 == 0
                sum = sum + e[k+1-i] * s[i+1]
            else
                sum = sum - e[k+1-i] * s[i+1]
            end
        end
        #println("sum=$sum")
        if k%2 == 0
            #s[k+1] = sum - k * e[k+1]
            s[k+1] = -(sum + k * e[k+1])
        else
            #s[k+1] = -1 * (sum - k * e[k+1])
            s[k+1] = sum + k * e[k+1]
        end


        pN = ZZ(k) * modulus[d+1-k]
        s[k+1] = mod(s[k+1], pN)
        #println(s[k+1])

        if ZZ(s[k+1])^2 > (ZZ(d)*ZZ(d)*ZZ(p)^(weight*k))
            s[k+1] = -mod(-s[k+1], pN)
        end
        #println(s[k+1])

        if k%2 == 0
            #e[k+1] = div(sum-s[k+1], k)
            e[k+1] = div(-s[k+1]-sum, k)
            cp_coeffs[d+1-k] = e[k+1]
        else
            #e[k+1] = div(sum+s[k+1], k)
            e[k+1] = div(s[k+1]-sum, k)
            cp_coeffs[d+1-k] = -e[k+1]
        end
        #println("k=$k")
        #println("e=$e")
        #println("s=$s")

    end
    return cp_coeffs
end

"""
    compute_Lpolynomial(weight, d, p, polygon, relative_precision, cp_coeffs)
INPUTS:
* "weight" -- integer, motivic weight of the cohomology group
* "p" -- integer, prime number
* "polygon" -- A SlopesPolygon struct describing a polygon whose values describe the
divisibility of the roots.
* "relative_precision" -- list, list of relative precisions, can be computed by "calculate_relative_precision"
* "cp_coeffs" -- list, list of integers representing the coefficient of the characteristic polynomial of the
                        Frobenius matrix computed over the integers

"""
function compute_Lpolynomial(weight, p, polygon, relative_precision, cp_coeffs, verbose = 0)
    p = ZZ(p)
    dimension = weight
    Lpoly_coeffs = [ZZ(x) for x in cp_coeffs]
    (1 == verbose) && println("initial coefficients = $Lpoly_coeffs")
    precision_vec = prec_vec(polygon, relative_precision)
    (1 == verbose) && println("prec_vec = $precision_vec")
    d = length(cp_coeffs) - 1  # degree of characteristic polynomial
    modulus = [ZZ(0) for i = 1:(d+1)]
    for i = 1:(d+1)
        modulus[i] = p^precision_vec[i]
        Lpoly_coeffs[i] = mod(Lpoly_coeffs[i], modulus[i])
    end
    (1 == verbose) && println("coefficients after moding by prec = $Lpoly_coeffs")
    #println("coefficients after moding by prec = $Lpoly_coeffs")

    Lpoly_coeffs[d+1] = 1  # set leading coefficient to 1
    sign = sign_fe(weight, d, p, precision_vec, Lpoly_coeffs)  # determine the sign of the functional equation
    (1 == verbose) && println("sign of functional equation = $sign")
    Lpoly_coeffs[1] = sign * p^(div(d*dimension, 2))
    #println(Lpoly_coeffs[1])
    #println("prec_vec=$precision_vec")
    Lpoly_coeffs = apply_symmetry(weight, d, p, precision_vec, Lpoly_coeffs, modulus, sign)
    (1 == verbose) && println("coefficients after apply symmetry = $Lpoly_coeffs")
    #println("prec_vec=$precision_vec")
    #println("coefficients after apply symmetry = $Lpoly_coeffs")
    Lpoly_coeffs = apply_newton_identity(weight, d, p, precision_vec, Lpoly_coeffs, modulus)
    (1 == verbose) &&
        println("coeffficients after applying Newton's identity = $Lpoly_coeffs")

    return Lpoly_coeffs
end
