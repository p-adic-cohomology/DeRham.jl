
"""
    computeD(N, m)

Returns a list of length N where D_{j,m} = sum_{i=j}^{N-1} (-1)^{i+j}binom{-m}{i}binom{i}{j}

INPUTS:
* "N" -- integer
* "m" -- integer
"""
function computeD(N, m)
    D = zeros(Int, N)
    for j = 0:(N-1)
        D[j+1] = sum((-1)^(i+j)*binomial(-m, i)*binomial(i, j) for i = j:(N-1))
    end
    return D
end

function _generalized_binomial(a::QQFieldElem, i::Integer)
    value = QQ(1)
    for t in 0:(i-1)
        value *= (a - t) / QQ(t + 1)
    end
    return value
end

"""
    computeD_half_integral(N, k)

Return the length-`N` Frobenius series coefficients for the pole exponent `k`.
The coefficients are exact rationals and use the generalized binomial
coefficient with exponent `-k/2`.
"""
function computeD_half_integral(N, k::Integer)
    a = QQ(-BigInt(k), BigInt(2))
    D = Vector{QQFieldElem}(undef, N)
    generalized_binomial_coefficients = [_generalized_binomial(a, i) for i = 0:(N-1)]
    for j = 0:(N-1)
        D[j+1] = sum(
            (-1)^(i+j) * generalized_binomial_coefficients[i+1] * binomial(i, j) for i = j:(N-1);
            init = QQ(0),
        )
    end
    return D
end

"""
    applyFrobeniusToMon(n,d,f,N,p,beta,m,R,PR)

Computes the power series expansion of p^{m-n-1}sigma(x^{beta}Omega/f^m)
using formula (1.10) in Costa's thesis


INPUTS:
* "n" -- integer, dimension of ambient projective space
* "d" -- integer, degree of the hypersurface f
* "f" -- polynomial, defining homogeneous equation of the hypersurface lifted to characteristic 0
* "N" -- integer, series precision
* "p" -- integer, a prime number that is the characteristic of the base field of the hypersurface
* "beta" -- vector, representing the exponents in the monomial of the basis element
* "m" -- integer, pole order of the basis element
* "R" -- ring, precision ring
* "PR" -- ring, polynomial ring with coefficients in R
* "termorder" -- the term ordering that should be used in vector representations
* "vars_reversed" -- reverses the order of basis vectors at various places
>>>if you don't know what this is, ignore it.
"""
function applyFrobeniusToMon(
    n,
    d,
    f,
    N,
    p,
    beta,
    m,
    R,
    PR,
    termorder,
    vars_reversed;
    verbose = false,
)
    (9 < verbose) && println("Applying Frobenius to $beta with pole order $m")
    (9 < verbose) && println("N=$N, m=$m")
    if N == 0
        return []
    end
    (9 < verbose) && println("Scaling by factorial of: ", p * (N + m - 1) - 1)
    Factorial = factorial(ZZ(p * (N + m - 1) - 1))
    o = ones(Int64, n+1)
    B = MPolyBuildCtx(PR)
    push_term!(B, R(1), o)
    X1 = finish(B)
    D = computeD(N, m)
    result = []
    for j = 0:(N-1)
        e = j + m
        factorial_e = R(ZZ(Factorial//factorial(ZZ(p * e - 1))))
        ev = gen_exp_vec(n+1, d*j, termorder)
        fj = f^j
        sum = 0
        for alpha in ev
            B = MPolyBuildCtx(PR)
            push_term!(B, R(1), Int64(p) * (beta + alpha + o))
            monomial = div(finish(B), X1)
            coefficient = R(factorial_e * (D[j+1] * (coeff(fj, alpha)))) # needs to act by Frobenius on coeff() if base field is not FF_p
            sum = sum + coefficient * monomial
            if (9 < verbose) && coefficient != 0
                Djm = D[j+1]
                C_jalpha = coeff(fj, alpha)
                println("Djm=$Djm, C_jalpha=$C_jalpha\n")
            end
        end
        push!(result, [sum, p*(m+j)])
    end
    return result
end

function _residue_ring_rational(R, q)
    q = QQ(q)
    return R(Oscar.numerator(q)) * inv(R(Oscar.denominator(q)))
end

function _frobenius_prescale_factors(fe::ZetaFrontend, m, N)
    N == 0 && return QQFieldElem[]
    terminal_index = BigInt(frobenius_term_pole_order(fe, m, N - 1))
    terminal_prescale = reduction_prescale(fe, terminal_index)
    factors = Vector{QQFieldElem}(undef, N)
    for j = 0:(N-1)
        pole_order = BigInt(frobenius_term_pole_order(fe, m, j))
        factors[j+1] = terminal_prescale / reduction_prescale(fe, pole_order)
    end
    return factors
end

function _frobenius_prescale_factors(fe::DoubleCover, m, N)
    N == 0 && return QQFieldElem[]
    P = presentation(fe)
    pole_orders = BigInt[frobenius_term_pole_order(fe, m, j) for j = 0:(N-1)]
    factors = Vector{QQFieldElem}(undef, N)
    scale = QQ(1)
    for i in N:-1:1
        factors[i] = scale
        if i > 1
            for pole_order in pole_orders[i-1]:(pole_orders[i]-1)
                scale *= inv(reduction_scalar(P, pole_order))
            end
        end
    end
    return factors
end

function applyFrobeniusToMon(
    n,
    d,
    f,
    fe::ZetaFrontend,
    N,
    p,
    beta,
    m,
    R,
    PR,
    termorder,
    vars_reversed;
    verbose = false,
    series = nothing,
    prescale_factors = nothing,
)
    (9 < verbose) && println("Applying Frobenius to $beta with pole order $m")
    (9 < verbose) && println("N=$N, m=$m")
    if N == 0
        return []
    end
    terminal_index = frobenius_term_pole_order(fe, m, N - 1)
    (9 < verbose) && println("Scaling by prescale at pole order: $terminal_index")
    rescaling = prescale_factors === nothing ? _frobenius_prescale_factors(fe, m, N) : prescale_factors
    o = ones(Int64, n + 1)
    B = MPolyBuildCtx(PR)
    push_term!(B, R(1), o)
    X1 = finish(B)
    D = series === nothing ? frobenius_series(fe, N, m) : series
    result = []
    for j = 0:(N-1)
        e = j + m
        pole_order = frobenius_term_pole_order(fe, m, j)
        factorial_e = _residue_ring_rational(R, rescaling[j+1])
        series_coefficient = _residue_ring_rational(R, D[j+1])
        ev = gen_exp_vec(n + 1, d * j, termorder)
        fj = f^j
        sum = 0
        for alpha in ev
            B = MPolyBuildCtx(PR)
            push_term!(B, R(1), Int64(p) * (beta + alpha + o))
            monomial = div(finish(B), X1)
            coefficient = factorial_e * series_coefficient * coeff(fj, alpha)
            sum = sum + coefficient * monomial
            if (9 < verbose) && coefficient != 0
                println("Djm=$(D[j+1]), C_jalpha=$(coeff(fj, alpha))\n")
            end
        end
        push!(result, [sum, pole_order])
    end
    return result
end

"""
    applyFrobeniusToBasis(Basis,n,d,f,N,p,R,PR)

Applies the frobenius to all the elements of Basis

INPUTS:
* "Basis" -- array of basis elmenets
* "n" -- number of variables minus 1
* "f" -- polynomial which is the denominator of poles (lifted version)
* "N_m" -- series precision
* "p" -- the prime
* "params" - the ControlledReductionParamaters
"""
function applyFrobeniusToBasis(Basis, f, N_m, p, params)
    termorder = params.termorder
    vars_reversed = params.vars_reversed
    verbose = params.verbose

    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    result = []
    for b in Basis
        Fmon = applyFrobeniusToMon(
            n,
            d,
            f,
            N_m[b[2]],
            p,
            exponent_vector(b[1], 1),
            b[2],
            R,
            PR,
            termorder,
            vars_reversed,
            verbose = verbose,
        )
        (9 < verbose) && println("Applying Frobenius to $b gives $Fmon")
        push!(result, Fmon)
    end
    return result
end

function applyFrobeniusToBasis(Basis, fe::ZetaFrontend, N_m, p, params)
    isempty(Basis) && return []

    termorder = params.termorder
    vars_reversed = params.vars_reversed
    verbose = params.verbose

    f = polynomial(fe)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(Basis[1][1])
    R = coefficient_ring(PR)
    fLift = liftCoefficients(R, PR, f)
    result = []
    series_cache = Dict{Int,Any}()
    prescale_cache = Dict{Int,Vector{QQFieldElem}}()
    if fe isa DoubleCover
        for b in Basis
            m = b[2]
            N = N_m[m]
            if N != 0 && !haskey(series_cache, m)
                series_cache[m] = frobenius_series(fe, N, m)
            end
        end
    end
    for b in Basis
        m = b[2]
        N = N_m[m]
        if N != 0 && !haskey(prescale_cache, m)
            prescale_cache[m] = _frobenius_prescale_factors(fe, m, N)
        end
    end
    for b in Basis
        m = b[2]
        Fmon = applyFrobeniusToMon(
            n,
            d,
            fLift,
            fe,
            N_m[b[2]],
            p,
            exponent_vector(b[1], 1),
            b[2],
            R,
            PR,
            termorder,
            vars_reversed,
            verbose = verbose,
            series = get(series_cache, m, nothing),
            prescale_factors = get(prescale_cache, m, nothing),
        )
        (9 < verbose) && println("Applying Frobenius to $b gives $Fmon")
        push!(result, Fmon)
    end
    return result
end
