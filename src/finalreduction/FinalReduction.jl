"""
    computeT(fe,Basis,M,params)

Computes the "T matrix" as in the notation in Costa's thesis.

fields
------
fe - zeta frontend
Basis - output of compute_monomial_bases
M - algorithm precision
params - the ControlledReductionParamaters
cache - the GradedExpCache used for this controlled reduction
"""
function computeT(fe::ZetaFrontend, Basis, M, params, cache)
    f = polynomial(fe)
    P = presentation(fe)
    p = characteristic(parent(f))
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i = 0:n])
    partials_index = [i for i = 1:(n+1)]

    precisionring, pi = residue_ring(ZZ, p^M)
    precisionringpoly, pvars = polynomial_ring(precisionring, ["x$i" for i = 0:n])
    f_lift = liftCoefficients(precisionring, precisionringpoly, f)

    terminal_degree = numerator_degree(P, terminal_pole_order(P))
    exp_vec = cache[terminal_degree]

    monomials = gen_mon(exp_vec, precisionringpoly)

    len = length(monomials)

    partials = [derivative(f_lift, i) for i = 1:(n+1)]

    T = zero_matrix(precisionring, 0, len)

    for m in reverse(pole_orders(P))
        l = numerator_degree(P, m)
        basis = Basis[m]
        len_basis = length(basis)
        scale = _residue_ring_rational(precisionring, reduction_prescale(fe, BigInt(m)))
        if l >= 0
            exp_vec = cache[l]
            if (l-(d-1)) >= 0
                monomials_domain = compute_monomials(
                    n+1,
                    l-(d-1),
                    precisionringpoly,
                    params.termorder,
                    cache,
                )
                len_domain = length(monomials_domain)
                change_basis, change_basis_inverse =
                    monomial_change_basis_inverse_lifted(f, l, basis, M, params, cache)
                tmp = zero_matrix(precisionring, len_basis, len)

                for j = 1:len # indexing over monomials
                    row_vec = convert_p_to_m([monomials[j]], exp_vec)
                    vec = change_basis_inverse * transpose(row_vec)

                    for k = 1:len_basis
                        tmp[k, j] = scale * vec[k, 1]
                    end

                    # standard reduction step applying Griffiths-Dwork relation
                    term = 0
                    for t = 1:(len_domain*(n+1))
                        t_partial = ceil(Int, t/len_domain)
                        t_domain = t%(len_domain)
                        if t_domain == 0
                            t_domain = len_domain
                        end
                        term =
                            term +
                            vec[t+len_basis, 1] * derivative(
                                monomials_domain[t_domain],
                                partials_index[t_partial],
                            )
                    end
                    monomials[j] = term
                end
                T = vcat(tmp, T)
            else
                tmp = zero_matrix(precisionring, len_basis, len)
                for j = 1:len # indexing over monomials
                    for k = 1:length(basis)
                        tmp[k, j] = coeff(
                            monomials[j],
                            liftCoefficients(precisionring, precisionringpoly, basis[k]),
                        )
                    end
                end
                T = vcat(tmp, T)
            end
        end
    end

    return T
end

computeT(f::MPolyRingElem, Basis, M, params, cache) =
    computeT(ProjectiveHypersurface(f), Basis, M, params, cache)

"""
    reduce_order_n_to_basis

Takes a polynomial with terminal pole order, and converts
it to a vector in terms of the basis of cohomology

The Frobenius power of p is `weight(fe)`. Applying sigma to Omega gives
`sigma(Omega) = p^n (x_0...x_n)^(p-1) Omega` on P^n for both frontends.
For a hypersurface, the residue twist divides by p, giving p^(n-1) =
p^weight; for a double cover there is no twist, giving p^n = p^weight.

fe - zeta frontend
e - the original pole order of the element before the frobenius was applied
precision - the series precision for the vector
VS - a matrix space for where the vector will end up
poly_to_reduce - a polynomial with pole that will

"""
function reduce_order_n_to_basis(
    fe::ZetaFrontend,
    p,
    e,
    precision,
    VS,
    poly_to_reduce,
    T,
    R,
    cache,
    ev;
    output_ring = R,
    final_scale = nothing,
)

    N = precision

    if N == 0
        return zero_matrix(output_ring, nrows(T), 1)
    end

    pole_order = frobenius_term_pole_order(fe, e, N - 1)
    ff = reduction_prescale(fe, BigInt(pole_order))
    if final_scale === nothing
        val_ff = valuation(ff, p)
        final_val = weight(fe) - val_ff
    else
        final_val = final_scale
        val_ff = weight(fe) - final_val
    end
    ff_invertible = _residue_ring_rational(R, ff / ZZ(p)^val_ff)

    inverse_ff = inv(ff_invertible)

    temp = VS()
    temp2 = convert_p_to_m([poly_to_reduce[1]], ev, vars_reversed = cache.vars_reversed)
    for i = 1:length(ev)
        temp[i, 1] = R(temp2[i])
    end
    temp = T * temp
    output = zero_matrix(output_ring, length(temp), 1)
    for i = 1:length(temp)
        ele = inverse_ff * temp[i]
        if 0 ≤ final_val
            output[i, 1] = output_ring(lift(ZZ, ele)) * p^final_val
        else
            lifted = divexact(lift(ZZ, ele), ZZ(p)^(-final_val))
            output[i, 1] = output_ring(lifted)
        end
    end

    output

end

function _frobenius_column_scale_exponents(fe::ZetaFrontend, p, N_m, Basis)
    scale_exponents = Int[]
    for i = 1:length(Basis)
        e = Basis[i][2]
        N = N_m[e]
        if N == 0
            push!(scale_exponents, 0)
        else
            pole_order = BigInt(frobenius_term_pole_order(fe, e, N - 1))
            if fe isa ProjectiveHypersurface && pole_order - 1 <= typemax(Int)
                factorial_valuation = _ibm_vp_factorial(Int(pole_order - 1), p)
            else
                prescale = reduction_prescale(fe, pole_order)
                factorial_valuation = valuation(prescale, p)
            end
            push!(scale_exponents, weight(fe) - factorial_valuation)
        end
    end
    scale_exponents
end

function _frobenius_output_ring(p, base_modulus, scale_exponent)
    output_modulus = BigInt(base_modulus)
    if scale_exponent >= 0
        output_modulus *= BigInt(p)^scale_exponent
    else
        output_modulus, remainder = divrem(output_modulus, BigInt(p)^(-scale_exponent))
        iszero(remainder) || error("Frobenius precision is not p-adic")
    end
    output_modulus > 1 || error("Frobenius precision must be positive")
    output_ring, = residue_ring(ZZ, output_modulus)
    output_ring
end

_frobenius_matrix_base_ring(FM::MatElem) = base_ring(FM)
_frobenius_matrix_base_ring(FM::AbstractMatrix) = parent(FM[1, 1])
_as_frobenius_matrix(FM::MatElem) = FM
_as_frobenius_matrix(FM::AbstractMatrix) = matrix(_frobenius_matrix_base_ring(FM), FM)

"""
    compute_frobenius_matrix(fe,Reductions,T)

Computes Frobenius Matrix

INPUTS:
* "fe" -- zeta frontend
* "p" -- characteristic
* "N_m" -- series precision by pole order
* "Reductions" -- output of computeReductionOfTransformLA
* "T" -- output of computeT
* "Basis" -- ?????
* "params" -- the ControlledReductionParamaters
* "cache" -- the GradedExpCache used for this controlled reduction
"""
function compute_frobenius_matrix(fe::ZetaFrontend, p, N_m, Reductions, T, Basis, params, cache)
    verbose = params.verbose
    termorder = params.termorder
    (9 < verbose) && println("Terms after controlled reduction: $Reductions")
    R = parent(T[1, 1])
    scale_exponents = _frobenius_column_scale_exponents(fe, p, N_m, Basis)
    output_ring = _frobenius_output_ring(p, modulus(R), maximum(scale_exponents))
    frob_mat_temp = []
    denomArray = []
    P = presentation(fe)
    ev = cache[numerator_degree(P, terminal_pole_order(P))]
    VS = matrix_space(R, length(ev), 1)
    for i in eachindex(Reductions)
        e = Basis[i][2] # pole order of basis element
        N = N_m[e]
        poly_to_reduce = Reductions[i][1]
        column_ring = _frobenius_output_ring(p, modulus(R), scale_exponents[i])
        temp = reduce_order_n_to_basis(
            fe,
            p,
            e,
            N,
            VS,
            poly_to_reduce,
            T,
            R,
            cache,
            ev;
            output_ring = column_ring,
            final_scale = scale_exponents[i],
        )

        push!(frob_mat_temp, temp)
    end
    frob_mat = matrix_space(output_ring, size(T, 1), length(frob_mat_temp))()
    for column in eachindex(frob_mat_temp)
        for row = 1:size(T, 1)
            frob_mat[row, column] = output_ring(lift(ZZ, frob_mat_temp[column][row, 1]))
        end
    end

    return frob_mat
end

#end
