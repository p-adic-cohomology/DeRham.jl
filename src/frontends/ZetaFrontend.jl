abstract type ZetaFrontend end

function polynomial end
function cohomology_quotient end
function weight end
function check_input end
function reduction_prescale end
function _supports_reduction_algorithm end
function _reduction_allowed_policies end
function reduction_ruv_function end
function reduction_directions end
function varbyvar_reduction_setup end
function reduction_input_keys end
function reduction_input_vector! end
function reduction_terminal_polynomial end
function reduction_add_terminal_chain end
function frobenius_series end
function frobenius_term_pole_order end

function reduction_prescale(P::CohomologyQuotient, m::Integer)
    m >= 1 || throw(ArgumentError("the pole order must be positive"))
    scale = QQ(1)
    for j = 1:(m-1)
        scale *= inv(reduction_scalar(P, j))
    end
    return scale
end

_reduction_allowed_policies(::ZetaFrontend) = (:depthfirst,)

_supports_reduction_algorithm(fe::ZetaFrontend, algorithm) =
    algorithm in (:default, _reduction_allowed_policies(fe)...)

"""
    reduction_ruv_function(fe, S, f, pseudoInverseMat, cache, params)

Return the function that computes the R_uv matrices for each reduction direction.
"""
function reduction_ruv_function(fe::ZetaFrontend, S, f, pseudoInverseMat, cache, params)
    P = cohomology_quotient(fe)
    return V -> computeRuvS(V, S, f, P, pseudoInverseMat, cache, params)
end

"""
    reduction_directions(fe, S)

Return a predicate `(u, k, v) -> Bool` for membership in the direction collection
V_S(u, k) of `fe`.
"""
function reduction_directions(fe::ZetaFrontend, S)
    P = cohomology_quotient(fe)
    d = P.d
    return function (u, k, v)
        length(u) == P.n + 1 || return false
        length(v) == length(u) || return false
        k >= 1 || return false
        sum(v) == d || return false
        all(x -> x >= 0, v) || return false
        all(i -> u[i] - k * v[i] >= 0, eachindex(u)) || return false
        all(i -> 0 <= i <= P.n && (u[i+1] == 0 || v[i+1] != 0), S)
    end
end

function varbyvar_reduction_setup(fe::ZetaFrontend, S)
    required_S = [cohomology_quotient(fe).n]
    return :depthfirst, S, nothing, required_S
end

"""
    reduction_input_keys(fe, frobenius_output, S)

Return the distinct `(u, order)` chain keys with their input terms, grouped by
pole order. The exponent vector `u` is shifted by `x^S` and tweaked to the
reduction degree.
"""
function reduction_input_keys(fe::ZetaFrontend, frobenius_output, S)
    P = cohomology_quotient(fe)
    reduction_degree = _reduction_g_degree(P)
    Stilda = zeros(Int, P.n + 1)
    for i in S
        Stilda[i+1] = 1
    end

    groups = Tuple{Int,Vector{Any}}[]
    group_indices = Dict{Int,Int}()
    key_indices = Dict{Any,Tuple{Int,Int}}()
    for pol in frobenius_output
        order = Int(pol[2])
        for monomial in terms(pol[1])
            U = exponent_vector(monomial, 1) + Stilda
            u = tweak(U, reduction_degree)
            key = (order, Tuple(u))
            location = get(key_indices, key, nothing)
            if location === nothing
                group_index = get(group_indices, order, 0)
                if group_index == 0
                    push!(groups, (order, Any[]))
                    group_index = length(groups)
                    group_indices[order] = group_index
                end
                push!(groups[group_index][2], (key = (u, order), terms = [monomial]))
                key_indices[key] = (group_index, length(groups[group_index][2]))
            else
                group_index, entry_index = location
                push!(groups[group_index][2][entry_index].terms, monomial)
            end
        end
    end
    return groups
end

function _zero_reduction_vector(g::Vector{ZZRingElem})
    return [ZZ(0) for _ in eachindex(g)]
end

function _zero_reduction_vector(g)
    result = similar(g)
    my_zero!(result)
    return result
end

function _copy_reduction_vector(g)
    result = _zero_reduction_vector(g)
    my_copy!(result, g)
    return result
end

"""
    reduction_input_vector!(g, fe, key, S, cache, params)

Fill `g` with the summed contribution of the input terms belonging to `key`.
"""
function reduction_input_vector!(g, fe::ZetaFrontend, input, S, cache, params)
    P = cohomology_quotient(fe)
    reduction_degree = _reduction_g_degree(P)
    basis_lookup = cache[reduction_degree, :reverse]
    Stilda = zeros(Int, P.n + 1)
    for i in S
        Stilda[i+1] = 1
    end
    my_zero!(g)
    coefficients = Dict{Int,Any}()
    for monomial in input.terms
        U = exponent_vector(monomial, 1) + Stilda
        exponent = U - input.key[1]
        index = get(basis_lookup, exponent, nothing)
        index === nothing && continue
        coefficient = coeff(monomial, 1)
        coefficients[index] = get(coefficients, index, zero(coefficient)) + coefficient
    end
    for (index, coefficient) in coefficients
        CUDA.@allowscalar g[index] = lift(coefficient)
    end
    return g
end

"""
    reduction_terminal_polynomial(fe, chains, PR, S, params)

Convert terminal reduction chains to the polynomial output expected by the
final reduction stage. The chain iterator is consumed once and may be lazy.
"""
function reduction_terminal_polynomial(fe::ZetaFrontend, chains, PR, S, params)
    result = PR()
    for chain in chains
        result = reduction_add_terminal_chain(fe, result, chain, PR, params)
    end
    return _finish_reduction_terminal_polynomial(fe, result, PR, S)
end

"""
    reduction_add_terminal_chain(fe, result, chain, PR, params)

Add one finished terminal chain to a partially assembled reduced polynomial.
"""
function reduction_add_terminal_chain(fe::ZetaFrontend, result, chain, PR, params)
    P = cohomology_quotient(fe)
    reduction_degree = _reduction_g_degree(P)
    u, g_vec = chain
    cpu_g = Array(g_vec)
    CUDA.@allowscalar g =
        vector_to_polynomial(cpu_g, P.n, reduction_degree, PR, params.termorder)
    vars = gens(PR)
    return result + prod(vars .^ u) * g
end

function _finish_reduction_terminal_polynomial(fe::ZetaFrontend, result, PR, S)
    P = cohomology_quotient(fe)
    vars = gens(PR)
    XS = prod(PR(vars[i+1]) for i in S; init = PR(1))
    return [[div(result, XS), terminal_pole_order(P)]]
end
