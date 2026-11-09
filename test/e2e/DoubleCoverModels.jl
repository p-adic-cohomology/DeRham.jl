"""
Build finite-field branch polynomials from zeta_test_suite double-cover models.
"""
module DoubleCoverModels

using Oscar
import ..Common
import ..DeRham
import ..ExampleCase

function _invalid(example::ExampleCase, message::AbstractString)
    throw(ArgumentError("example $(repr(example.id)): $(message)"))
end

function _model(example::ExampleCase)
    model = get(example.variety, "model", nothing)
    model isa AbstractDict || _invalid(example, "variety.model must be an object")
    get(example.variety["coeff_domain"], "kind", nothing) == "integer" ||
        _invalid(example, "double-cover models require coeff_domain.kind == \"integer\"")
    return model
end

function _prime(example::ExampleCase, p::Integer)
    p isa Bool && _invalid(example, "p must be a non-boolean odd prime, got $(repr(p))")
    p > 2 || _invalid(example, "double-cover models require odd prime p; got p=$(p)")
    DeRham.is_prime(BigInt(p)) ||
        _invalid(example, "double-cover models require odd prime p; got p=$(p)")
    p <= typemax(Int) || _invalid(example, "p=$(p) is too large for a finite-field polynomial ring")
    return Int(p)
end

function _integer_coefficients(
    example::ExampleCase,
    model::AbstractDict,
    key::String;
    empty_is_zero::Bool = false,
)
    raw = get(model, key, nothing)
    raw isa AbstractVector && (empty_is_zero || !isempty(raw)) ||
        _invalid(example, "model.$(key) must be a nonempty array of integer coefficients")
    isempty(raw) && return BigInt[0]
    coefficients = BigInt[]
    for (index, value) in enumerate(raw)
        coefficient = try
            Common._to_bigint(value)
        catch
            _invalid(example, "model.$(key)[$(index)] must be an integer coefficient")
        end
        push!(coefficients, coefficient)
    end
    return coefficients
end

function _trim_trailing_zeros!(coefficients::Vector{BigInt})
    while length(coefficients) > 1 && iszero(last(coefficients))
        pop!(coefficients)
    end
    return coefficients
end

function _poly_degree(coefficients::Vector{BigInt})
    return length(coefficients) == 1 && iszero(only(coefficients)) ? -1 : length(coefficients) - 1
end

function _hyperelliptic(example::ExampleCase, model::AbstractDict, p::Int)
    get(example.variety, "dim", nothing) == 1 ||
        _invalid(example, "hyperelliptic models require variety.dim == 1")
    genus = get(example.variety, "genus", nothing)
    (genus isa Integer && !(genus isa Bool) && genus >= 0) ||
        _invalid(example, "hyperelliptic models require a nonnegative integer variety.genus")
    genus <= (typemax(Int) - 2) ÷ 2 || _invalid(example, "variety.genus is too large")
    target_degree = 2 * Int(genus) + 2

    f = _integer_coefficients(example, model, "f_coeffs_asc")
    h = _integer_coefficients(example, model, "h_coeffs_asc"; empty_is_zero = true)
    h_squared = zeros(BigInt, length(h) * 2 - 1)
    for i in eachindex(h), j in eachindex(h)
        h_squared[i+j-1] += h[i] * h[j]
    end

    coefficients = zeros(BigInt, max(length(f), length(h_squared)))
    for i in eachindex(f)
        coefficients[i] += 4 * f[i]
    end
    for i in eachindex(h_squared)
        coefficients[i] += h_squared[i]
    end
    _trim_trailing_zeros!(coefficients)

    branch_degree = _poly_degree(coefficients)
    branch_degree <= target_degree || _invalid(
        example,
        "degree(4f + h^2)=$(branch_degree) exceeds the degree $(target_degree) prescribed by genus $(genus)",
    )
    branch_degree >= target_degree - 1 || _invalid(
        example,
        "degree(4f + h^2)=$(branch_degree) is below degree $(target_degree - 1) required by genus $(genus)",
    )

    reduced = _trim_trailing_zeros!(mod.(coefficients, BigInt(p)))
    reduced_degree = _poly_degree(reduced)
    reduced_degree >= target_degree - 1 || _invalid(
        example,
        "reduction modulo p=$(p) lowers degree(4f + h^2) below $(target_degree - 1); " *
        "the homogeneous branch form is singular at infinity",
    )

    append!(reduced, zeros(BigInt, target_degree + 1 - length(reduced)))
    ring, _ = Oscar.polynomial_ring(Oscar.GF(p), [:x, :z])
    x, z = Oscar.gens(ring)
    branch = zero(ring)
    for degree in 0:target_degree
        coefficient = reduced[degree+1]
        iszero(coefficient) || (branch += ring(coefficient) * x^degree * z^(target_degree-degree))
    end
    return branch
end

function _variables(example::ExampleCase, raw, key::String, count::Int)
    raw isa AbstractVector && length(raw) == count ||
        _invalid(example, "model.$(key) must contain exactly $(count) variable names")
    all(name -> name isa AbstractString && !isempty(strip(name)), raw) ||
        _invalid(example, "model.$(key) must contain nonempty string variable names")
    names = String.(raw)
    length(unique(names)) == count ||
        _invalid(example, "model.$(key) must contain distinct variable names")
    return names
end

function _monomials(example::ExampleCase, raw, key::String, variable_count::Int)
    raw isa AbstractVector && !isempty(raw) ||
        _invalid(example, "model.$(key) must be a nonempty array of monomials")
    monomials = Vector{Vector{Int}}()
    for (index, monomial) in enumerate(raw)
        monomial isa AbstractVector && length(monomial) == variable_count ||
            _invalid(example, "model.$(key)[$(index)] must have $(variable_count) exponents")
        valid_exponents = all(
            exponent -> exponent isa Integer && !(exponent isa Bool) && 0 <= exponent <= typemax(Int),
            monomial,
        )
        valid_exponents ||
            _invalid(example, "model.$(key)[$(index)] must contain nonnegative integer exponents")
        push!(monomials, Int.(monomial))
    end
    return monomials
end

function _double_cover_P2(example::ExampleCase, model::AbstractDict, p::Int)
    get(example.variety, "dim", nothing) == 2 ||
        _invalid(example, "double_cover_P2 models require variety.dim == 2")
    variables = _variables(example, get(model, "branch_vars", nothing), "branch_vars", 3)
    monomials = _monomials(example, get(model, "branch_monomials", nothing), "branch_monomials", 3)
    coefficients = _integer_coefficients(example, model, "branch_coeffs")
    length(monomials) == length(coefficients) ||
        _invalid(example, "model.branch_monomials and model.branch_coeffs must have the same length")

    degrees = sum.(monomials)
    degree = first(degrees)
    degree > 0 && iseven(degree) ||
        _invalid(example, "double_cover_P2 branch polynomial must have positive even degree")
    all(==(degree), degrees) ||
        _invalid(example, "model.branch_monomials must all have the same total degree")

    ring, _ = Oscar.polynomial_ring(Oscar.GF(p), Symbol.(variables))
    return Common._multivariable_poly(monomials, coefficients, ring)
end

function build_input(example::ExampleCase, raw_prime::Integer)
    p = _prime(example, raw_prime)
    model = _model(example)
    kind = get(model, "kind", nothing)
    if kind == "hyperelliptic"
        return _hyperelliptic(example, model, p)
    elseif kind == "double_cover_P2"
        return _double_cover_P2(example, model, p)
    end
    _invalid(
        example,
        "unsupported double-cover model.kind $(repr(kind)); expected \"hyperelliptic\" or \"double_cover_P2\"",
    )
end

end
