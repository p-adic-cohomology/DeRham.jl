"""
Shared example, record, and capability helpers for E2E executors.
"""
module Common

using CUDA: CUDA
using JSON
using Oscar
import ..ExampleCase
import ..ExecutionRecord

_to_bigint(c::Integer) = BigInt(c)
_to_bigint(c::AbstractString) = parse(BigInt, c)

function _find_result(example::ExampleCase, p::Integer)
    matches = [r for r in example.results if _to_bigint(r["p"]) == BigInt(p)]
    isempty(matches) && error("example $(repr(example.id)) has no result for p=$(p)")
    length(matches) > 1 &&
        error("example $(repr(example.id)) has more than one result for p=$(p)")
    return matches[1]
end

function resolve_prime(example::ExampleCase, record::ExecutionRecord)
    if record.p !== nothing
        return record.p
    end
    length(example.results) == 1 || error(
        "execution record $(repr(record.name)) references example $(repr(example.id)), which has $(length(example.results)) results; 'p' is required to disambiguate",
    )
    return Int(_to_bigint(example.results[1]["p"]))
end

function stored_lpoly_descending(example::ExampleCase, record::ExecutionRecord)
    p = resolve_prime(example, record)
    result = _find_result(example, p)
    haskey(result, "Lpoly") || error(
        "example $(repr(example.id)), p=$(p): zeta_coefficients requires a Lpoly result, got keys $(collect(keys(result)))",
    )
    lpoly = result["Lpoly"]
    lpoly isa AbstractDict ||
        error("example $(repr(example.id)), p=$(p): Lpoly must be an object")
    haskey(lpoly, "coeffs_asc") ||
        error("example $(repr(example.id)), p=$(p): Lpoly is missing coeffs_asc")
    raw_coeffs = lpoly["coeffs_asc"]
    raw_coeffs isa AbstractVector ||
        error("example $(repr(example.id)), p=$(p): Lpoly.coeffs_asc must be an array")
    coeffs_asc = _to_bigint.(raw_coeffs)
    isempty(coeffs_asc) &&
        error("example $(repr(example.id)), p=$(p): Lpoly.coeffs_asc must be nonempty")
    if coeffs_asc[1] != 1
        if coeffs_asc[end] == 1
            error(
                "example $(repr(example.id)), p=$(p): coefficient-orientation mistake — " *
                "Lpoly.coeffs_asc looks reversed (descending); the schema requires ascending order with coeffs_asc[1] == 1",
            )
        else
            error(
                "example $(repr(example.id)), p=$(p): invalid Lpoly — coeffs_asc[1] must equal 1, got $(coeffs_asc[1])",
            )
        end
    end
    return reverse(coeffs_asc)
end

function validate_geometry(
    example::ExampleCase,
    invariant::AbstractString;
    model_kind = nothing,
    dimension = nothing,
    variable_count = nothing,
    degree = nothing,
)
    variety = example.variety
    model = variety["model"]
    kind = get(model, "kind", nothing)
    kind in ("plane_curve", "projective_hypersurface") || error(
        "example $(repr(example.id)): $(invariant) does not support model.kind $(repr(kind))",
    )
    if dimension !== nothing && get(variety, "dim", nothing) != dimension
        error(
            "example $(repr(example.id)): $(invariant) requires dimension $(dimension), got $(repr(get(variety, "dim", nothing)))",
        )
    end
    model_kind !== nothing &&
        kind != model_kind &&
        error(
            "example $(repr(example.id)): $(invariant) requires model.kind $(repr(model_kind)), got $(repr(kind))",
        )
    variables = get(model, "vars", nothing)
    if variable_count !== nothing &&
       (!(variables isa AbstractVector) || length(variables) != variable_count)
        error(
            "example $(repr(example.id)): $(invariant) requires $(variable_count) variables",
        )
    end
    if degree !== nothing
        monomials = get(model, "monomials", nothing)
        monomials isa AbstractVector ||
            error("example $(repr(example.id)): $(invariant) requires monomial model data")
        all(mono -> sum(Int, mono) == degree, monomials) || error(
            "example $(repr(example.id)): $(invariant) requires a homogeneous degree-$(degree) polynomial",
        )
    end
    return nothing
end

function _multivariable_poly(monomials, coeffs, R)
    length(monomials) == length(coeffs) ||
        error("monomials and coefficients must have the same length")
    gens_R = Oscar.gens(R)
    nv = length(gens_R)
    out = zero(R)
    for (mono, c) in zip(monomials, coeffs)
        length(mono) == nv ||
            error("monomial length $(length(mono)) does not match $(nv) variables")
        term = R(_to_bigint(c))
        for (j, e) in enumerate(mono)
            term *= gens_R[j]^Int(e)
        end
        out += term
    end
    return out
end

function build_input(example::ExampleCase, record::ExecutionRecord)
    p = resolve_prime(example, record)
    model = example.variety["model"]
    R, _ = Oscar.polynomial_ring(Oscar.GF(Int(p)), Symbol.(model["vars"]))
    return _multivariable_poly(model["monomials"], model["coeffs"], R)
end

function parse_kwargs(params::AbstractDict)
    kwargs = Dict{Symbol,Any}()
    for (k, v) in params
        sym = Symbol(k)
        kwargs[sym] = if sym in (:algorithm, :termorder)
            Symbol(v)
        elseif sym == :S
            Int.(v)
        else
            v
        end
    end
    return kwargs
end

function required_capabilities(record::ExecutionRecord)
    caps = Set{Symbol}()
    if "gpu" in record.tags || get(record.parameters, "use_gpu", false) == true
        push!(caps, :gpu)
    end
    return caps
end

capability_available(::Val{:gpu}) = CUDA.functional()
capability_available(cap::Symbol) = capability_available(Val(cap))

function parse_notes_invariants(notes::AbstractString)
    parsed_tag = nothing
    for (line_number, line) in enumerate(split(notes, '\n'))
        stripped = strip(line)
        startswith(stripped, "derham-invariants") || continue
        match_tag = match(r"^derham-invariants\s*:(.*)$", stripped)
        match_tag === nothing && error(
            "malformed derham-invariants tag on notes line $(line_number): expected 'derham-invariants: <JSON object>'",
        )
        parsed_tag === nothing ||
            error("malformed derham-invariants tag: only one tag line is allowed")
        payload = strip(match_tag.captures[1])
        isempty(payload) && error(
            "malformed derham-invariants tag on notes line $(line_number): JSON object is missing",
        )
        raw_tag = try
            JSON.parse(payload)
        catch e
            error(
                "malformed derham-invariants tag on notes line $(line_number): $(sprint(showerror, e))",
            )
        end
        raw_tag isa AbstractDict || error(
            "malformed derham-invariants tag on notes line $(line_number): expected a JSON object",
        )
        values = Dict{String,Dict{String,Any}}()
        for (prime_string, raw_values) in raw_tag
            prime = tryparse(BigInt, prime_string)
            (prime !== nothing && prime > 1 && string(prime) == prime_string) || error(
                "malformed derham-invariants tag on notes line $(line_number): prime key $(repr(prime_string)) must be a canonical prime string",
            )
            raw_values isa AbstractDict || error(
                "malformed derham-invariants tag on notes line $(line_number): values for prime $(prime_string) must be a JSON object",
            )
            values[prime_string] = Dict{String,Any}(String(k) => v for (k, v) in raw_values)
        end
        parsed_tag = values
    end
    return parsed_tag === nothing ? Dict{String,Dict{String,Any}}() : parsed_tag
end

function notes_invariant_value(example::ExampleCase, p::Integer, invariant::AbstractString)
    values = parse_notes_invariants(example.notes)
    prime_values = get(values, string(p), nothing)
    prime_values === nothing && error(
        "example $(repr(example.id)) notes have no derham-invariants values for p=$(p)",
    )
    haskey(prime_values, invariant) || error(
        "example $(repr(example.id)) notes have no value for invariant $(repr(invariant)) at p=$(p)",
    )
    return prime_values[invariant]
end

end
