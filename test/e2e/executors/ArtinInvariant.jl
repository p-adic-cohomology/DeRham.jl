"""
Executor for Artin-invariant records whose expected values are stored in notes.
"""
module ArtinInvariantExecutor

import ..DeRham
import ..ExampleCase
import ..ExecutionRecord
import ..Common

const INVARIANT_NAME = "artin_invariant"
const EXPECTED_SOURCE = :notes

function validate_parameters(name::AbstractString, params::AbstractDict)
    isempty(params) || error("execution record $(repr(name)): artin_invariant accepts no parameters")
    return nothing
end

function expected(example::ExampleCase, record::ExecutionRecord)
    Common.validate_geometry(
        example,
        INVARIANT_NAME;
        model_kind = "projective_hypersurface",
        dimension = 2,
        variable_count = 4,
        degree = 4,
    )
    get(example.variety["coeff_domain"], "kind", nothing) == "integer" || error(
        "example $(repr(example.id)): artin_invariant only supports coeff_domain.kind == \"integer\"",
    )

    p = Common.resolve_prime(example, record)
    p > 2 && DeRham.is_prime(p) || error(
        "example $(repr(example.id)): artin_invariant requires an odd prime field, got p=$(p)",
    )
    coeffs = Common.stored_lpoly_descending(example, record)
    height = DeRham.artin_mazur_height(coeffs, p)
    height == Inf || error(
        "example $(repr(example.id)): artin_invariant requires a supersingular quartic K3; stored Lpoly gives finite Artin-Mazur height $(height)",
    )

    value = Common.notes_invariant_value(example, p, INVARIANT_NAME)
    value isa Integer && !(value isa Bool) || error(
        "example $(repr(example.id)) notes value for artin_invariant at p=$(p) must be an integer",
    )
    1 <= value <= 10 || error(
        "example $(repr(example.id)) notes value for artin_invariant at p=$(p) must be in 1:10",
    )
    return Int(value)
end

function execute(f, params::AbstractDict)
    isempty(params) || error("artin_invariant accepts no execution parameters")
    return DeRham.artin_invariant(f)
end

function normalize(raw)
    raw isa Integer && !(raw isa Bool) || error("artin_invariant returned a non-integer value")
    return Int(raw)
end

compare(expected::Integer, actual::Integer) = expected == actual

end
