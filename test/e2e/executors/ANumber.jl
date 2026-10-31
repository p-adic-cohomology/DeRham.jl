"""
Executor for a-number records whose expected values are stored in notes.
"""
module ANumberExecutor

import ..DeRham
import ..ExampleCase
import ..ExecutionRecord
import ..Common

const INVARIANT_NAME = "a_number"
const EXPECTED_SOURCE = :notes

function validate_parameters(name::AbstractString, params::AbstractDict)
    for parameter in keys(params)
        error("execution record $(repr(name)): unknown parameter $(repr(parameter))")
    end
    return nothing
end

function expected(example::ExampleCase, record::ExecutionRecord)
    Common.validate_geometry(
        example,
        INVARIANT_NAME;
        model_kind = "plane_curve",
        dimension = 1,
        variable_count = 3,
    )
    get(example.variety["coeff_domain"], "kind", nothing) == "integer" || error(
        "example $(repr(example.id)): a_number only supports coeff_domain.kind == \"integer\"",
    )
    p = Common.resolve_prime(example, record)
    value = Common.notes_invariant_value(example, p, INVARIANT_NAME)
    value isa Integer && !(value isa Bool) || error(
        "example $(repr(example.id)) notes value for a_number at p=$(p) must be an integer",
    )
    value >= 0 || error(
        "example $(repr(example.id)) notes value for a_number at p=$(p) must be nonnegative",
    )
    return Int(value)
end

function execute(f, params::AbstractDict)
    isempty(params) || error("a_number accepts no execution parameters")
    return DeRham.a_number(f)
end

function normalize(raw)
    raw isa Integer && !(raw isa Bool) || error("a_number returned a non-integer value")
    return Int(raw)
end

compare(expected::Integer, actual::Integer) = expected == actual

end # module ANumberExecutor
