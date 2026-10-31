module PointCountsExecutor

import ..DeRham
import ..ExampleCase
import ..ExecutionRecord
import ..Common

const INVARIANT_NAME = "point_counts"
const EXPECTED_SOURCE = :lpoly

function extension_degree(name::AbstractString, params::AbstractDict)
    haskey(params, "k") ||
        error("execution record $(repr(name)): point_counts requires parameter 'k'")
    k = params["k"]
    (k isa Integer && !(k isa Bool) && k > 0) ||
        error("execution record $(repr(name)): parameter 'k' must be a positive integer")
    return Int(k)
end

function zeta_parameters(params::AbstractDict)
    return Dict{String,Any}(key => value for (key, value) in params if key != "k")
end

function validate_parameters(name::AbstractString, params::AbstractDict)
    extension_degree(name, params)
    Common.validate_common_parameters(name, zeta_parameters(params))
    return nothing
end

function expected(example::ExampleCase, record::ExecutionRecord)
    coeffs = Common.stored_lpoly_descending(example, record)
    p = Common.resolve_prime(example, record)
    m = Int(example.variety["dim"])
    k = extension_degree(record.name, record.parameters)
    return BigInt.(DeRham.point_counts(coeffs, p, m, k))
end

function execute(F, params::AbstractDict)
    k = extension_degree("record", params)
    kwargs = Common.parse_kwargs(zeta_parameters(params))
    return DeRham.point_counts(F, k; kwargs...)
end

normalize(raw) = raw === false ? false : BigInt.(raw)
compare(expected, actual) = expected == actual

end
