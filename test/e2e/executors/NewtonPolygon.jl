module NewtonPolygonExecutor

import ..DeRham
import ..ExampleCase
import ..ExecutionRecord
import ..Common
import ..ZetaCoefficientsExecutor

const INVARIANT_NAME = "newton_polygon"
const EXPECTED_SOURCE = :lpoly

validate_parameters(name::AbstractString, params::AbstractDict) =
    ZetaCoefficientsExecutor.validate_parameters(name, params)

function expected(example::ExampleCase, record::ExecutionRecord)
    coeffs = ZetaCoefficientsExecutor.expected(example, record)
    p = Common.resolve_prime(example, record)
    return DeRham.newton_polygon(coeffs, p)
end

function execute(F, params::AbstractDict)
    kwargs = Common.parse_kwargs(params)
    return DeRham.newton_polygon(F; kwargs...)
end

normalize(raw) = raw
compare(expected, actual) = expected == actual
compare(expected::DeRham.SlopesPolygon, actual::DeRham.SlopesPolygon) =
    expected.slopes == actual.slopes &&
    expected.slopelengths == actual.slopelengths &&
    expected.values == actual.values &&
    expected.slopesbefore == actual.slopesbefore

end
