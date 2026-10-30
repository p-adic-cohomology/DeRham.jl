module K3GeometricPicardRankExecutor

import ..DeRham
import ..ExampleCase
import ..ExecutionRecord
import ..Common
import ..ZetaCoefficientsExecutor

const INVARIANT_NAME = "k3_geometric_picard_rank"
const EXPECTED_SOURCE = :lpoly

validate_parameters(name::AbstractString, params::AbstractDict) =
    ZetaCoefficientsExecutor.validate_parameters(name, params)

function expected(example::ExampleCase, record::ExecutionRecord)
    Common.validate_geometry(
        example,
        INVARIANT_NAME;
        model_kind = "projective_hypersurface",
        dimension = 2,
        variable_count = 4,
        degree = 4,
    )
    coeffs = ZetaCoefficientsExecutor.expected(example, record)
    p = Common.resolve_prime(example, record)
    return DeRham.k3_geometric_picard_rank(coeffs, p)
end

function execute(F, params::AbstractDict)
    kwargs = Common.parse_kwargs(params)
    return DeRham.k3_geometric_picard_rank(F; kwargs...)
end

normalize(raw) = raw
compare(expected, actual) = expected == actual

end
