module PicardRankBoundExecutor

import ..DeRham
import ..ExampleCase
import ..ExecutionRecord
import ..Common
import ..ZetaCoefficientsExecutor

const INVARIANT_NAME = "picard_rank_bound"
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
    )
    coeffs = Common.stored_lpoly_descending(example, record)
    p = Common.resolve_prime(example, record)
    return DeRham.picard_rank_bound(coeffs, p)
end

function execute(F, params::AbstractDict)
    kwargs = Common.parse_kwargs(params)
    return DeRham.picard_rank_bound(F; kwargs...)
end

normalize(raw) = raw
compare(expected, actual) = expected == actual

end
