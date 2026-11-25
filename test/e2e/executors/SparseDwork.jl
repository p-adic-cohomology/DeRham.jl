module SparseDworkExecutor

import ..DeRham
import ..ZetaCoefficientsExecutor

const INVARIANT_NAME = "sparse_dwork_coefficients"
const EXPECTED_SOURCE = :lpoly

function validate_parameters(name::AbstractString, params::AbstractDict)
    isempty(params) || error("sparse Dwork record $name takes no parameters")
end

expected(example, record) = ZetaCoefficientsExecutor.expected(example, record)
execute(F, ::AbstractDict) = DeRham.sparse_dwork_coefficients(F)
normalize(raw) = BigInt.(raw)
compare(expected::Vector{BigInt}, actual::Vector{BigInt}) = expected == actual

end
