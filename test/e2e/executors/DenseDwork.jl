module DenseDworkExecutor

import ..DeRham
import ..Common

const INVARIANT_NAME = "dense_dwork_zeta"
const EXPECTED_SOURCE = :lpoly

function validate_parameters(name, parameters)
    isempty(parameters) || error("$name: dense Dwork accepts no parameters")
    return nothing
end

function expected(example, record)
    dimension = example.variety["dim"]
    Common.validate_geometry(
        example,
        INVARIANT_NAME;
        dimension = dimension,
        variable_count = dimension+2,
        degree = dimension == 2 ? 4 : 3,
    )
    return DeRham.zeta_function(
        Common.stored_lpoly_descending(example, record),
        Common.resolve_prime(example, record),
        dimension,
    )
end

execute(f, parameters) = DeRham.dense_dwork_zeta(f)
normalize(raw) = raw
compare(expected, actual) = expected == actual

end
