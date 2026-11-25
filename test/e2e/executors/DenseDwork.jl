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
    Common.validate_geometry(
        example,
        INVARIANT_NAME;
        dimension = 2,
        variable_count = 4,
        degree = 4,
    )
    return DeRham.zeta_function(
        Common.stored_lpoly_descending(example, record),
        Common.resolve_prime(example, record),
        2,
    )
end

execute(f, parameters) = DeRham.dense_dwork_zeta(f)
normalize(raw) = raw
compare(expected, actual) = expected == actual

end
