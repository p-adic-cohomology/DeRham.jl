"""
Executor for tests that compare experimental AKR Frobenius matrices with the
default reduction at a requested relative precision.
"""
module AKRFrobeniusMatrixAgreementExecutor

import ..DeRham
import ..ExampleCase
import ..ExecutionRecord
import ..Common
import Oscar

const INVARIANT_NAME = "akr_frobenius_matrix"
const EXPECTED_SOURCE = :notes
const SUPPORTED_ALGORITHMS = ("akr", "akr_naive", "akr_levelbylevel")

function validate_parameters(name::AbstractString, params::AbstractDict)
    Common.validate_common_parameters(name, params)
    algorithm = get(params, "algorithm", nothing)
    algorithm in SUPPORTED_ALGORITHMS || error(
        "execution record $(repr(name)): algorithm must be one of $(SUPPORTED_ALGORITHMS)",
    )
    get(params, "relative_precision", nothing) == 1 ||
        error("execution record $(repr(name)): relative_precision must be 1")
    get(params, "changef", nothing) === false ||
        error("execution record $(repr(name)): changef must be false")
    for key in ("algorithm_precision", "series_precision", "givefrobmat")
        haskey(params, key) &&
            error("execution record $(repr(name)): parameter $(repr(key)) is unsupported")
    end
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
    p = Common.resolve_prime(example, record)
    value = Common.notes_invariant_value(example, p, INVARIANT_NAME)
    value isa AbstractVector || error(
        "example $(repr(example.id)) notes value for $(INVARIANT_NAME) at p=$(p) must be a matrix",
    )
    dimension = length(value)
    dimension > 0 &&
    all(row -> row isa AbstractVector && length(row) == dimension, value) || error(
        "example $(repr(example.id)) notes value for $(INVARIANT_NAME) at p=$(p) must be a nonempty square matrix",
    )
    return [[parse(BigInt, string(entry)) for entry in row] for row in value]
end

function execute(F, params::AbstractDict)
    kwargs = Common.parse_kwargs(params)
    precision = pop!(kwargs, :relative_precision)
    algorithm = pop!(kwargs, :algorithm)
    return DeRham.frobenius_matrix_with_precision(
        F,
        precision;
        algorithm = algorithm,
        kwargs...,
    )
end

function normalize(raw)
    return [
        [parse(BigInt, string(Oscar.lift(Oscar.ZZ, raw[i, j]))) for j = 1:size(raw, 2)] for
        i = 1:size(raw, 1)
    ]
end

compare(expected::Vector{Vector{BigInt}}, actual::Vector{Vector{BigInt}}) =
    expected == actual

end
