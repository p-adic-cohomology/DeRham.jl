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
const REFERENCE_TERM_ORDER = "invlex"
const REFERENCE_RELATIVE_PRECISION = 1
const PRIMITIVE_DIMENSION = 21

function validate_parameters(name::AbstractString, params::AbstractDict)
    Common.validate_common_parameters(name, params)
    algorithm = get(params, "algorithm", nothing)
    algorithm in SUPPORTED_ALGORITHMS || error(
        "execution record $(repr(name)): algorithm must be one of $(SUPPORTED_ALGORITHMS)",
    )
    get(params, "relative_precision", nothing) == REFERENCE_RELATIVE_PRECISION || error(
        "execution record $(repr(name)): relative_precision must be $(REFERENCE_RELATIVE_PRECISION)",
    )
    get(params, "termorder", REFERENCE_TERM_ORDER) == REFERENCE_TERM_ORDER || error(
        "execution record $(repr(name)): termorder must be $(repr(REFERENCE_TERM_ORDER))",
    )
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
    length(value) == PRIMITIVE_DIMENSION &&
    all(row -> row isa AbstractVector && length(row) == PRIMITIVE_DIMENSION, value) ||
        error(
            "example $(repr(example.id)) notes value for $(INVARIANT_NAME) at p=$(p) must be a 21 × 21 matrix",
        )
    return [[parse(BigInt, string(entry)) for entry in row] for row in value]
end

function execute(F, params::AbstractDict)
    kwargs = Common.parse_kwargs(params)
    delete!(kwargs, :relative_precision)
    kwargs[:termorder] = Symbol(REFERENCE_TERM_ORDER)
    basis = DeRham.get_basis_of_cohomology(F)
    matrix = DeRham.frobenius_matrix_with_precision(
        F,
        REFERENCE_RELATIVE_PRECISION;
        basis,
        kwargs...,
    )
    p = BigInt(Oscar.characteristic(parent(F)))
    n = Oscar.nvars(parent(F)) - 1
    column_moduli = [p^(REFERENCE_RELATIVE_PRECISION + n - term[2]) for term in basis]
    return (; matrix, column_moduli)
end

function normalize(raw)
    (; matrix, column_moduli) = raw
    rows = [
        [BigInt(Oscar.lift(Oscar.ZZ, matrix[i, j])) for j in axes(matrix, 2)] for
        i in axes(matrix, 1)
    ]
    return (; rows, column_moduli)
end

function compare(expected, actual)
    (; rows, column_moduli) = actual
    size_matches =
        length(column_moduli) == PRIMITIVE_DIMENSION &&
        length(rows) == PRIMITIVE_DIMENSION &&
        all(row -> length(row) == PRIMITIVE_DIMENSION, rows)
    return size_matches && all(
        iszero(mod(rows[i][j] - expected[i][j], column_moduli[j])) for
        i in eachindex(expected) for j in eachindex(expected[i])
    )
end

end
