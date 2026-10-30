"""
Executor for zeta-coefficient records.
"""
module ZetaCoefficientsExecutor

import ..DeRham
import ..ExampleCase
import ..ExecutionRecord
import ..Common

const SUPPORTED_MODEL_KINDS = Set(["plane_curve", "projective_hypersurface"])
const INVARIANT_NAME = "zeta_coefficients"
const EXPECTED_SOURCE = :lpoly

const KNOWN_PARAMETERS = Set([
    "S",
    "verbose",
    "changef",
    "givefrobmat",
    "algorithm",
    "termorder",
    "vars_reversed",
    "fastevaluation",
    "always_use_bigints",
    "use_gpu",
    "use_threads",
])
const KNOWN_ALGORITHMS = Set(["default", "depthfirst", "varbyvar", "pchunk"])

function validate_parameters(name::AbstractString, params::AbstractDict)
    for k in keys(params)
        k in KNOWN_PARAMETERS ||
            error("execution record $(repr(name)): unknown parameter $(repr(k))")
    end
    if haskey(params, "algorithm")
        params["algorithm"] isa AbstractString && params["algorithm"] in KNOWN_ALGORITHMS ||
            error(
                "execution record $(repr(name)): unknown algorithm $(repr(get(params, "algorithm", nothing)))",
            )
    end
    if haskey(params, "termorder") && !(params["termorder"] isa AbstractString)
        error("execution record $(repr(name)): termorder must be a string")
    end
    if haskey(params, "S") && !(params["S"] isa AbstractVector)
        error("execution record $(repr(name)): S must be an array of integers")
    end
    return nothing
end

function expected(example::ExampleCase, record::ExecutionRecord)
    variety = example.variety
    kind = variety["model"]["kind"]
    kind in SUPPORTED_MODEL_KINDS || error(
        "example $(repr(example.id)): zeta_coefficients does not support model.kind $(repr(kind))",
    )
    get(variety["coeff_domain"], "kind", nothing) == "integer" || error(
        "example $(repr(example.id)): zeta_coefficients only supports coeff_domain.kind == \"integer\"",
    )
    nmf = variety["non_middle_factors"]
    get(nmf, "kind", nothing) == "projective_lefschetz" ||
        error("example $(repr(example.id)): unsupported non_middle_factors $(nmf)")
    middle_factor_content = get(nmf, "middle_factor_content", nothing)
    dimension = get(variety, "dim", nothing)
    supported_content =
        middle_factor_content == "primitive" ||
        (dimension isa Integer && isodd(dimension) && middle_factor_content == "full")
    supported_content || error(
        "example $(repr(example.id)): DeRham.zeta_coefficients returns the primitive middle factor; got middle_factor_content=$(repr(middle_factor_content)) for dim=$(repr(dimension))",
    )

    return Common.stored_lpoly_descending(example, record)
end

function execute(F, params::AbstractDict)
    kwargs = Common.parse_kwargs(params)
    return DeRham.zeta_coefficients(F; kwargs...)
end

normalize(raw) = BigInt.(raw)

compare(expected::Vector{BigInt}, actual::Vector{BigInt}) = expected == actual

end # module ZetaCoefficientsExecutor
