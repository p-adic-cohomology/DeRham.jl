"""
Executor for zeta-coefficient records.
"""
module ZetaCoefficientsExecutor

import ..DeRham
import ..ExampleCase
import ..ExecutionRecord
import ..Common

const SUPPORTED_MODEL_KINDS =
    Set(["plane_curve", "projective_hypersurface", "hyperelliptic", "double_cover_P2"])
const INVARIANT_NAME = "zeta_coefficients"
const EXPECTED_SOURCE = :lpoly

validate_parameters(name::AbstractString, params::AbstractDict) =
    Common.validate_common_parameters(name, params)

function _projective_precision_override(kwargs)
    has_algorithm_precision = haskey(kwargs, :algorithm_precision)
    has_other_precision =
        haskey(kwargs, :relative_precision) || haskey(kwargs, :series_precision)
    if has_algorithm_precision
        algorithm = get(kwargs, :algorithm, :default)
        algorithm in (:akr, :akr_naive, :akr_levelbylevel) || throw(
            ArgumentError(
                "projective zeta records can override algorithm_precision only for AKR",
            ),
        )
        has_other_precision && throw(
            ArgumentError(
                "projective zeta records with algorithm_precision use automatic relative and series precisions",
            ),
        )
        kwargs[:algorithm_precision] > 0 ||
            throw(ArgumentError("projective AKR algorithm_precision must be positive"))
    elseif has_other_precision
        throw(
            ArgumentError(
                "projective zeta records only accept the AKR algorithm_precision override",
            ),
        )
    end
    return has_algorithm_precision
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
    if kind in ("plane_curve", "projective_hypersurface")
        _projective_precision_override(Common.parse_kwargs(record.parameters))
    end

    return Common.stored_lpoly_descending(example, record)
end

function execute(F, params::AbstractDict)
    kwargs = Common.parse_kwargs(params)
    if F isa DeRham.DoubleCover || !_projective_precision_override(kwargs)
        return DeRham.zeta_coefficients(F; kwargs...)
    end

    algorithm_precision = pop!(kwargs, :algorithm_precision)::Int
    hodge_polygon, relative_precision, series_precision, _ =
        DeRham.calculate_precision_information(F)
    return Base.invokelatest(
        DeRham.zeta_coefficients_with_precision,
        F,
        hodge_polygon.slopelengths,
        relative_precision,
        series_precision,
        algorithm_precision;
        kwargs...,
    )
end

normalize(raw) = BigInt.(raw)

compare(expected::Vector{BigInt}, actual::Vector{BigInt}) = expected == actual

end # module ZetaCoefficientsExecutor
