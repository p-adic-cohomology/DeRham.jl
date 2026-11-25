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

function algorithm_precision_kwargs(kwargs)
    return (
        S = Int.(get(kwargs, :S, [-1])),
        verbose = Int(get(kwargs, :verbose, 0)),
        changef = Bool(get(kwargs, :changef, true)),
        givefrobmat = Bool(get(kwargs, :givefrobmat, false)),
        algorithm = Symbol(get(kwargs, :algorithm, :default)),
        termorder = Symbol(get(kwargs, :termorder, :invlex)),
        vars_reversed = Bool(get(kwargs, :vars_reversed, false)),
        fastevaluation = Bool(get(kwargs, :fastevaluation, true)),
        always_use_bigints = Bool(get(kwargs, :always_use_bigints, false)),
        use_gpu = Bool(get(kwargs, :use_gpu, false)),
        use_threads = Bool(get(kwargs, :use_threads, false)),
    )
end

function execute(F, params::AbstractDict)
    kwargs = Common.parse_kwargs(params)
    if haskey(kwargs, :algorithm_precision) && !(F isa DeRham.DoubleCover)
        algorithm_precision = pop!(kwargs, :algorithm_precision)::Int
        if haskey(kwargs, :relative_precision) || haskey(kwargs, :series_precision)
            throw(
                ArgumentError(
                    "projective zeta records with algorithm_precision use automatic relative and series precisions",
                ),
            )
        end
        hodge_polygon, relative_precision, series_precision, _ =
            DeRham.calculate_precision_information(F)
        if isempty(kwargs)
            return Base.invokelatest(
                DeRham.zeta_coefficients_with_precision,
                F,
                hodge_polygon.slopelengths,
                relative_precision,
                series_precision,
                algorithm_precision;
                algorithm = :default,
            )
        elseif length(kwargs) == 1 && haskey(kwargs, :algorithm)
            algorithm = kwargs[:algorithm]::Symbol
            return Base.invokelatest(
                DeRham.zeta_coefficients_with_precision,
                F,
                hodge_polygon.slopelengths,
                relative_precision,
                series_precision,
                algorithm_precision;
                algorithm = algorithm,
            )
        end
        precision_kwargs = algorithm_precision_kwargs(kwargs)
        return Base.invokelatest(
            DeRham.zeta_coefficients_with_precision,
            F,
            hodge_polygon.slopelengths,
            relative_precision,
            series_precision,
            algorithm_precision;
            precision_kwargs...,
        )
    end
    return DeRham.zeta_coefficients(F; kwargs...)
end

normalize(raw) = BigInt.(raw)

compare(expected::Vector{BigInt}, actual::Vector{BigInt}) = expected == actual

end # module ZetaCoefficientsExecutor
