# E2E.jl — the declarative E2E catalogue. See test/e2e/README.md for the
# contributor workflow (adding examples, execution records, running by
# name/tag/workflow).

module E2E

using Test
using TOML
using JSON
using Printf
using Pkg.Artifacts
import Main.DeRham

include("Model.jl")
include("Config.jl")
include("Loader.jl")
include("executors/Common.jl")
include("executors/ZetaCoefficients.jl")
include("executors/NewtonPolygon.jl")
include("executors/PointCounts.jl")
include("executors/PicardRankBound.jl")
include("executors/K3GeometricPicardRank.jl")
include("executors/K3PicardRank.jl")
include("executors/K3PicardRealizationDegree.jl")
include("executors/ArtinMazurHeight.jl")

const EXECUTOR_REGISTRY = Dict{String,Module}(
    "zeta_coefficients" => ZetaCoefficientsExecutor,
    "newton_polygon" => NewtonPolygonExecutor,
    "point_counts" => PointCountsExecutor,
    "picard_rank_bound" => PicardRankBoundExecutor,
    "k3_geometric_picard_rank" => K3GeometricPicardRankExecutor,
    "k3_picard_rank" => K3PicardRankExecutor,
    "k3_picard_realization_degree" => K3PicardRealizationDegreeExecutor,
    "artin_mazur_height" => ArtinMazurHeightExecutor,
)
const INVARIANT_REGISTRY = let registry = Dict{String,Module}()
    for executor in values(EXECUTOR_REGISTRY)
        isdefined(executor, :INVARIANT_NAME) ||
            error("executor $(nameof(executor)) is missing INVARIANT_NAME")
        isdefined(executor, :EXPECTED_SOURCE) ||
            error("executor $(nameof(executor)) is missing EXPECTED_SOURCE")
        invariant = String(getproperty(executor, :INVARIANT_NAME))
        haskey(registry, invariant) &&
            error("duplicate executor registration for invariant $(repr(invariant))")
        getproperty(executor, :EXPECTED_SOURCE) in (:lpoly, :notes) ||
            error("executor $(nameof(executor)) has an unsupported EXPECTED_SOURCE")
        registry[invariant] = executor
    end
    registry
end

include("Validate.jl")
include("Select.jl")
include("Run.jl")

# Entry point called from test/runtests.jl with Base.ARGS (i.e. Pkg.test's
# test_args). Bare Pkg.test() passes no args, which defaults to --workflow=ci.
function run!(args::AbstractVector{<:AbstractString} = ARGS)
    return run_e2e(@__DIR__, args)
end

end # module E2E
