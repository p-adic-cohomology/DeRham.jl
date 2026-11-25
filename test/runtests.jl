const TEST_ENTRY_STARTED_AT = time_ns()

using Test
using Oscar

# source files, in the future replace this with `using DeRham`
using DeRham
using CUDA

include("e2e/E2E.jl")
using .E2E

println(
    "[tests] load: wall=$(round((time_ns() - TEST_ENTRY_STARTED_AT) / 1.0e9, digits = 3))s " *
    "Julia=$(VERSION) source=$(pathof(DeRham)) project=$(Base.active_project())",
)
flush(stdout)

"""
Unit test files (paths relative to `test/`) and the workflow tier each runs
in. `presubmit` files run in every workflow; `ci` files run in `ci` and
`full`. The quality gates run in `ci` and `full`.
"""
const UNIT_TEST_TIERS = [
    "CoreExamples.jl" => :ci,
    "invariants/HodgePolygon.jl" => :presubmit,
    "invariants/TateTwist.jl" => :presubmit,
    "invariants/PointCounts.jl" => :presubmit,
    "invariants/NewtonPolygon.jl" => :ci,
    "invariants/HasseWitt.jl" => :ci,
    "invariants/PicardRank.jl" => :presubmit,
    "invariants/ArtinMazurHeight.jl" => :ci,
    "invariants/ArtinInvariant.jl" => :ci,
    "invariants/CharPolyFrob.jl" => :presubmit,
    "invariants/Jacobian.jl" => :ci,
    "quotients/CohomologyQuotient.jl" => :presubmit,
    "frontends/ProjectiveHypersurface.jl" => :presubmit,
    "frontends/DoubleCover.jl" => :ci,
    "frontends/Refusals.jl" => :presubmit,
    "reduction/Selection.jl" => :presubmit,
    "reduction/AKR.jl" => :ci,
    "reduction/FrontendMethods.jl" => :presubmit,
    "reduction/EngineInterfaces.jl" => :presubmit,
    "reduction/AccountingEngine.jl" => :presubmit,
    "reduction/MatrixEngine.jl" => :presubmit,
    "reduction/RuvRebinding.jl" => :presubmit,
    "reduction/DepthFirstPolicy.jl" => :presubmit,
    "reduction/VarByVarPolicy.jl" => :presubmit,
    "reduction/PChunkDirection.jl" => :presubmit,
    "reduction/PChunkLimits.jl" => :presubmit,
    "reduction/PChunkThreading.jl" => :ci,
    "graded_ring/GradedExpCache.jl" => :presubmit,
    "e2e/test/runtests.jl" => :ci,
]

function split_unit_test_args(args::AbstractVector{<:AbstractString})
    units = String[]
    rest = String[]
    for arg in args
        if startswith(arg, "--unit=")
            push!(units, arg[(length("--unit=")+1):end])
        else
            push!(rest, arg)
        end
    end
    return units, rest
end

function include_unit_test(file::AbstractString)
    path = joinpath(@__DIR__, file)
    isfile(path) || error(
        "unknown unit test file $(repr(file)); pass a path relative to test/, " *
        "for example --unit=reduction/AKR.jl",
    )
    Base.include(Main, path)
    return nothing
end

function unit_test_workflow(args::AbstractVector{<:AbstractString})
    selection = E2E.parse_selection(args)
    return selection.mode == :workflow ? selection.value : :ci
end

function run_de_rham_test_suite!(args::AbstractVector{<:AbstractString})
    workflow = unit_test_workflow(args)
    primed = String[]
    priming_timing = @timed @testset "E2E priming" begin
        append!(primed, E2E.prime!(args))
    end
    println(
        "[tests] priming: wall=$(round(priming_timing.time, digits = 3))s " *
        "compile=$(round(priming_timing.compile_time, digits = 3))s",
    )
    flush(stdout)

    started_at = time_ns()
    unit_timing = @timed for (file, tier) in UNIT_TEST_TIERS
        if tier == :presubmit || workflow in (:ci, :full)
            timing = @timed include_unit_test(file)
            println(
                "[tests] unit $(file): wall=$(round(timing.time, digits = 3))s " *
                "compile=$(round(timing.compile_time, digits = 3))s",
            )
            flush(stdout)
        end
    end

    e2e_timing = @timed @testset "E2E" begin
        E2E.run!(args; primed = primed)
    end

    elapsed = (time_ns() - started_at) / 1.0e9
    println(
        "[tests] after priming: wall=$(round(elapsed, digits = 3))s " *
        "units=$(round(unit_timing.time, digits = 3))s " *
        "e2e=$(round(e2e_timing.time, digits = 3))s " *
        "compile=$(round(unit_timing.compile_time + e2e_timing.compile_time, digits = 3))s",
    )
    flush(stdout)

    workflow in (:ci, :full) && include_unit_test("quality_gates.jl")
    return nothing
end

function run_unit_test_selection!(
    units::AbstractVector{<:AbstractString},
    rest::AbstractVector{<:AbstractString},
)
    isempty(rest) ||
        E2E.record_selector_only(rest) ||
        error(
            "--unit= can be combined only with --name= and --tag=; got $(join(rest, " "))",
        )
    foreach(include_unit_test, units)
    isempty(rest) || @testset "E2E" begin
        E2E.run!(rest)
    end
    return nothing
end

let (units, rest) = split_unit_test_args(ARGS)
    if isempty(units)
        E2E.dispatch_test_entrypoint(rest, E2E.run!, () -> run_de_rham_test_suite!(rest))
    else
        run_unit_test_selection!(units, rest)
    end
end
