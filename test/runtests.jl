const E2E_RECORD_SELECTOR_ONLY =
    !isempty(ARGS) &&
    all(arg -> startswith(arg, "--name=") || startswith(arg, "--tag="), ARGS)

using Test
using Oscar

# source files, in the future replace this with `using DeRham`
using DeRham
using CUDA

include("e2e/E2E.jl")
using .E2E

if E2E_RECORD_SELECTOR_ONLY
    E2E.run!(ARGS)
else
    include("FirstEllipticCurveExample.jl")
    include("Precision.jl")

    @testset "DeRham.jl" begin

        @testset "The curve y^2 - x^3 - x - 1 = 0, reproducing Costa's results" begin
            testEllCurve1_7() #TODO: why is this failing?
            testMonomialBasis()
            testLinAlgProb()
            testFrobTrans()
            testT()
        end

        @testset "Precision information" begin
            test_user_relative_precision_bypasses_default_table()
            test_max_auto_precision_without_default_table()
        end

    end

    include("invariants/HodgePolygon.jl")
    include("invariants/TateTwist.jl")
    include("invariants/PointCounts.jl")
    include("invariants/NewtonPolygon.jl")
    include("invariants/HasseWitt.jl")
    include("invariants/PicardRank.jl")
    include("invariants/ArtinMazurHeight.jl")
    include("invariants/ArtinInvariant.jl")
    include("invariants/Jacobian.jl")

    # A --name= or --tag=-only invocation runs just matching E2E records.
    # Bare Pkg.test() and --workflow= runs still execute conventional tests,
    # E2E harness tests, and quality gates.
    include("e2e/test/runtests.jl")

    @testset "E2E" begin
        E2E.run!(ARGS)
    end

    include("quality_gates.jl")
end
