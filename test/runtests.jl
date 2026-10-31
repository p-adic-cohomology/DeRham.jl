using Test
using Oscar

# source files, in the future replace this with `using DeRham`
using DeRham
using CUDA

include("e2e/E2E.jl")
using .E2E

const E2E_TEST_ENTRYPOINT_ROUTE = E2E.test_entrypoint_route(ARGS)

if E2E_TEST_ENTRYPOINT_ROUTE == :record_only
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

    include("e2e/test/runtests.jl")

    @testset "E2E" begin
        E2E.run!(ARGS)
    end

    include("quality_gates.jl")
end
