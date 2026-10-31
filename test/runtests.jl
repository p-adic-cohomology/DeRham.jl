using Test
using Oscar

# source files, in the future replace this with `using DeRham`
using DeRham
using CUDA

include("e2e/E2E.jl")
using .E2E

function run_de_rham_full_test_suite!()
    include(joinpath(@__DIR__, "FirstEllipticCurveExample.jl"))
    include(joinpath(@__DIR__, "Precision.jl"))

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

    include(joinpath(@__DIR__, "invariants/HodgePolygon.jl"))
    include(joinpath(@__DIR__, "invariants/TateTwist.jl"))
    include(joinpath(@__DIR__, "invariants/PointCounts.jl"))
    include(joinpath(@__DIR__, "invariants/NewtonPolygon.jl"))
    include(joinpath(@__DIR__, "invariants/HasseWitt.jl"))
    include(joinpath(@__DIR__, "invariants/PicardRank.jl"))
    include(joinpath(@__DIR__, "invariants/ArtinMazurHeight.jl"))
    include(joinpath(@__DIR__, "invariants/ArtinInvariant.jl"))
    include(joinpath(@__DIR__, "invariants/Jacobian.jl"))

    include(joinpath(@__DIR__, "e2e/test/runtests.jl"))

    @testset "E2E" begin
        E2E.run!(ARGS)
    end

    include(joinpath(@__DIR__, "quality_gates.jl"))
end

E2E.dispatch_test_entrypoint(ARGS, E2E.run!, run_de_rham_full_test_suite!)
