using Test
using Oscar

# source files, in the future replace this with `using DeRham`
using DeRham
using CUDA

include("e2e/E2E.jl")
using .E2E

function run_de_rham_full_test_suite!()
    test_dir = @__DIR__
    return Base.eval(Main, quote
        include(joinpath($test_dir, "FirstEllipticCurveExample.jl"))
        include(joinpath($test_dir, "Precision.jl"))

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

        include(joinpath($test_dir, "invariants/HodgePolygon.jl"))
        include(joinpath($test_dir, "invariants/TateTwist.jl"))
        include(joinpath($test_dir, "invariants/PointCounts.jl"))
        include(joinpath($test_dir, "invariants/NewtonPolygon.jl"))
        include(joinpath($test_dir, "invariants/HasseWitt.jl"))
        include(joinpath($test_dir, "invariants/PicardRank.jl"))
        include(joinpath($test_dir, "invariants/ArtinMazurHeight.jl"))
        include(joinpath($test_dir, "invariants/ArtinInvariant.jl"))
        include(joinpath($test_dir, "invariants/Jacobian.jl"))

        include(joinpath($test_dir, "e2e/test/runtests.jl"))

        @testset "E2E" begin
            E2E.run!(ARGS)
        end

        include(joinpath($test_dir, "quality_gates.jl"))
    end)
end

E2E.dispatch_test_entrypoint(ARGS, E2E.run!, run_de_rham_full_test_suite!)
