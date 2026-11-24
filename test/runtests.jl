using Test
using Oscar

# source files, in the future replace this with `using DeRham`
using DeRham
using CUDA

include("e2e/E2E.jl")
using .E2E

function run_de_rham_full_test_suite!()
    test_dir = @__DIR__
    return Base.eval(
        Main,
        quote
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
                    test_explicit_precision_information()
                    test_max_auto_precision_without_default_table()
                    test_certified_series_precision()
                end

                @testset "AKR lattice bound" begin
                    test_akr_lattice_bound()
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
            include(joinpath($test_dir, "invariants/CharPolyFrob.jl"))
            include(joinpath($test_dir, "invariants/Jacobian.jl"))
            include(joinpath($test_dir, "invariants/QuadricPencil.jl"))
            include(joinpath($test_dir, "quotients/CohomologyQuotient.jl"))
            include(joinpath($test_dir, "frontends/ProjectiveHypersurface.jl"))
            include(joinpath($test_dir, "frontends/DoubleCover.jl"))
            include(joinpath($test_dir, "frontends/Refusals.jl"))
            include(joinpath($test_dir, "reduction/Selection.jl"))
            include(joinpath($test_dir, "reduction/AKR.jl"))
            include(joinpath($test_dir, "reduction/FrontendMethods.jl"))
            include(joinpath($test_dir, "reduction/EngineInterfaces.jl"))
            include(joinpath($test_dir, "reduction/AccountingEngine.jl"))
            include(joinpath($test_dir, "reduction/MatrixEngine.jl"))
            include(joinpath($test_dir, "reduction/RuvRebinding.jl"))
            include(joinpath($test_dir, "reduction/DepthFirstPolicy.jl"))
            include(joinpath($test_dir, "reduction/VarByVarPolicy.jl"))
            include(joinpath($test_dir, "reduction/PChunkDirection.jl"))
            include(joinpath($test_dir, "reduction/PChunkLimits.jl"))
            include(joinpath($test_dir, "reduction/PChunkThreading.jl"))
            include(joinpath($test_dir, "graded_ring/GradedExpCache.jl"))

            include(joinpath($test_dir, "e2e/test/runtests.jl"))

            @testset "E2E" begin
                E2E.run!(ARGS)
            end

            include(joinpath($test_dir, "quality_gates.jl"))
        end,
    )
end

E2E.dispatch_test_entrypoint(ARGS, E2E.run!, run_de_rham_full_test_suite!)
