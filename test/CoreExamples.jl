include("FirstEllipticCurveExample.jl")
include("Precision.jl")

@testset "DeRham.jl" begin

    @testset "The curve y^2 - x^3 - x - 1 = 0, reproducing Costa's results" begin
        testEllCurve1_7()
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
