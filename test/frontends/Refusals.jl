struct StubFrontend <: DeRham.ZetaFrontend end

@testset "Frontend algorithm refusals" begin
    for algorithm in (:pchunk, :varbyvar, :akr)
        @test_throws ArgumentError DeRham.zeta_coefficients(StubFrontend(); algorithm)
    end
    @test_throws ArgumentError DeRham.zeta_coefficients(StubFrontend(); use_gpu = true)
    @test_throws ArgumentError DeRham.zeta_coefficients(StubFrontend(); context = :custom)
    @test_throws ArgumentError DeRham.zeta_coefficients(StubFrontend(); series_precision = 3)
    @test_throws ArgumentError DeRham.zeta_coefficients(StubFrontend(); algorithm_precision = 3)
end
