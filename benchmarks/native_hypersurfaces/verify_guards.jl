using DeRham, Oscar, Test
include("fixtures.jl")
@testset "native hypersurface public refusals and frontends" begin
    for (n, p, outside) in ((4, 151, 257), (5, 13, 29), (6, 11, 173))
        case=only(filter(c->c.n==n && c.p==p, NATIVE_HYPERSURFACE_FIXTURES))
        f=native_fixture(case)
        xs=gens(parent(f))
        d=case.d
        @test_throws ArgumentError DeRham.dense_dwork_zeta(sum(xs)^d)
        @test_throws ArgumentError DeRham.dense_dwork_zeta(sum(x^d for x in xs))
        @test_throws ArgumentError DeRham.dense_dwork_zeta(sum(xs)^(d+1))
        R, ys=polynomial_ring(GF(outside), ["x$i" for i = 1:n])
        @test_throws ArgumentError DeRham.dense_dwork_zeta(sum(ys)^d)
        R, ys=polynomial_ring(GF(11, 2), ["x$i" for i = 1:n])
        @test_throws ArgumentError DeRham.dense_dwork_zeta(sum(ys)^d)
        @test_throws ErrorException DeRham._dense_dwork_matrix(f, p, n==4 ? 5 : 6)
        expected=DeRham.zeta_function(reverse(case.expected), p, n-2)
        @test DeRham.dense_dwork_zeta(DeRham.ProjectiveHypersurface(f))==expected
    end
end
