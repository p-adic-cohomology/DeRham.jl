using Test
using DeRham

@testset "dense Dwork reconstruction guards" begin
    # No congruence information cannot produce a certified answer.
    @test DeRham._dense_dwork_lift(fill((big(0), big(1)), 22), 11) === nothing
    # A determinant known to be zero contradicts determinant = +/- p^21.
    residues = vcat([(big(1), big(1))], fill((big(0), big(11)^100), 21))
    @test_throws ErrorException DeRham._dense_dwork_lift(residues, 11)
    # Exact rational cutoff test for the complete accepted parameter set.
    for p = 11:101
        DeRham.is_prime(p) || continue
        for N = 3:(p<=47 ? 5 : 4)
            bound = ((p*N-1)*(p-1))//p^2 + (1-N)//(p-1)
            @test bound > N-1
            @test big(p)^N < big(2)^28
            @test big(p)^(2N-3) < big(2)^39
        end
    end
end
