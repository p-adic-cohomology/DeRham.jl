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

@testset "dense hypersurface precision ranges and functional equations" begin
    for (n, largest) in ((4, 251), (5, 23), (6, 167)), p = 11:largest
        DeRham.is_prime(p) || continue
        for N in (n==4 ? (3:(p<=47 ? 5 : 4)) : (5:5))
            h=n==4 ? 1 : 2
            bound=((p*N-h)*(p-1))//p^2+(h-N)//(p-1)
            @test bound>N-1
            @test big(p)^(N-(n==6))<big(2)^32
            @test big(p)^(2N-3)<big(2)^52
        end
    end
    R, T=DeRham.polynomial_ring(DeRham.ZZ, "T")
    for (dimension, P) in (
        (2, (1-11T)^21),
        (3, (1+11^3*T^2)^5),
        (4, (1-11^2*T)^22),
        (4, (1-11^2*T)^21*(1+11^2*T)),
    )
        c=BigInt[DeRham.coeff(P, i) for i = 0:DeRham.Oscar.degree(P)]
        modulus=big(11)^100
        data=[(mod(x, modulus), modulus) for x in c]
        @test DeRham._dense_dwork_lift(data, 11; dimension)==c
    end
end
