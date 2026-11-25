# julia --project=. benchmarks/dense_dwork/verify.jl
using DeRham, Oscar, Test
include("fixtures.jl")

function small_dense(p)
    R, xs = polynomial_ring(GF(p), ["x", "y", "z", "w"])
    # Coefficients are deterministic; geometric smoothness checked by both solvers.
    S = [[a, b, c, 4-a-b-c] for a = 0:4 for b = 0:(4-a) for c = 0:(4-a-b)]
    f = sum(R(mod(i^3+3i+1, p-1)+1)*prod(xs .^ e) for (i, e) in enumerate(S))
    return f
end

function direct_projective_count(f, p)
    K = coefficient_ring(parent(f))
    count = 0
    for first = 1:4
        tails = Iterators.product(ntuple(_ -> 0:(p-1), 4-first)...)
        for tail in tails
            point = vcat(zeros(Int, first-1), [1], collect(tail))
            count += iszero(evaluate(f, K.(point)))
        end
    end
    return count
end

@testset "four archived primitive polynomials and ambient factors" begin
    for case in filter(c -> c.p >= 37, DENSE_DWORK_FIXTURES)
        f = dense_fixture(case)
        cs = DeRham.dense_dwork_coefficients(f)
        @test reverse(BigInt.(cs)) == case.expected
        Z = DeRham.dense_dwork_zeta(f)
        R, T = polynomial_ring(ZZ, "T")
        unit = coeff(denominator(Z), 0)
        @test numerator(Z) == unit
        @test denominator(Z) == unit*R(case.expected)*(1-T)*(1-case.p*T)*(1-case.p^2*T)
        @test DeRham.dense_dwork_zeta(DeRham.ProjectiveHypersurface(f)) == Z
    end
end

@testset "small primes: direct counts and independent controlled reduction" begin
    for p in (11, 13)
        f = small_dense(p)
        @test length(f) == 35
        cs = DeRham.dense_dwork_coefficients(f)
        points = direct_projective_count(f, p)
        @test points == 1+p+p^2-BigInt(cs[end-1])
        independent =
            DeRham.zeta_coefficients(f; use_gpu = false, use_threads = false, verbose = 0)
        @test cs == independent
        println("small_p=", p, " points=", points, " coefficients=", BigInt.(reverse(cs)))
        flush(stdout)
    end
end

@testset "domain and singularity" begin
    R, xs = polynomial_ring(GF(11), ["x", "y", "z", "w"])
    singular = sum(xs)^4
    @test length(singular) == 35
    @test_throws ArgumentError DeRham.dense_dwork_zeta(singular)
    @test_throws ArgumentError DeRham.dense_dwork_zeta(sum(x^4 for x in xs))
    for p in (7, 257)
        @test_throws ArgumentError DeRham.dense_dwork_zeta(small_dense(p))
    end
    Rq, ys = polynomial_ring(GF(11, 2), ["x", "y", "z", "w"])
    @test_throws ArgumentError DeRham.dense_dwork_zeta(sum(ys)^4)
end

@testset "maximum native precision" begin
    for p in (11, 13)
        f = small_dense(p)
        A = DeRham._dense_dwork_matrix(f, p, 5)
        cs = DeRham._dense_dwork_lift(DeRham._dense_dwork_residues(A, p, 5), p)
        @test cs !== nothing
        @test ZZ.(reverse(cs)) == DeRham.dense_dwork_coefficients(f)
    end
    f = small_dense(47)
    A = DeRham._dense_dwork_matrix(f, 47, 5)
    @test length(A) == 21
    cs = DeRham._dense_dwork_lift(DeRham._dense_dwork_residues(A, 47, 5), 47)
    @test cs !== nothing
    @test ZZ.(reverse(cs)) == DeRham.dense_dwork_coefficients(f)
    @test_throws ErrorException DeRham._dense_dwork_matrix(f, 47, 6)
    @test_throws ErrorException DeRham._dense_dwork_matrix(small_dense(53), 53, 5)
end
