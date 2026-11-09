@testset "cohomology presentations" begin
    for (n, d) in ((1, 3), (1, 4), (2, 4), (3, 3))
        presentation = DeRham.HypersurfacePresentation(n, d)
        @test presentation isa DeRham.CohomologyPresentation
        @test DeRham.pole_orders(presentation) == 1:n
        @test DeRham.terminal_pole_order(presentation) == n
        @test DeRham.numerator_degree(presentation, n) == d * n - (n + 1)
        for m in 1:n
            @test DeRham.numerator_degree(presentation, m) == d * m - (n + 1)
            @test DeRham.pole_exponent(presentation, m) == m
            @test DeRham.reduction_scalar(presentation, m) == Oscar.QQ(1, m)
        end
    end
end
