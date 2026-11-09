@testset "graded exponent cache and presentation degrees" begin
    cache = DeRham.PolyExpCache(3, :invlex, [2], [2])
    @test !isempty(cache[2])
    @test isempty(cache[-1])
    @test cache[0] == [[0, 0, 0]]
    @test isempty(cache[-1, :reverse])
    @test cache[0, :reverse] == [[0, 0, 0]]
    missing_forward = try
        cache[1]
        nothing
    catch exception
        exception
    end
    @test missing_forward isa ArgumentError
    @test occursin("degree=1", sprint(showerror, missing_forward))
    @test_throws ArgumentError cache[1, :forward]
    missing_reverse = try
        cache[1, :reverse]
        nothing
    catch exception
        exception
    end
    @test missing_reverse isa ArgumentError
    @test occursin("degree=1", sprint(showerror, missing_reverse))
    @test_throws BoundsError cache[3]
    @test_throws BoundsError cache[3, :reverse]

    forward_cache = DeRham.controlled_reduction_cache(
        2,
        3,
        [0, 1, 2],
        DeRham.default_params(),
        DeRham.HypersurfacePresentation(2, 3),
    )
    @test DeRham._reduction_g_degree(DeRham.HypersurfacePresentation(2, 3)) == 4
    @test DeRham.cached_degrees(forward_cache) == [1, 2, 3, 4]
    @test DeRham.cached_reverse_lookups(forward_cache) == [3, 4]

    variable_cache = DeRham.controlled_reduction_cache(
        3,
        4,
        [3],
        DeRham.default_params(),
        DeRham.HypersurfacePresentation(3, 4),
    )
    @test DeRham._reduction_g_degree(DeRham.HypersurfacePresentation(3, 4)) == 9
    @test DeRham.cached_degrees(variable_cache) == [1, 4, 5, 8, 9, 12]
    @test DeRham.cached_reverse_lookups(variable_cache) == [4, 9, 12]

    R, (x, y, z) = polynomial_ring(GF(2), ["x", "y", "z"])
    f = x^2*y + y^2*z + z^2*x
    P = DeRham.HypersurfacePresentation(2, 3)
    basis = [[x, 1], [y, 2], [z, 2]]
    @test DeRham._hodge_polygon(basis, P).slopelengths == [2, 1]
    @test DeRham.get_basis_of_cohomology(DeRham.ProjectiveHypersurface(f)) ==
          DeRham.get_basis_of_cohomology(f)
end
