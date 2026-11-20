@testset "pchunk direction selection" begin
    R, variables = polynomial_ring(GF(5), ["x$i" for i = 0:4])
    fe = DeRham.ProjectiveHypersurface(sum(variable^3 for variable in variables))
    S = [0, 1, 2]
    directions = DeRham.ReductionDirectionSet(DeRham.reduction_directions(fe, S))
    key = DeRham.ReductionChainKey([11, 8, 13, 9, 6], 7)

    direction = DeRham._pchunk_choose_direction(key, 3, S)
    @test direction == [1, 1, 1, 0, 0]
    chunk = DeRham._depthfirst_chunk_length(collect(key.u), direction, 3)
    @test chunk == 3
    @test directions(key.u, chunk, direction)
    @test any(!iszero(rem(exponent, 2)) for exponent in key.u)

    invalid_key = DeRham.ReductionChainKey([0, 0, 0, 0, 0], 7)
    @test_throws ArgumentError DeRham._pchunk_choose_direction(invalid_key, 3, S)
end
